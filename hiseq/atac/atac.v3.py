"""
ATAC-seq pipeline

1. cutadapt, recursive trimming 

2. Alignment
--very-sensitive -X 2000

3. Mito%

pysam.idxstats(), chrM

4. Remove PCR duplicates

$ sambamba markdup

5. unique alignments

$ samtools view -q 30


6. proper mapped

$ samtools view -F 1804 

0x4 Read unmapped
0x8 Mate unmapped 
0x100 Not primary alignment 
0x200 Read fails platform/vendor quality checks
0x400 Read is PCR or optical duplicate

7. merge replicate

$ samtools merge 

8. Shifting reads (For TF footprinting)
For Tn5 insertion: +4, -5

$ samtools sort -n 
$ bedtools bamtobed -i in.bam -bedpe | awk -v OFS="\t" '{($9=="+"){print $1,$2+4,$6+4} \
  ($9=="-"){print $1,$2-5,$6-5}}' > fragments.bed   


9. Peak calling (MACS2)

$ macs2 callpeak --nomodel -f BAMPE --keep-dup all --cutoff-analysis 
$ macs2 callpeak 

10. bigWig

$ bamCoverage --binSize 10 --normalizeUsing RPGC 


# use --ATACshift
alignmentSieve --numberOfProcessors 8 --ATACshift --bam sample1.bam -o sample1.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ 8 -O bam -o sample1.shifted.bam sample1.tmp.bam
samtools index -@ 8 sample1.shifted.bam
rm sample1.tmp.bam

Working mode:

1. ATACseqSingle(): for single sample/fastq  
2. ATACseqMultiple(): for multiple sample/fastq files
3. ATACseqSample(): for single sample (could have multiple replicates)
4. ATACseq():

bam -> rmdup -> proper_paired -> peak/bw/...

##
AtacSingle()
AtacMultiple()
AtacReplicate()
Atac()

Input config file:

genome    outdir    name    fq1    fq2

Output:

- raw_data 
- clean_data
- align  
- bam_files 
- bw_files
- peak 
- motif 
- report
  - reads num
  - peaks num 
  - peaks annotation 
  - qc 
 
####

compare replicates:

- bam_files
- bw_files
- peak
- motif
- report

Mission:
atac_r1: single
atac_rn: merge replicates
atac_rx: multiple groups

"""

import os
import re
import random
import string
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trimmer
from hiseq.align.alignment import Alignment
from hiseq.peak.call_peak import Macs2
from hiseq.atac.atac_utils import *
from hiseq.fragsize.fragsize import BamPEFragSize
import hiseq
import collections


def gen_random_string(slen=10):
    return ''.join(random.sample(string.ascii_letters + string.digits, slen))


def len2(x):
    n = 0
    if x is None:
        pass
    elif isinstance(x, str):
        n += 1
    elif isinstance(x, list):
        n += len(x)
    elif isinstance(x, dict):
        n += len(x)
    else:
        print('unknown type')
    return n


def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def check_fq(fq):
    """
    Make sure
    fq: str or list, or None
    """
    if fq is None:
        # raise ValueError('fq1 required, got None')
        pass
    elif isinstance(fq, list):
        fq = file_abspath(fq)
    elif isinstance(fq, str):
        fq = [file_abspath(fq)]
    else:
        log.error('fq failed, Nont, str, list expected, got {}'.format(type(fq).__name__))

    return fq


def fq_paired(fq1, fq2):
    """
    Make sure fq1 and fq2, proper paired
    """
    fq1 = check_fq(fq1)
    fq2 = check_fq(fq2)

    if isinstance(fq1, str) and isinstance(fq2, str):
        return distance(fq1, fq2) == 1
    elif isinstance(fq1, list) and isinstance(fq2, list):
        return [distance(i, j) == 1 for i, j in zip(fq1, fq2)]
    else:
        log.warning('fq not paired: {}, {}'.format(fq1, fq2))
        return False


def init_cpu(threads=1, parallel_jobs=1):
    """
    threads, CPUs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

    max_jobs = int(n_cpu / 4.0)
    ## check parallel_jobs (max: 1/4 of n_cpus)
    if parallel_jobs > max_jobs: 
        log.warning('Too large, change parallel_jobs from {} to {}'.format(
            parallel_jobs, max_jobs))
        parallel_jobs = max_jobs

    ## check threads
    max_threads = int(0.8 * n_cpu / parallel_jobs)
    if threads * parallel_jobs > 0.8 * n_cpu:
        log.warning('Too large, change threads from {} to {}'.format(
            threads, max_threads))
        threads = max_threads

    return (threads, parallel_jobs)


def update_obj(obj, d, force=True, remove=False):
    """
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    # fresh start
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)

    return obj


class AtacReader(object):
    """
    Read config.txt/pickle
    """
    def __init__(self, x):
        self.x = x
        self.get_config() # update
        self.atacseq_type = self.args.get('atacseq_type', None)

        # auto
        self.is_atac_s1r1 = self.atacseq_type == 'atacseq_s1r1'
        self.is_atac_s1rn = self.atacseq_type == 'atacseq_s1rn'
        self.is_atac_snrn = self.atacseq_type == 'atacseq_snrn'


    def get_config(self):
        """
        locate config.pickle file
        """
        config_pickle = os.path.join(self.x, 'config', 'arguments.pickle')
        if file_exists(config_pickle):
            with open(config_pickle, 'rb') as fh:
                self.args = pickle.load(fh)
        else:
            self.args = None


class Atac(object):
    """
    Main port for ATACseq analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        # self.atacseq_type = 'atacseq_rx' # force update


    def run_atacseq_design(self):
        """
        Create design
        """
        AtacDesign(**self.__dict__).save()


    def run(self):
        """
        Run all
        """
        func = {
            'build_design': AtacRd,
            'atacseq_r1': AtacR1,
            'atacseq_rn': AtacRn,
            'atacseq_rx': AtacRx
        }
        atac = func.get(self.atacseq_type, None)
        
        if atac is None:
           raise ValueError('unknown: {}'.format(self.atacseq_type))

        # print('!AAAA-2')
        # print_dict(self.__dict__)
        # sys.exit()

        # run analysis
        atac(**self.__dict__).run()

        ## run report
        if not self.atacseq_type in ['build_design']:
            ## save arguments
            self.atacseq_type = 'atacseq_rx' # force
            chk0 = args_checker(self.__dict__, self.config_pickle)
            chk1 = args_logger(self.__dict__, self.config_txt)
    
            ## report
            args_local = self.__dict__
            args_local.update({
                'atacseq_type': 'atacseq_rx',
                'project_dir': self.outdir})
            AtacRt(**args_local).run()
        
        # print('!BBBB-1', self.atacseq_type)
        # if self.atacseq_type == 'build_design':
        #     self.run_atacseq_design()
        #     AtacRd(**self__dict__).run()
        # elif self.atacseq_type == 'atacseq_r1':
        #     AtacR1(**self.__dict__).run()
        # elif self.atacseq_type == 'atacseq_rn':
        #     AtacRn(**self.__dict__).run()
        # elif self.atacseq_type == 'atacseq_rx':
        #     AtacRx(**self.__dict__).run()
        # else:
        #     raise ValueError('unknown: {}'.format(self.atacseq_type))


class AtacR1(object):
    """
    Run ATACseq basic for single fastq file

    align - bam - rmdup - proper_paired - peak - bw - qc - report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    #######################################
    ## main pipeline
    def prep_raw(self, copy=False):
        """
        Copy raw data to dest dir
        if not: create a symlink
        
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list
        # copy
        if copy is True:
            shutil.copy(self.fq1, raw_fq1)
            shutil.copy(self.fq2, raw_fq2)
        else:
            symlink(self.fq1, raw_fq1, absolute_path=True)
            symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        args_trim = self.__dict__.copy()
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # update args
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.clean_dir,
            'library_type': 'Nextera'
        }
        args_trim.update(args_init)
        
        if trimmed is True:
            ## fq1
            if is_gz(fq1):
                symlink(fq1, clean_fq1, absolute_path=True)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    symlink(fq2, clean_fq2, absolute_path=True)
                else:
                    gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)
        else:
            if check_file(self.clean_fq_list[0]):
                log.info('trim() skipped, file exists: {}'.format(
                    self.clean_fq_list))
            else:
                trimmer = Trimmer(**args_trim)
                trimmer.run()
                fq1, fq2 = trimmer.out_files
                symlink(fq1, clean_fq1)
                symlink(fq2, clean_fq2)


    def align(self):
        """
        Alignment reads to reference genome, using bowtie2

        extra: -X 2000, ...
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'fq': fq1,
            'outdir': self.align_dir,
            'extra_para': '-X 2000 --very-sensitive-local --no-mixed', # 
            'aligner': 'bowtie2'
        }
        args_local.update(args_init)

        # output
        if check_file(args_local.get('bam_rmdup', None)):
            log.info('align() skipped, file exists: {}'.format(
                args_local.get('bam_rmdup', None)))
        else:
            Alignment(**args_local).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bam_dir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        bamlist = listfile(self.align_dir, '*.bam', recursive=True)
        bamlist = sorted(bamlist)
        bamlist = [i for i in bamlist if not i.endswith('.raw.bam')] # remove raw.bam
        # [chrM, genome]
        # 
        return(bamlist[-1]) # the last one


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        bam_raw = self.get_raw_bam()

        # rmdup
        if rmdup is True:
            Bam(bam_raw).rmdup(self.bam_rmdup)
        else:
            shutil.copy(bam_raw, self.bam_rmdup)

        # bam_proper_pair = bam
        if check_file(self.bam):
            log.info('bam_rmdup() skipped, file exists: {}'.format(
                self.bam))
        else:
            # save proper pair reads
            # Bam(self.bam_rmdup).proper_pair(self.bam)
            Bam(self.bam_rmdup).sort(self.bam)
            # index
            Bam(self.bam).index()


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def call_peak(self):
        """
        Call peaks using MACS2
        """
        args_peak = self.__dict__.copy()
        args_peak['genome_size'] = getattr(self, 'genome_size', 0)
        bam = self.bam
        bed = os.path.splitext(bam)[0] + '.bed'
        Bam(bam).to_bed(bed)
        args_peak.pop('genome', None)
        args_peak.pop('outdir', None)
        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(self.peak))
        else:
            Macs2(bed, self.genome, self.peak_dir, self.project_name, 
                atac=False).callpeak()


    def qc_lendist(self):
        """
        Create length distribution, txt
        """
        if os.path.exists(self.lendist_txt) and self.overwrite is False:
            log.info('lendist() skipped: file exists: {}'.format(
                self.lendist_txt))
        else:
            # _ = frag_length(self.bam, self.lendist_txt)
            BamPEFragSize(self.bam).saveas(self.lendist_txt)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            frip, n, total = peak_FRiP(self.peak, self.bam, threads=1)

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            # n.append('self.config.fqname')
            n.append(self.project_name)
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')


    def qc_mito(self):
        """
        Mito reads, percentage
        """
        pass


    def qc_tss(self):
        """
        TSS enrichment
        require:
        BAM, peak, TSS, ...
        """
        pass


    def qc(self):
        """
        QC for ATACseq sample
        FRiP
        Fragment size (length distribution)
        TSS
        """
        self.qc_frip()
        self.qc_lendist()
        self.qc_mito()
        self.qc_tss()


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atacseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'ATACseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(atac_report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            run_shell_cmd(cmd) 
        # try:
        #     run_shell_cmd(cmd)
        # except:
        #     log.warning('report() failed.')


    def run(self):
        """
        Run all steps for ATACseq pipeline
        """
        # init dir
        args = self.__dict__.copy()

        trimmed = args.get('trimmed', False)
        copy_raw_fq = args.get('copy_raw_fq', False)
        rmdup = args.get('rmdup', True) # remove PCR dup

        # 1. copy raw data
        self.prep_raw(copy_raw_fq)

        # 2. trim
        self.trim(trimmed)

        # 3. align
        self.align()

        # 4. rmdup
        self.get_bam_rmdup(rmdup)

        # 5. bw
        self.bam_to_bw()

        # 5. peak 
        self.call_peak()
        
        # 6. motif
        # self.motif()

        # 7. qc
        self.qc()

        # 8.report
        self.report()


class AtacRn(object):
    """
    Run ATACseq for n replicates
    input: rep_list/fq1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacRnConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def get_bam_list(self):
        """
        get the proper_mapped.bam files
        """
        return [AtacReader(i).args.get('bam', None) for i in self.rep_list]


    def get_peak_list(self):
        """
        get the .narrowPeak files
        """
        return [AtacReader(i).args.get('peak', None) for i in self.rep_list]


    def merge_bam(self):
        """
        Merge replicates, BAM
        """
        self.bam_list = self.get_bam_list()

        cmd = ' '.join([
            'samtools merge {}'.format(self.bam + '.tmp'),
            ' '.join(self.bam_list),
            '&& samtools sort -o {} {}'.format(self.bam, self.bam + '.tmp'),
            '&& samtools index {}'.format(self.bam)])

        if os.path.exists(self.bam):
            log.info('merge_bam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def call_peak(self):
        """
        Call peaks using MACS2
        ...
        """
        args_peak = self.__dict__.copy()
        args_peak['genome_size'] = getattr(self, 'genome_size', 0)

        
        bam = self.bam
        bed = os.path.splitext(bam)[0] + '.bed'
        Bam(bam).to_bed(bed)
        genome = args_peak.pop('genome', None)
        output = args_peak.pop('peak_dir', None)
        prefix = args_peak.pop('smp_name', None)

        if check_file(self.peak):
            log.info('call_peak() skipped, file exists: {}'.format(
                self.peak))
        else:
            Macs2(bed, genome, output, prefix, atac=True).callpeak()


    def get_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        args_local = {
            'bam': self.get_bam_list(),
            'outdir': self.qc_cor_dir,
            'threads': self.threads,
            'overwrite': self.overwrite,
            'binsize': self.binsize
        }
        Bam2cor(**args_local).run()


    def get_bam_fingerprint(self):
        """
        Compute the fringerprint for bam files

        plotFingerPrint (deeptools)
        """
        args_local = {
            'bam_list': self.get_bam_list(),
            'outdir': self.qc_fp_dir,
            'threads': self.threads
        }
        Bam2fingerprint(**args_local).run()


    def get_peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_overlap_dir,
            'overwrite': self.overwrite
        }   

        BedOverlap(**args_local).run()


    def get_peak_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_idr_dir,
            'overwrite': self.overwrite
        }

        PeakIDR(**args_local).run()


    def get_peak_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # get list
        self.frip_list = []
        for fq in self.fq1:
            fq_dir = os.path.join(self.outdir, fq_name(fq, pe_fix=True))
            self.frip_list.append(
                AtacReader(fq_dir).args.get('frip_txt', None))

        if not os.path.exists(self.frip_txt):
            with open(self.frip_txt, 'wt') as w:
                for f in self.frip_list:
                    with open(f) as r:
                        for line in r:
                            w.write(line)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            frip, n, total = peak_FRiP(self.peak, self.bam, threads=1)

            hd = ['FRiP', "peak_reads", "total_reads", "id"]
            n = list(map(str, [frip, n, total]))
            n.append(self.project_name)
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(hd) + '\n')
                w.write('\t'.join(n) + '\n')


    def get_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def get_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        pass


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atacseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'ATACseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        if file_exists(atac_report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            run_shell_cmd(cmd) 


    def pick_fq_samples(self, i):
        """
        Pick samples for each group
        """
        if isinstance(i, int):
            # fq1
            if isinstance(self.fq1, list):
                fq1 = self.fq1[i]
            else:
                fq1 = None

            # fq2
            if isinstance(self.fq2, list):
                fq2 = self.fq2[i]
            else:
                fq2 = None

            # rep_list
            rep_list = None

            args_fq = {
                'fq1': fq1,
                'fq2': fq2,
                'rep_list': rep_list}

            return(args_fq)
        else:
            raise ValueError('int expected, {} found'.format(type(g).__name__))


    def run_fq_single(self, i):
        """
        Run for single group
        """
        args_tmp = self.__dict__.copy()
        # required args
        args_required = ['align_to_chrM', 'aligner', 'fq1', 'fq2', 'genome', 
            'genome_size', 'trimmed', 'outdir', 'overwrite', 'parallel_jobs', 
            'threads']
        args_local = dict((k, args_tmp[k]) for k in args_required 
            if k in args_tmp)

        args_input = self.pick_fq_samples(i)
        args_init = {
            'build_design': None,
            'design': None
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) #

        obj_local = AtacR1Config(**args_local)
        AtacR1(**obj_local.__dict__).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # run each fq
        if isinstance(self.fq1, list):
            i_list = list(range(len(self.fq1)))
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_fq_single, i_list)

            # # # alternative
            # for i in i_list:
            #     self.run_fq_single(i)

        if len(self.rep_list) > 1:
            # run
            self.merge_bam()
            self.bam_to_bw()
            self.call_peak()
            # qc
            self.get_bam_cor()
            self.get_bam_fingerprint()
            self.get_peak_overlap()
            self.get_peak_idr()
            self.get_peak_frip()
            self.qc_frip()
            self.get_align_txt()
            self.get_tss_enrich()
            self.report()
        elif len(self.rep_list) == 1:
            log.warning('merge() skipped, Only 1 replicate detected')
            # copy files: bam, bw, peak
            rep_dict = AtacReader(self.rep_list[0]).args
            symlink(rep_dict.get('bam', None), self.bam)
            symlink(rep_dict.get('bw', None), self.bw)
            symlink(rep_dict.get('peak', None), self.peak)
        else:
            log.error('merge() failed, no rep detected')
            raise ValueError('merge() failed, no rep detected')


class AtacRx(object):
    """
    Run ATACseq for multiple groups
    input: groups
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacRxConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## group


    def run(self):
        """
        Run n group, n samples / group > 1
        """
        for k, v in self.fq_groups.items():
            v['build_design'] = v.get('build_design', False)
            rn_args = self.__dict__
            rn_args.update(v) # update fq1, fq2
            AtacRn(**rn_args).run()


class AtacRt(object):
    """
    Run ATACseq for multiple samples
    input: project_dir
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacRtConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atacseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'ATACseq_report.html')

        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir)

        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
    
        if file_exists(atac_report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                atac_report_html))
        else:
            run_shell_cmd(cmd) 


    def run(self):
        """
        Create report for multiple samples (n groups)
        """
        self.report()


class AtacRd(object):
    """
    Generate design.json for ATACseq
    Prepare fq for ATACseq analysis
    append/update

    format:
    input: fq1, fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_design()


    def init_args(self):
        obj_local = AtacRdConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        self.atacseq_type = 'atacseq_rd'
        self.build_design = False # force
        
        ## check
        self.design = file_abspath(self.design)
        if not isinstance(self.design, str):
            raise ValueError('--design, required')
            

    def init_design(self):
        """
        Save args in Json format
        """
        # design: init
        design_dict = {} # init
        
        if file_exists(self.design) and self.append:
            design_dict = Json(self.design).reader()
        
        # update
        # check file exists
        log.info('Check file exitence')
        for g, d in self.fq_groups.items():
            self.fq_input(d.get('fq1', None))
            self.fq_input(d.get('fq2', None))

            if d in list(design_dict.values()):
                log.warning('design exists, skipped ...')
                continue
            
            if g in design_dict:
                g = gen_random_string()
                
            # update
            design_dict[g] = d
        
        # arrange
        i = 0
        dout = {}
        for k, v in design_dict.items():
            i += 1
            key = 'atacseq_{:03d}'.format(i)
            dout[key] = v
            
        # update
        self.fq_groups = dout


    def fq_input(self, fq):
        """
        Check input file, str or list
        
        convert to absolute path
        """
        if isinstance(fq, str):
            fq = [fq] # convert to list

        if isinstance(fq, list):
            for f in fq:
                log.info('{} : {}'.format(file_exists(f), f))
                
            chk1 = file_exists(fq)
            
            if all(chk1):
                return file_abspath(fq)
            else:
                log.error('fastq files not exists')
        else:
            raise ValueError('unknown fastq, list expected, got {}'.format(
              type(fq).__name__))


    def run(self):
        # save to file
        Json(self.fq_groups).writer(self.design)


class AtacConfig(object):
    """
    Global config for ATACseq analysis

    input:

    {design|fq1,fq2|rep_list} -> {R1|Rn} -> {Rx}
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()
        self.init_mission()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'build_design': False,
            'design': None,
            'fq_dir': None,
            'fq1': None,
            'fq2': None,
            'rep_list': None,
            'smp_name': None,
            'group': None,
            'genome': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        
        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        
        # 1st level
        if not self.fq1 is None:
            self.fq1, self.fq2 = self.check_fq_paired(self.fq1, self.fq2)

        # 2nd level
        if isinstance(self.rep_list, list):
            self.rep_list = file_abspath(self.rep_list)
        elif isinstance(self.rep_list, str):
            self.rep_list = [file_abspath(self.rep_list)]
        else:
            pass

        # 3rd level
        # update group
        if file_exists(self.design):
            self.design = file_abspath(self.design)
            self.fq_groups = Json(self.design).reader()
        else:
            self.fq_groups = self.group_fq()

        self.group = list(self.fq_groups.keys())


    def group_fq(self):
        """
        separate fastq files into groups, based on filename
        """
        if isinstance(self.fq_dir, str):
            f1 = listfile(self.fq_dir, "*.fastq.gz")
            f2 = listfile(self.fq_dir, "*.fq.gz")
            f3 = listfile(self.fq_dir, "*.fastq")
            f4 = listfile(self.fq_dir, "*.fq")
            f_list = f1 + f2 + f3 + f4 # all fastq files
        elif isinstance(self.fq1, list):
            f_list = self.fq1 + self.fq2
        elif isinstance(self.fq1, str):
            f_list = [self.fq1, self.fq2]
        else:
            f_list = []
        g_list = fq_name_rmrep(f_list)
        g_list = sorted(list(set(g_list))) # unique

        # split into groups
        d = {} # 
        for g in g_list:
            # fq files for one group
            g_fq = [i for i in f_list if g in os.path.basename(i)]
            # for fq1, fq2
            r1 = re.compile('1.f(ast)?q(.gz)?')
            r2 = re.compile('2.f(ast)?q(.gz)?')
            g_fq1 = [i for i in g_fq if r1.search(i)]
            g_fq2 = [i for i in g_fq if r2.search(i)]
            # save to dict
            d[g] = {'fq1': g_fq1,
                    'fq2': g_fq2,
                    'group': g}

        return d


    def check_fq(self, fq):
        """
        Make sure
        fq: str or list, or None
        """
        if fq is None:
            # raise ValueError('fq1 required, got None')
            pass
        elif isinstance(fq, list):
            fq = file_abspath(fq)
        elif isinstance(fq, str):
            fq = [file_abspath(fq)]
        else:
            log.error('fq failed, Nont, str, list expected, got {}'.format(type(fq).__name__))

        return fq


    def check_fq_paired(self, fq1, fq2=None):
        """
        Make sure fq1 and fq2, proper paired
        """
        fq1 = self.check_fq(fq1)
        fq2 = self.check_fq(fq2)

        # fq1 = fq2
        if isinstance(fq2, list):
            if not len(fq1) == len(fq2):
                raise ValueError('fq1, fq2, files not paired correctly.')

            # check
            q_tag = 0
            for q1, q2 in zip(fq1, fq2):
                if not distance(q1, q2) == 1:
                    log.error('PE fq files failed, {}, {}'.format(q1, q2))
                    q_tag += 1

            if q_tag:
                raise ValueError('fastq files not paired correctly.')

        return (fq1, fq2)


    def init_dirs(self, create_dirs=True):
        """
        prepare directories
        """
        self.project_dir = os.path.join(self.outdir)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.config_txt = os.path.join(self.config_dir, 'arguments.txt')
        self.config_pickle = os.path.join(self.config_dir, 'arguments.pickle')
        self.config_json = os.path.join(self.config_dir, 'arguments.json')

        if create_dirs:
            check_path(self.config_dir)


    def init_mission(self):
        """
        Determine the type of ATACseq analysis
        1. build_design
        2. Single: single replicate
        3. Merge: n replicates
        4. Multiple: multiple samples/groups
        """
        self.atacseq_type = None
        # 1st level
        if self.build_design:
            self.atacseq_type = 'build_design'
        elif len(self.group) > 1: # !!!!
            self.atacseq_type = 'atacseq_rx'
        elif len(self.group) == 1:
            self.atacseq_type = 'atacseq_rn'
        elif isinstance(self.fq1, str):
            self.atacseq_type = 'atacseq_r1'
        else:
            raise ValueError('unknown atacseq type')

        # check
        if self.atacseq_type is None:
            raise ValueError('unknown atacseq_type')


class AtacR1Config(object):
    """
    Prepare directories for R1 single replicate

    require: fq1/fq2, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'trimmed': False,
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_r1' #

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # check fastq files
        if not isinstance(self.fq1, str):
            raise ValueError('AtacR1Config failed, fq1 str expected, got {}'.format(type(self.fq1).__name__))

        if not isinstance(self.fq2, str):
            raise ValueError('AtacR1Config failed, fq2 str expected, got {}'.format(type(self.fq1).__name__))

        if not fq_paired(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not paired properly: {}, {}'.format(self.fq1, self.fq2))

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, self.parallel_jobs)



    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        auto_files = {
            ## dirs
            'raw_dir': os.path.join(self.project_dir, 'raw_data'),
            'clean_dir': os.path.join(self.project_dir, 'clean_data'),
            'align_dir': os.path.join(self.project_dir, 'align'),
            'bam_dir': os.path.join(self.project_dir, 'bam_files'),
            'bw_dir': os.path.join(self.project_dir, 'bw_files'),
            'peak_dir': os.path.join(self.project_dir, 'peak'),
            'motif_dir': os.path.join(self.project_dir, 'motif'),
            'qc_dir': os.path.join(self.project_dir, 'qc'),
            'qc_lendist_dir': os.path.join(self.project_dir, 'qc', 'lendist'),
            'qc_frip_dir': os.path.join(self.project_dir, 'qc', 'FRiP'),
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'align_stat': os.path.join(self.project_dir, 'align', self.smp_name + '.align.txt'),
            # 'bam_raw': from align_dir, to-do
            'bam_rmdup': os.path.join(self.project_dir, 'bam_files', self.smp_name + '.rmdup.bam'),
            'bam_proper_pair': os.path.join(self.project_dir, 'bam_files', self.smp_name + '.proper_pair.bam'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'bed': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bed'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            # qc
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP', 'FRiP.txt'),
            }
        self = update_obj(self, auto_files, force=True) # key

        ## raw data
        self.raw_fq_list = [os.path.join(self.raw_dir, os.path.basename(self.fq1))]
        fq2_raw = None if self.fq2 is None else os.path.join(self.raw_dir, os.path.basename(self.fq2))
        self.raw_fq_list.append(fq2_raw)

        ## clean data
        self.clean_fq_list = [os.path.join(self.clean_dir, fq_name(self.fq1, pe_fix=False) + '.fq.gz')]
        fq2_clean = None if self.fq2 is None else os.path.join(self.clean_dir, fq_name(self.fq2, pe_fix=False) + '.fq.gz')
        self.clean_fq_list.append(fq2_clean)

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.raw_dir,
                self.clean_dir,
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir,
                os.path.join(self.qc_dir, 'lendist'),
                os.path.join(self.qc_dir, 'FRiP') ])


class AtacRnConfig(object):
    """
    Prepare directories/args for n replicates

    require: rep_list/fq1-list, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'smp_name': None,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'genome': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'align_to_chrM': True,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rn'

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())            
        self.outdir = file_abspath(self.outdir)

        # check fastq files/build rep_list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]

        if isinstance(self.fq1, list):
            self.rep_list = [os.path.join(self.outdir, fq_name(i, pe_fix=True)) for i in self.fq1]

        if isinstance(self.rep_list, list):
            self.smp_name = fq_name_rmrep(self.rep_list)
            self.smp_name = self.smp_name.pop() # to str !!! to-do, check unique
        else:
            raise ValueError('AtacRnConfig failed, fq1/rep_list expect list, got {}'.format(type(self.rep_list).__name__))

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, self.parallel_jobs)



    def init_dirs(self, create_dirs=True):
        """
        for n replicate fastq file
        prepare directories
        """
        self.rep_list = file_abspath(self.rep_list)

        # path, files
        self.project_name = self.smp_name #
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        auto_files = {
            ## dirs
            'align_dir': os.path.join(self.project_dir, 'align'),
            'bam_dir': os.path.join(self.project_dir, 'bam_files'),
            'bw_dir': os.path.join(self.project_dir, 'bw_files'),
            'peak_dir': os.path.join(self.project_dir, 'peak'),
            'motif_dir': os.path.join(self.project_dir, 'motif'),
            'qc_dir': os.path.join(self.project_dir, 'qc'),
            'qc_lendist_dir': os.path.join(self.project_dir, 'qc', 'lendist'),
            'qc_frip_dir': os.path.join(self.project_dir, 'qc', 'FRiP'),
            'qc_idr_dir': os.path.join(self.project_dir, 'qc', 'IDR'),
            'qc_cor_dir': os.path.join(self.project_dir, 'qc', 'cor'),
            'qc_fp_dir': os.path.join(self.project_dir, 'qc', 'fingerprint'),
            'qc_overlap_dir': os.path.join(self.project_dir, 'qc', 'overlap'),
            'report_dir': os.path.join(self.project_dir, 'report'),
            ## files
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'bam': os.path.join(self.project_dir, 'bam_files', self.project_name + '.bam'),
            'peak': os.path.join(self.project_dir, 'peak', self.project_name + '_peaks.narrowPeak'),
            'bw': os.path.join(self.project_dir, 'bw_files', self.project_name + '.bigWig'),
            ## qc files
            'lendist_txt': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.txt'),
            'lendist_pdf': os.path.join(self.project_dir, 'qc', 'lendist', 'length_distribution.pdf'),
            'frip_txt': os.path.join(self.project_dir, 'qc', 'FRiP', 'FRiP.txt'),
            'cor_npz': os.path.join(self.project_dir, 'qc', 'cor', 'cor.bam.npz'),
            'cor_counts': os.path.join(self.project_dir, 'qc', 'cor', 'cor.bam.counts.tab'),
            'idr_txt': os.path.join(self.project_dir, 'qc', 'IDR', 'dir.txt'),
            'idr_log': os.path.join(self.project_dir, 'qc', 'IDR', 'dir.log'),
            'peak_overlap_pdf': os.path.join(self.project_dir, 'qc', 'overlap', 'peak_overlap_pdf')        
            }
        self = update_obj(self, auto_files, force=True) # key

        if create_dirs:
            check_path([
                self.project_dir, 
                self.config_dir, 
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.peak_dir, 
                self.qc_dir, 
                self.report_dir,
                os.path.join(self.qc_dir, 'lendist'),
                os.path.join(self.qc_dir, 'FRiP'),
                os.path.join(self.qc_dir, 'cor'),
                os.path.join(self.qc_dir, 'IDR'),
                os.path.join(self.qc_dir, 'overlap')])


class AtacRxConfig(object):
    """
    Run ATACseq for n groups
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'outdir': None}
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rx'

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # groups
        if not isinstance(self.group, list):
            raise ValueError('group, list expected, got {}'.format(
                type(self.group).__name__))


class AtacRtConfig(object):
    """
    Give report for ATACseq directories

    output: config, report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'project_dir': None}
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rt'


        # get sample list from project_dir !!!!
        if not isinstance(self.project_dir, str):
            raise ValueError('project_dir, str expected, got {}'.format(
                type(self.project_dir).__name__))


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        # path, files
        self.report_dir = os.path.join(self.project_dir, 'report')

        if create_dirs:
            check_path(self.report_dir)


class AtacRdConfig(object):
    """
    Generate ATACseq design.json

    output: fq_dir/fq1,fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_dirs()


    def init_args(self):
        """
        required arguments for ATACseq design
        """
        args_init = {
          'fq_dir': None,
          'fq1': None,
          'fq2': None,
          'design': None,
          'append': True}
        self = update_obj(self, args_init, force=False)
        
        # fastq files
        self.fq_groups = self.group_fq()


    def init_dirs(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        pass
            

    def group_fq(self):
        """
        separate fastq files into groups, based on filename
        """
        if isinstance(self.fq_dir, str):
            f1 = listfile(self.fq_dir, "*.fastq.gz")
            f2 = listfile(self.fq_dir, "*.fq.gz")
            f3 = listfile(self.fq_dir, "*.fastq")
            f4 = listfile(self.fq_dir, "*.fq")
            f_list = f1 + f2 + f3 + f4 # all fastq files
        elif isinstance(self.fq1, list):
            f_list = self.fq1 + self.fq2
        elif isinstance(self.fq1, str):
            f_list = [self.fq1, self.fq2]
        else:
            f_list = []
            
        # check fastq files
        if len(f_list) < 1:
            raise ValueError('fq_dir, fq1,fq2, fq files not found')
            
        g_list = fq_name_rmrep(f_list)
        g_list = sorted(list(set(g_list))) # unique

        # split into groups
        d = {} # 
        for g in g_list:
            g_fq = [i for i in f_list if g in os.path.basename(i)]
            # for fq1, fq2
            r1 = re.compile('1.f(ast)?q(.gz)?')
            r2 = re.compile('2.f(ast)?q(.gz)?')
            g_fq1 = [i for i in g_fq if r1.search(i)]
            g_fq2 = [i for i in g_fq if r2.search(i)]
            ## pairing
            if len(g_fq2) > 0: # PE reads
                if not all(fq_paired(g_fq1, g_fq2)):
                    raise ValueError('fq not paired')
            # save to dict
            d[g] = {'fq1': g_fq1,
                    'fq2': g_fq2,
                    'group': g}

        return d

