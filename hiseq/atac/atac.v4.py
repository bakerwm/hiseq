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
import glob
import random
import string
import collections
import hiseq
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trim
from hiseq.align.alignment import AlignIndex
# from hiseq.peak.call_peak import Macs2
from hiseq.atac.atac_utils import *
from hiseq.fragsize.fragsize import BamPEFragSize


def gen_random_string(slen=10):
    return ''.join(random.sample(string.ascii_letters + string.digits, slen))


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


class Atac(object):
    """
    Main port for ATACseq analysis

    if fq_groups >= 1:
        (AtacRx <- AtacRn <- AtacR1) <- AtacRp

    if fq1, fq2: (str)
        AtacR1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)

        ## save arguments
        # chk0 = args_checker(self.__dict__, self.config_pickle)
        # chk1 = args_logger(self.__dict__, self.config_txt)


    def run(self):
        """
        Run all
        """
        if self.atacseq_type == 'build_design':
            AtacRd(**self.__dict__).run()
        elif self.atacseq_type == 'atacseq_r1':
            AtacR1(**self.__dict__).run()
        elif self.atacseq_type == 'atacseq_rn':
            AtacRn(**self.__dict__).run()
        elif self.atacseq_type == 'atacseq_rx':
            AtacRx(**self.__dict__).run()
        else:
            raise ValueError('unknown atacseq_type: {}'.format(self.atacseq_type))


class AtacRx(object):
    """
    Run ATACseq for multiple groups
    input: self.fq_groups
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


    def run_group_single(self, d):
        """
        d dict, contains fq1, fq2
        the other args were saved in self (global)
        """
        args_local = self.__dict__.copy()

        args_required = {
            'parallel_jobs': 1, # force
            'build_design': False,            
            'design': None,
            'fq_dir': None,
            'rep_list': d.get('rep_list', None),
            'fq1': d.get('fq1', None),
            'fq2': d.get('fq2', None)
        }
        args_local.update(args_required) # force
        AtacRn(**args_local).run()


    def run(self):
        """
        self.fq_groups
          - key: string
          - value: fq1/fq2/rep_list
        """
        if self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_group_single, self.fq_groups.values())
        else:
            for d in list(self.fq_groups.values()):
                self.run_group_single(d)

        # run report for all samples
        args_local = self.__dict__.copy()
        AtacRp(**args_local).run()


class AtacRn(object):
    """
    Run ATACseq for multiple replicates, merge replicates
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
        return [AtacReader(i).args.get('bam_rmdup', None) for i in self.rep_list]


    def get_bw_list(self):
        """
        get the proper_mapped.bam files
        """
        return [AtacReader(i).args.get('bw', None) for i in self.rep_list]


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
            'samtools merge -',
            ' '.join(self.bam_list),
            '| samtools sort -o {} -'.format(self.bam),
            '&& samtools index {}'.format(self.bam)])

        if os.path.exists(self.bam):
            log.info('merge_bam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        # remove dup
        if rmdup:
            if file_exists(self.bam):
                if file_exists(self.bam_rmdup) and not self.overwrite:
                    log.info('rmdup() skipped, file exists: {}'.format(
                        self.bam_rmdup))
                else:
                    Bam(self.bam).rmdup(self.bam_rmdup)
                    Bam(self.bam_rmdup).index()
        else:
            file_symlink(self.bam, self.bam_rmdup)

        ## calculate norm scale
        self.align_scale = self.cal_norm_scale(self.bam_rmdup)
        with open(self.align_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.align_scale))


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam_rmdup,
            'scaleFactor': self.align_scale,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 0,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def bam_to_bg(self):
        """
        Convert bam to bedgraph
        
        norm
        bedtools genomecov -bg -scale $scale_factor -ibam bam > bg
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedtools')),
            'genomecov -bg -scale {}'.format(self.align_scale),
            '-ibam {}'.format(self.bam_rmdup),
            '| sort -k1,1 -k2,2n > {}'.format(self.bg)
            ])

        if file_exists(self.bg) and not self.overwrite:
            log.info('bam_to_bg() skipped, file exists:{}'.format(self.bg))
        else:
            cmd_txt = self.bg_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bam_to_bg() failed')


    def bg_to_bw(self):
        """
        Create bigWig
        bedgraph to bigWig
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedGraphToBigWig')),
            '{} {} {}'.format(self.bg, self.genome_size_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_txt = self.bw_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def call_peak(self):
        """
        Call peaks using MACS2/SEACR
        """
        args_local = {
            'ip': self.bam_rmdup,
            'input': None,
            'outdir': self.peak_dir,
            'prefix': self.smp_name,
            'genome': self.genome,
            'genome_size': self.genome_size,
            'genome_size_file': self.genome_size_file
        }
        CallPeak(method='macs2', **args_local).run()
        CallPeak(method='seacr', **args_local).run()


    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor
        
        Bam().count_reads()
        """
        bam = Bam(bam)
        is_pe = bam.isPaired()
        n_map = bam.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    ## quality control ##
    def qc_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        args_local = {
            'bam': self.get_bam_list(),
            'outdir': self.qc_dir,
            'prefix': '06.bam_cor',
            'threads': self.threads,
            'overwrite': self.overwrite,
            'binsize': self.binsize
        }
        Bam2cor(**args_local).run()


    def qc_peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_dir,
            'prefix': '08.peak_overlap',
            'overwrite': self.overwrite
        }   

        BedOverlap(**args_local).run()


    def qc_peak_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        args_local = {
            'peak': self.get_peak_list(),
            'outdir': self.qc_dir,
            'prefix': '07.peak_idr',
            'overwrite': self.overwrite
        }

        PeakIDR(**args_local).run()


    def get_peak_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # get list
        self.frip_list = [AtacReader(i).args.get('frip_txt', None) for \
            i in self.rep_list]
        self.frip_list = [i for i in self.frip_list if file_exists(i)]

        if len(self.frip_list) > 0:
            with open(self.frip_txt, 'wt') as w:
                w.write('\t'.join(['total', 'peak_reads', 'FRiP', 'id']) + '\n')
                for f in self.frip_list:
                    with open(f) as r:
                        for line in r:
                            if 'FRiP' in line:
                                continue
                            w.write(line)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            # total, n, frip
            s = PeakFRiP(peak=self.peak, bam=self.bam_rmdup, method='featureCounts').run()
            s.append(self.project_name)
            s = list(map(str, s))

            hd = ['total', 'peak_reads', 'FRiP', 'id']
            with open(self.frip_txt, 'wt') as w:
                w.write('\n'.join([
                    '\t'.join(hd),
                    '\t'.join(s)
                    ])+'\n')


    def qc_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def qc_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        bw_list = self.get_bw_list()
        bw_list.append(self.bw)
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
            '-o {}'.format(self.tss_enrich_matrix),
            '--referencePoint TSS',
            '-b 2000 -a 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss_enrich() skipped, file exists: {}'.format(
                self.tss_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.tss_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_tss_enrich() failed, see: {}'.format(
                        self.tss_enrich_matrix_log))


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        bw_list = self.get_bw_list()
        bw_list.append(self.bw)
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
            '-o {}'.format(self.genebody_enrich_matrix),
            '-b 2000 -a 2000 --regionBodyLength 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.genebody_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.genebody_enrich_matrix),
            '-o {}'.format(self.genebody_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.genebody_enrich_png) and not self.overwrite:
            log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
                self.genebody_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.genebody_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')
                    
                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_genebody_enrich() failed, see: {}'.format(
                        self.genebody_enrich_matrix_log))


    def qc_bam_fingerprint(self):
        """
        Calculate fingerprint for bam files
        """
        Bam2fingerprint(
            bam_list=self.get_bam_list(),
            prefix='09.fingerprint',
            outdir=self.qc_dir,
            threads=self.threads).run()


    def report(self):
        """
        Create report for one sample
        html
        """
        # run report for all samples
        args_local = self.__dict__.copy()
        AtacRp(**args_local).run()

        # pkg_dir = os.path.dirname(hiseq.__file__)
        # qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        # hiseq_report_html = os.path.join(
        #     self.report_dir, 
        #     'HiSeq_report.html')

        # cmd = ' '.join([
        #     'Rscript',
        #     qc_reportR,
        #     self.project_dir,
        #     self.report_dir,
        #     '1>',
        #     self.report_log
        #     ])

        # cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        # with open(cmd_txt, 'wt') as w:
        #     w.write(cmd + '\n')

        # if file_exists(hiseq_report_html):
        #     log.info('report() skipped, file exists: {}'.format(
        #         hiseq_report_html))
        # else:
        #     run_shell_cmd(cmd) 


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
        args_required = ['aligner', 'fq1', 'fq2', 'genome', 'gene_bed',
            'genome_size', 'is_ip', 'trimmed', 'outdir', 'overwrite', 
            'parallel_jobs', 'threads', 'spikein', 'spikein_index',
            'extra_index', 'cut_to_length', 'recursive']
        args_local = dict((k, args_tmp[k]) for k in args_required 
            if k in args_tmp)

        args_input = self.pick_fq_samples(i)
        args_init = {
            'build_design': None,
            'design': None,
            'gene_bed': self.gene_bed
        }
        args_local.update(args_input) # update fq1/rep_list/group
        args_local.update(args_init) #

        AtacR1(**args_local).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        # run each fq
        if isinstance(self.fq1, list):
            i_list = list(range(len(self.fq1)))
            # in parallel
            if self.parallel_jobs > 1:
                with Pool(processes=self.parallel_jobs) as pool:
                    pool.map(self.run_fq_single, i_list)
            else:
                for i in i_list:
                    self.run_fq_single(i)

            # # # alternative
            # for i in i_list:
            #     self.run_fq_single(i)

        if len(self.rep_list) > 1:
            # run
            self.merge_bam()
            self.get_bam_rmdup()
            # self.bam_to_bw()
            self.bam_to_bg()
            self.bg_to_bw()
            self.call_peak()
            # qc
            self.qc_bam_cor()
            self.qc_peak_overlap()
            self.qc_peak_idr()
            # self.get_peak_frip()
            self.qc_frip()
            self.qc_align_txt()
            self.qc_tss_enrich()
            self.qc_genebody_enrich()
            self.qc_bam_fingerprint()
            self.report()
        elif len(self.rep_list) == 1:
            # new files, symlinks
            log.warning('merge() skipped, Only 1 replicate detected')
            # copy files: bam, bw, peak
            rep = AtacReader(self.rep_list[0]).args
            bam = rep.get('bam_rmdup', None)
            bg = rep.get('bg', None)
            bw = rep.get('bw', None)
            peak = rep.get('peak', None)
            peak_seacr = rep.get('seacr', None)

            align_scale_txt = rep.get('align_scale_txt', None)
            align_json = rep.get('align_json')

            file_symlink(align_scale_txt, self.align_scale_txt)
            file_symlink(align_json, self.align_json)
            file_symlink(bam, self.bam_rmdup)
            Bam(self.bam_rmdup).index()
            file_symlink(bg, self.bg)
            file_symlink(bw, self.bw)
            file_symlink(peak, self.peak)
            # file_symlink(peak_seacr, self.peak_seacr)

            # self.get_peak_frip()
            self.qc_frip()
            self.qc_tss_enrich()
            self.qc_genebody_enrich()
            self.qc_bam_fingerprint()
            self.report()

        else:
            log.error('merge() failed, no rep detected')
            raise ValueError('merge() failed, no rep detected')


class AtacR1(object):
    """
    Run ATACseq basic for single fastq file (ip/input)

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


    # main pipeline #
    def prep_raw(self):
        """
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list
        file_symlink(self.fq1, raw_fq1, absolute_path=True)
        file_symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Cut the reads to specific length, from the 3' end

        For ATACseq, CUT&RUN, CUT&TAG experiment

        suggest=50

        if not: create a symlink

        Trim reads:
        Trim(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if not:
            copy/links
        """
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # whether to trim or not
        if trimmed:
            ## copy files
            file_symlink(fq1, clean_fq1)
            file_symlink(fq2, clean_fq2)        
        else:
            # update args
            args_local = self.__dict__.copy()
            args_init = {
                'fq1': fq1,
                'fq2': fq2,
                'outdir': self.clean_dir,
                'library_type': None,
                'len_min': 20,
                'cut_to_length': self.cut_to_length,
                'recursive': self.recursive,
                'parallel_jobs': 1 # do not allowed > 1
            }
            args_local.update(args_init)
            trim = Trim(**args_local)
            trim.run()

            ## copy files
            file_symlink(trim.clean_fq1, clean_fq1)
            file_symlink(trim.clean_fq2, clean_fq2)


    def align_genome(self):
        """
        Alignment PE reads to reference genome, using bowtie2
        
        --sensitive --local -X 2000 

        samtools view -bhS -f 2 -F 1804
        
        include: -f
            2    0x2    Read mapped in proper pair

        exclude: -F
            4    0x4    Read unmapped, 
            8    0x8    Mate unmapped,
            256  0x100  Not primary alignment, 
            512  0x200  Reads fails platform quality checks,
            1024 0x400  Read is PCR or optical duplicate
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.align_dir,
            'genome': self.genome,
            'extra_index': self.extra_index,
            'aligner': self.aligner,
            'keep_tmp': self.keep_tmp,
            'max_fragment': 2000 
        }
        args_local.update(args_init)

        # output
        bam = args_local.get('bam')
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            Align(**args_local).run()

        # copy files
        g = AtacReader(self.align_dir).args
        g_bam = g.get('bam', None)
        g_align_bam = g.get('bam', None)
        g_align_stat = g.get('align_stat', None)
        g_align_json = g.get('align_json', None)
        g_align_flagstat = g.get('align_flagstat', None)

        # copy files
        file_symlink(g_align_bam, self.bam)
        file_copy(g_align_stat, self.align_stat)
        file_copy(g_align_json, self.align_json)
        file_copy(g_align_flagstat, self.align_flagstat)


    def align_spikein(self):
        """
        Alignment PE reads, spikein
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.spikein_dir,
            'genome': None,
            'index': self.spikein_index,
            'threads': self.threads,
            'keep_tmp': self.keep_tmp
        }

        # output
        bam = args_local.get('spikein_bam')
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            Align(**args_sp).run()

        # copy files
        sp = AtacReader(self.spikein_dir).args
        sp_bam = sp.get('bam', None)
        sp_align_bam = sp.get('bam', None)
        sp_align_stat = sp.get('align_stat', None)
        sp_align_json = sp.get('align_json', None)
        sp_align_flagstat = sp.get('align_flagstat', None)

        # copy files
        file_copy(sp_align_bam, self.spikein_bam)
        file_copy(sp_align_stat, self.spikein_stat)
        file_copy(sp_align_json, self.spikein_json)
        file_copy(sp_align_flagstat, self.spikein_flagstat)

        ## calculate norm scale
        self.spikein_scale = self.cal_norm_scale(self.spikein_bam, 10000)
        with open(self.spikein_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.spikein_scale))


    # Deprecated
    def get_bam(self, dir):
        """
        Get the align bam file
        """
        bam_list = listfile(dir, '*.bam', recursive=True)
        bam_list = sorted(bam_list)
        bam_list = [i for i in bam_list if not i.endswith('.raw.bam')] #

        return bam_list[-1]


    def get_bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        # remove dup
        if rmdup:
            Bam(self.bam).rmdup(self.bam_rmdup)
            Bam(self.bam_rmdup).index()
        else:
            file_symlink(self.bam, self.bam_rmdup)

        ## calculate norm scale
        self.align_scale = self.cal_norm_scale(self.bam_rmdup)
        with open(self.align_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.align_scale))


    def call_peak(self):
        """
        Call peaks using MACS2/SEACR
        
        -f BAMPE

        """
        args_local = {
            'ip': self.bam_rmdup, # rm_dup
            'input': None,
            'outdir': self.peak_dir,
            'prefix': self.smp_name,
            'genome': self.genome,
            'genome_size': self.genome_size,
            'genome_size_file': self.genome_size_file
        }
        CallPeak(method='macs2', **args_local).run()
        CallPeak(method='seacr', **args_local).run()


    def bam_to_bg(self):
        """
        Convert bam to bedgraph
        
        norm
        bedtools genomecov -bg -scale $scale_factor -ibam bam > bg
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedtools')),
            'genomecov -bg -scale {}'.format(self.align_scale),
            '-ibam {}'.format(self.bam_rmdup),
            '| sort -k1,1 -k2,2n > {}'.format(self.bg)
            ])

        if file_exists(self.bg) and not self.overwrite:
            log.info('bam_to_bg() skipped, file exists:{}'.format(self.bg))
        else:
            cmd_txt = self.bg_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bam_to_bg() failed')


    def bg_to_bw(self):
        """
        Create bigWig
        bedgraph to bigWig
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedGraphToBigWig')),
            '{} {} {}'.format(self.bg, self.genome_size_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_txt = self.bw_dir + '/cmd.txt'
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')

    # quality control #
    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor
        
        Bam().count_reads()
        """
        bam_o = Bam(bam)
        is_pe = bam_o.isPaired()
        n_map = bam_o.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    def qc_trim_summary(self):
        """
        reads trim off
        
        # version-1
        #sample input   output  percent
        fq_rep1      2234501 2234276 99.99%

        # version-2
        #name   total   too_short       dup     too_short2      clean   percent
        """
        if file_exists(self.trim_stat_txt):
            with open(self.trim_stat_txt, 'rt') as r:
                lines = r.readlines()
            lines = [i for i in lines if not i.startswith('#')]
            line = lines.pop() # last line
            tabs = line.strip().split('\t')
            fqname, n_input = tabs[:2]
            n_output, n_pct = tabs[-2:]
            d = {
                'id': fqname,
                'input': n_input,
                'output': n_output,
                'out_pct': n_pct,
                'rm_pct': 100 - float(n_pct.strip())
            }
        else:
            d = {
                'id': self.smp_name,
                'input': 0,
                'output': 0,
                'out_pct': 100,
                'rm_pct': 0
            }
        Json(d).writer(self.trim_summary_json)


    def qc_align_summary(self):
        """
        Save summary to file
        """
        # align to genome
        align = Json(self.align_json).reader()

        # duplicates
        nodup = Bam(self.bam_rmdup)
        is_pe = nodup.isPaired()
        n_nodup = nodup.getNumberOfAlignments()
        if is_pe:
            n_nodup = n_nodup/2
        
        # spikein 
        if file_exists(self.spikein_bam):
            n_spikein = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_spikein = n_spikein/2
        else:
            n_spikein = 0

        # Json(align).writer()
        align['nodup'] = n_nodup
        align['spikein'] = n_spikein
        align['chrM'] = self.qc_mito()

        # save to file
        Json(align).writer(self.align_summary_json)


    def qc_lendist(self):
        """
        Create length distribution, txt
        """
        if os.path.exists(self.lendist_txt) and self.overwrite is False:
            log.info('lendist() skipped: file exists: {}'.format(
                self.lendist_txt))
        else:
            BamPEFragSize(self.bam_rmdup).saveas(self.lendist_txt)


    def qc_frip(self):
        """
        Compute FRiP
        """
        if check_file(self.frip_txt):
            log.info('qc_frip() skipped, file exists: {}'.format(
                self.frip_txt))
        else:
            try:
                # total, n, frip
                s = PeakFRiP(peak=self.peak, bam=self.bam_rmdup, 
                    method='featureCounts').run()
            except:
                log.error('PeakFRiP() failed, see ')
                s = [0, 0, 0]

            s.append(self.project_name)
            s = list(map(str, s))
            hd = ['total', 'peak_reads', 'FRiP', 'id']

            with open(self.frip_txt, 'wt') as w:
                w.write('\n'.join([
                    '\t'.join(hd),
                    '\t'.join(s)
                    ])+'\n')


    def qc_mito(self):
        """
        Mito reads, percentage
        """
        lines = pysam.idxstats(self.bam_rmdup).splitlines()
        # extract chrM, MT
        lines = [i for i in lines if re.search('^(chrM|MT)', i)]

        try:
            nreads = sum(
                map(int, [x.split("\t")[2]
                          for x in lines if not x.startswith("#")]))

        except IndexError as msg:
            raise IndexError(
                "can't get number of reads from bamfile, msg=%s, data=%s" %
                (msg, lines))

        return nreads


    def qc_tss(self):
        """
        TSS enrichment

        computeMatrix referencepoint -b -R gene.bed -S in.bw -o mat.gz

        require:
        plotHeatmap()
        plotProfile()
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '--referencePoint TSS',
            '-b 2000 -a 2000 --binSize 10 --sortRegions descend --skipZeros',
            '--samplesLabel {}'.format(self.smp_name),
            '-p {}'.format(self.threads),
            '--regionsFileName {}'.format(self.gene_bed),
            '--scoreFileName {}'.format(self.bw),
            '-o {}'.format(self.tss_enrich_matrix),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss() skipped, file exists')
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                cmd_txt = self.qc_dir + '/tss_enrich.sh'
                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('failed to run computeMatrix, see: {}'.format(
                        self.tss_enrich_matrix_log ))


    def qc_tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'reference-point',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(self.bw),
            '-o {}'.format(self.tss_enrich_matrix),
            '--referencePoint TSS',
            '-b 2000 -a 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.tss_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.tss_enrich_matrix),
            '-o {}'.format(self.tss_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.tss_enrich_png) and not self.overwrite:
            log.info('qc_tss_enrich() skipped, file exists: {}'.format(
                self.tss_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.tss_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_tss_enrich() failed, see: ')


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(self.bw),
            '-o {}'.format(self.genebody_enrich_matrix),
            '-b 2000 -a 2000 --regionBodyLength 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.genebody_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.genebody_enrich_matrix),
            '-o {}'.format(self.genebody_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.genebody_enrich_png) and not self.overwrite:
            log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
                self.genebody_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.genebody_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')
                    
                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_genebody_enrich() failed, see: ')


    def pipe_genome(self):
        """
        Run for genome
        """
        self.prep_raw()
        self.trim(trimmed=self.trimmed)
        self.align_genome()
        self.get_bam_rmdup()
        self.call_peak()
        self.bam_to_bg()
        self.bg_to_bw()


    def pipe_spikein(self):
        """
        Run for spikein
        """
        if isinstance(self.spikein_index, str):
            self.align_spikein() # run alignment


    def pipe_qc(self):
        """
        Quality control
        """
        self.qc_trim_summary()
        self.qc_align_summary()
        self.qc_lendist()
        self.qc_frip()
        self.qc_mito()
        self.qc_tss_enrich()
        self.qc_genebody_enrich()


    def report(self):
        """
        Create report for one sample
        html
        """
        # run report for all samples
        args_local = self.__dict__.copy()
        AtacRp(**args_local).run()

        # pkg_dir = os.path.dirname(hiseq.__file__)
        # qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        # hiseq_report_html = os.path.join(
        #     self.report_dir, 
        #     'HiSeq_report.html')

        # cmd = ' '.join([
        #     'Rscript',
        #     qc_reportR,
        #     self.project_dir,
        #     self.report_dir,
        #     '1>',
        #     self.report_log
        #     ])

        # cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        # with open(cmd_txt, 'wt') as w:
        #     w.write(cmd + '\n')
    
        # if file_exists(hiseq_report_html):
        #     log.info('report() skipped, file exists: {}'.format(
        #         hiseq_report_html))
        # else:
        #     run_shell_cmd(cmd)


    def run(self):
        """
        Run all steps for ChIPseq pipeline
        """
        # init dir
        args = self.__dict__.copy()

        self.pipe_genome()
        self.pipe_spikein()
        self.pipe_qc()
        self.report()

        # # remove clean fastq files
        # del_list = self.clean_fq_list
        # if not self.keep_tmp:
        #     file_remove(del_list, ask=False)


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
        Save args in Json format， fq_groups
        """
        # design: init
        if file_exists(self.design):
            design_dict = Json(self.design).reader()
        else:
            design_dict = {}

        # update
        # check file exists
        log.info('Check file exitence')
        for g, d in self.fq_groups.items():
            self.fq_input(d.get('fq1', None))
            self.fq_input(d.get('fq2', None))

            if d in design_dict.values():
                log.warning('design exists, skipped ...')
                continue
            
            if g in design_dict:
                g = gen_random_string()
                
            # update
            design_dict[g] = d
        
        # assign group number
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


class AtacRp(object):
    """
    Run ATACseq for multiple samples
    input: project_dir
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = AtacRpConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        atac_report_html = os.path.join(
            self.report_dir, 
            'HiSeq_report.html')

        cmd = ' '.join([
            'Rscript',
            qc_reportR,
            self.project_dir,
            self.report_dir,
            '1>',
            self.report_log
            ])

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


class AtacConfig(object):
    """
    Global config for CnR analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_mission()


    def init_args(self):
        """
        required arguments for ChIPseq analysis
        """
        args_init = {
            'build_design': False,
            'design': None,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': 0,
            'gene_bed': None,
            'keep_tmp': None,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)

        # outdir
        if self.outdir is None:
          self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        # fq groups
        self.fq_groups = Json(self.design).reader() if \
            file_exists(self.design) else self.group_fq()


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_dir = self.outdir
        self.config_dir = self.project_dir + '/config'

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json'
        }
        self = update_obj(self, default_files, force=True) # key

        if create_dirs:
            check_path(self.config_dir)


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


    def init_mission(self):
        """
        Determine the type of ATACseq analysis
        1. build_design
        2. Single: IP/Input pair
        3. Multiple: IP/Input pair
        4. x, ip vs input
        """
        self.atacseq_type = None

        # 1st level
        if self.build_design:
            self.atacseq_type = 'build_design'
        elif len(self.fq_groups) >= 1:
            self.atacseq_type = 'atacseq_rx'
        # elif len(self.fq_groups) == 1:
        #     self.atacseq_type = 'atacseq_rn'
        elif isinstance(self.fq1, str):
            self.atacseq_type = 'atacseq_r1'
        else:
            raise ValueError('unknown ATACseq type')

        # check
        if self.atacseq_type is None:
            raise ValueError('unknown ATACseq')


class AtacRxConfig(object):
    """
    Prepare directories for Rx: N groups

    require: groups > 1
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        """
        Required:
          - design, str
          - fq_dir, str
          - rep_list, list
          - fq1, list
          - fq2, list
        """
        args_init = {
            'design': None,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': 0,
            'gene_bed': None,
            'keep_tmp': None,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rx' #

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # groups
        self.fq_groups = Json(self.design).reader() if \
            file_exists(self.design) else self.group_fq()

        if len(self.fq_groups) < 1:
            raise ValueError('fastq groups failed')


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
            raise ValueError('fq_dir, fq1, fq2, fastq files not found')
            
        g_list = fq_name_rmrep(f_list)
        g_list = sorted(list(set(g_list))) # unique group

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


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_dir = os.path.join(self.outdir)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.report_dir = os.path.join(self.project_dir, 'report')

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json'
        }
        self = update_obj(self, default_files, force=True) # key

        self.report_log = self.report_dir + '/report.log',  

        check_path([self.config_dir, self.report_dir])


class AtacRnConfig(object):
    """
    Prepare directories for R1 single replicates

    Required:
      - fq1: list
      - fq2: list
      - rep_list: list
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        """
        required arguments for CnR analysis
        """
        args_init = {
            'smp_name': None,
            'group': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': 0,
            'keep_tmp': None,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rn' # 

        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # smp_name
        self.group = fq_name_rmrep(self.rep_list).pop()
        self.smp_name = self.group

        # aligner
        if self.aligner is None:
            self.aligner = 'bowtie2'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads, 
            self.parallel_jobs)

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)


    def init_fq(self):
        """
        Required:
          - fq1, list
          - fq2, list
          - rep_list, list
        """
        # rep_list
        if isinstance(self.rep_list, list):
            pass
        else:
            # fq1
            if not isinstance(self.fq1, list):
                raise ValueError('--fq1, list expected, got {}'.format(
                    type(self.fq1).__name__))

            if not all(file_exists(self.fq1)):
                raise ValueError('--fq1, file not exists: {}'.format(
                    self.fq1))

            # fq2
            if not isinstance(self.fq2, list):
                raise ValueError('--fq2, list expected, got {}'.format(
                    type(self.fq2).__name__))

            if not all(file_exists(self.fq2)):
                raise ValueError('--fq2, file not exists: {}'.format(
                    self.fq2))

            if not all(fq_paired(self.fq1, self.fq2)):
                raise ValueError('fq1, fq2 not paired properly: {}, {}'.format(
                    self.fq1, self.fq2))

            # update rep_list
            fq_names = [fq_name(i, pe_fix=True) for i in self.fq1]
            self.rep_list = [self.outdir + '/' + i for i in fq_names]

            self.fq1 = file_abspath(self.fq1)
            self.fq2 = file_abspath(self.fq2)


    def init_index(self):
        """
        alignment index
        genome, extra_index, genome_index

        output: genome_index, spikein_index
        """
        # check genome size
        if isinstance(self.extra_index, str):
            self.genome_index = self.extra_index

        ## check genome index
        if isinstance(self.genome_index, str):
            ai = AlignIndex(index=self.genome_index)

            if self.genome_size < 1:
                self.genome_size = ai.index_size()

            if not isinstance(self.genome_size_file, str): 
                self.genome_size_file = ai.index_size(return_file=True)

        elif isinstance(self.genome, str):
            self.genome_index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group='genome')

            self.genome_size_file = Genome(genome=self.genome).get_fasize()

            with open(self.genome_size_file, 'rt') as r:
                s = [i.strip().split('\t')[1] for i in r.readlines()]
            
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))
            
        else:
            raise ValueError('index failed, extra_index, genome, genome_index')

        # check spikein index
        if isinstance(self.spikein_index, str):
            ai = AlignIndex(index=self.spikein_index, aligner=self.aligner)

            if not ai.is_index():
                raise ValueError('spikein_index failed, {}'.format(
                    self.spikein_index))

        elif isinstance(self.spikein, str):
            self.spikein_index = AlignIndex(aligner=self.aligner).search(
                genome=self.spikein, group='genome')

        else:
            self.spikein_index = None


    def init_files(self, create_dirs=True):
        """
        Prepare directories, files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        default_dirs = {
            'align_dir': 'align',
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'bg_dir': 'bg_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
            'report_log': self.report_dir + '/report.log',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',            
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'align_scale_txt': self.bam_dir + '/' + 'scale.txt',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.bowtie2.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.bowtie2.json',

            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_txt': self.qc_dir + '/03.FRiP.txt',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_idr_txt': self.qc_dir + '/07.peak_idr.txt',
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'peak_overlap_tiff': self.qc_dir + '/08.peak_overlap.tiff',
            'bam_fingerprint': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key

        if create_dirs:
            check_path([
                self.config_dir,
                self.align_dir,
                self.bam_dir, 
                self.bw_dir, 
                self.bg_dir,
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


class AtacR1Config(object):
    """
    Prepare directories for R1 single replicate

    Required
      - fq1, str
      - fq2, str
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        """
        required arguments for ATACseq analysis
        """
        args_init = {
            'smp_name': None,
            'group': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': 0,
            'keep_tmp': None,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_r1' #

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads, 
            self.parallel_jobs)


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (required)
        """
        # fq1
        if not isinstance(self.fq1, str):
            raise ValueError('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))

        # fq2
        if not isinstance(self.fq2, str):
            raise ValueError('--fq2, str expected, got {}'.format(
                type(self.fq2).__name__))

        if not fq_paired(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not paired properly: {}, {}'.format(
                self.fq1, self.fq2))            

        # abs
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)


    def init_index(self):
        """
        alignment index
        genome, extra_index, genome_index, spikein, spikein_index

        output: genome_index | spikein_index
        """
        # check genome size
        if isinstance(self.extra_index, str):
            self.genome_index = self.extra_index

        ## check genome index
        if isinstance(self.genome_index, str):
            ai = AlignIndex(index=self.genome_index)

            if self.genome_size < 1:
                self.genome_size = ai.index_size()

            if not isinstance(self.genome_size_file, str): 
                self.genome_size_file = ai.index_size(return_file=True)

        elif isinstance(self.genome, str):
            self.genome_index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group='genome')

            self.genome_size_file = Genome(genome=self.genome).get_fasize()

            with open(self.genome_size_file, 'rt') as r:
                s = [i.strip().split('\t')[1] for i in r.readlines()]
            
            if self.genome_size < 1:
                self.genome_size = sum(map(int, s))

        else:
            raise ValueError('index failed, extra_index, genome, genome_index')

        # check spikein index
        if isinstance(self.spikein_index, str):
            ai = AlignIndex(index=self.spikein_index, aligner=self.aligner)

            if not ai.is_index():
                raise ValueError('spikein_index failed, {}'.format(
                    self.spikein_index))

        elif isinstance(self.spikein, str):
            self.spikein_index = AlignIndex(aligner=self.aligner).search(
                genome=self.spikein, group='genome')

        else:
            self.spikein_index = None


    def init_files(self, create_dirs=True):
        """
        for single fastq file
        prepare directories
        """
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        default_dirs = {
            'raw_dir': 'raw_data',
            'clean_dir': 'clean_data',
            'align_dir': 'align',
            'spikein_dir': 'spikein',
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
            'report_log': self.report_dir + '/report.log',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',

            'align_scale_txt': self.bam_dir + '/' + 'scale.txt',
            'trim_stat_txt': self.clean_dir + '/' + self.project_name + '/' + self.project_name + '.trim.stat',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.bowtie2.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.bowtie2.json',
            'spikein_bam': self.spikein_dir + '/' + 'spikein.bam',
            'spikein_scale_txt': self.spikein_dir + '/' + 'scale.txt',
            'spikein_flagstat': self.spikein_dir + '/' + 'spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/' + 'spikein.align.stat',
            'spikein_json': self.spikein_dir + '/' + 'spikein.align.json',

            'trim_summary_json': self.qc_dir + '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_txt': self.qc_dir + '/03.FRiP.txt',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'bam_fingerprint': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key

        # raw data
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            self.raw_fq_list.append(self.raw_dir + '/' + \
                os.path.basename(self.fq2))

        ## clean data
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            self.clean_fq_list.append(self.clean_dir + '/' + \
                os.path.basename(self.fq2))

        if create_dirs:
            check_path([
                self.config_dir, 
                self.raw_dir,
                self.clean_dir,
                self.align_dir,
                self.spikein_dir,
                self.bam_dir, 
                self.bg_dir,
                self.bw_dir, 
                self.peak_dir, 
                self.motif_dir,
                self.qc_dir, 
                self.report_dir])


class AtacRdConfig(object):
    """
    Generate ATACseq design.json

    output: fq_dir/fq1,fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


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


class AtacRpConfig(object):
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
        self.report_log = os.path.join(self.report_dir, 'report.log')

        if create_dirs:
            check_path(self.report_dir)


class AtacReader(object):
    """
    Read config.txt/pickle
    """
    def __init__(self, x):
        self.x = x
        self.read()
        self.atacseq_type = self.args.get('atacseq_type', None)

        self.is_atacseq_r1 = self.atacseq_type == 'atacseq_r1'
        self.is_atacseq_rn = self.atacseq_type == 'atacseq_rn'
        self.is_atacseq_rx = self.atacseq_type == 'atacseq_rx'
        self.is_atacseq_rt = self.atacseq_type == 'atacseq_rt'
        self.is_atacseq_rd = self.atacseq_type == 'build_design'


        if isinstance(self.atacseq_type, str):
            self.is_hiseq = self.atacseq_type.startswith('atacseq_')
        else:
            self.is_hiseq = False


    def read(self):
        """
        # hiseq
        hiseq
          |-config
          |   |-config.pickle

        # alignment
        align_dir
          |- smp_nmae
          |    |- index
          |    |    |- config.pickle
        """
        p1x = os.path.join(self.x, 'config', 'config.pickle')
        p2x = os.path.join(self.x, '*', '*', 'config.pickle')
        p1 = glob.glob(p1x)
        p2 = glob.glob(p2x)

        # read config
        if len(p1) == 1:
            self.args = pickle2dict(p1[0])
        elif len(p2) == 1:
            self.args = pickle2dict(p2[0])
        else:
            self.args = {}


class Align(object):
    """
    Alignment reads for CUT&RUN, CUT&TAG
    suppose paired end reads
    --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700
    
    Input: fq, genome_index, outputdir, stat/log
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_index()
        self.init_files()        
        self.prep_cmd()


    def init_args(self):
        """
        ATACseq, max_fragment = 2000
        CUT&RUN, max_fragment = 700
        CUT&TAG, max_fragment = 700
        """
        default_args = {
            'fq1': None,
            'fq2': None,
            'index': None,
            'genome': None,
            'extra_index': None,
            'outdir': None,
            'aligner': 'bowtie2',
            'threads': 1,
            'overwrite': False,
            'max_fragment': 700,
            'remove_unmap': False,
            'hiseq_type': 'hiseq_alignment',
            'keep_tmp': False
            }
        self = update_obj(self, default_args, force=False)

        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
            log.warning('set pwd as outdir: {}'.format(self.outdir))
        self.outdir = file_abspath(self.outdir)
        # check_path(self.outdir)

        # fq1
        if isinstance(self.fq1, str):
            if check_file(self.fq1):
                self.fq1 = file_abspath(self.fq1)
                self.fq_check = 0
            else:
                log.error('--fq1, file not exists, {}'.format(self.fq1))
                self.fq_check = 1
        else:
            log.error('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))
            self.fq_check = 1

        # fq2
        if isinstance(self.fq2, str):
            if check_file(self.fq2):
                self.fq2 = file_abspath(self.fq2)
                self.fq_check = 0
            else:
                log.error('--fq2, file not exists, {}'.format(self.fq1))
                self.fq_check = 1
        else:
            log.error('--fq2, str expected, got {}'.format(
                type(self.fq2).__name__))
            self.fq_check = 1

        # prefix
        self.smp_name = fq_name(self.fq1, pe_fix=True)


    def init_index(self):
        """
        Determine the index
        """
        if isinstance(self.extra_index, str):
            self.index = self.extra_index
            self.index_check = 0 if AlignIndex(
                index=self.extra_index, aligner=self.aligner).is_index() else 1
        elif isinstance(self.index, str):
            self.index_check = 0 if AlignIndex(
                index=self.index, aligner=self.aligner).is_index() else 1
        elif isinstance(self.genome, str):
            self.index = AlignIndex(aligner=self.aligner).search(
                genome=self.genome, group = 'genome')
            self.index_check = 1 if self.index is None else 0
        else:
            log.error('-g, --index, failed')
            self.index_check = 1

        # index name
        if isinstance(self.index, str):
            self.index_name = AlignIndex(index=self.index).index_name()
        else:
            self.index_name = None


    def init_files(self):
        """
        default files for output
        bam, sam, log, unmap, ...
        """
        # subdir
        subdir = os.path.join(self.outdir, self.smp_name, self.index_name)
        check_path(subdir)

        # output files
        prefix = os.path.join(subdir, self.smp_name)
        default_files = {
            'subdir': subdir,
            'config_pickle': os.path.join(subdir, 'config.pickle'),
            'cmd_shell': os.path.join(subdir, 'cmd.txt'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.fastq',
            'unmap1': prefix + '.unmap.1.fastq',
            'unmap2': prefix + '.unmap.2.fastq',
            'align_log': prefix + '.bowtie2.log',
            'align_stat': prefix + '.bowtie2.stat',
            'align_json': prefix + '.bowtie2.json',
            'align_flagstat': prefix + '.flagstat',
            'rmdup_bam': prefix + '.rmdup.bam',
            'rmdup_matrix': prefix + '.rmdup.matrix.txt'
        }
        self = update_obj(self, default_files, force=True)


    def prep_cmd(self):
        """
        Alignment for CUT&RUN, CUT&TAG, ATACseq, ...
        
        samtools view -bhS -f 2 -F 1804

        include: -f
            2    0x2    Read mapped in proper pair

        exclude: -F
            4    0x4    Read unmapped, 
            8    0x8    Mate unmapped,
            256  0x100  Not primary alignment, 
            512  0x200  Reads fails platform quality checks,
            1024 0x400  Read is PCR or optical duplicate

        """
        self.cmd = ' '.join([
            '{}'.format(shutil.which('bowtie2')),
            '--mm -p {}'.format(self.threads),
            '--local --very-sensitive --no-mixed --no-discordant',
            '-I 10',
            '-X {}'.format(self.max_fragment),
            '-x {}'.format(self.index),
            '-1 {} -2 {}'.format(self.fq1, self.fq2),
            '--un-conc {}'.format(self.unmap),
            '1> {} 2> {}'.format(self.sam, self.align_log),
            '&& samtools view -@ {} -bS'.format(self.threads),
            '<(samtools view -H {} ;'.format(self.sam),
            "samtools view -q 30 -f 2 -F 1804 {} | grep 'YT:Z:CP')".format(self.sam),
            '| samtools sort -@ {} -o {} -'.format(self.threads, self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(self.bam, self.align_flagstat),
            '&& [[ -f {} ]] && rm {}'.format(self.sam, self.sam)
            ])
        
        # save cmd txt
        cmd_txt = self.subdir + '/cmd.txt'
        with open(cmd_txt, 'wt') as w:
            w.write(self.cmd + '\n')


    def parse_align(self, to_json=None):
        """
        Parsing the log of bowtie2
        """
        """
        Wrapper bowtie2 log
        Bowtie2:

        SE:
        10000 reads; of these:
          10000 (100.00%) were unpaired; of these:
            166 (1.66%) aligned 0 times
            2815 (28.15%) aligned exactly 1 time
            7019 (70.19%) aligned >1 times
        98.34% overall alignment rate

        PE:
        100000 reads; of these:
          100000 (100.00%) were paired; of these:
            92926 (92.93%) aligned concordantly 0 times
            5893 (5.89%) aligned concordantly exactly 1 time
            1181 (1.18%) aligned concordantly >1 times
            ----
            92926 pairs aligned concordantly 0 times; of these:
              1087 (1.17%) aligned discordantly 1 time
            ----
            91839 pairs aligned 0 times concordantly or discordantly; of these:
              183678 mates make up the pairs; of these:
                183215 (99.75%) aligned 0 times
                101 (0.05%) aligned exactly 1 time
                362 (0.20%) aligned >1 times
        8.39% overall alignment rate

        unique, multiple, unmap, map, total
        """
        self.unique_only = False

        dd = {}
        se_tag = 1 #
        with open(self.align_log, 'rt') as r:
            for line in r:
                value = line.strip().split(' ')[0]
                if '%' in value:
                    continue
                if line.strip().startswith('----'):
                    continue
                value = int(value)

                ## paired tag
                if 'were paired; of these' in line:
                    dd['total'] = value
                    se_tag = 0
                elif 'aligned concordantly 0 times' in line:
                    dd['unmap'] = value
                    se_tag = 0
                elif 'aligned concordantly exactly 1 time' in line:
                    dd['unique'] = value
                    se_tag = 0
                elif 'aligned concordantly >1 times' in line:
                    dd['multiple'] = value
                    se_tag = 0
                elif 'reads; of these' in line and se_tag:
                    dd['total'] = value
                elif 'aligned 0 times' in line and se_tag:
                    dd['unmap'] = value
                elif 'aligned exactly 1 time' in line and se_tag:
                    dd['unique'] = value
                elif 'aligned >1 times' in line and se_tag:
                    dd['multiple'] = value
                else:
                    pass

        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname,
        dd['fqname'] = self.smp_name
        dd['index_name'] = self.index_name
        self.log_dict = dd

        # save dict to plaintext file
        with open(self.align_stat, 'wt') as w:
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 'fqname', 
                'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        ## save to json
        if to_json:
            Json(dd).writer(self.align_json)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        if self.fq_check > 0:
            raise ValueError('--fq1, --fq2, illegal')
        if self.index_check > 0:
            raise ValueError('--index, --genome, illegal')

        # prepare dirs
        check_path(self.subdir)

        # save config
        dict2pickle(self.__dict__, self.config_pickle)

        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists:{}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))

        # read log
        if file_exists(self.align_log):
            self.parse_align(to_json=True)

        # temp files
        del_list = [self.sam, self.unmap1, self.unmap2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


class CallPeak(object):
    """
    Call peaks using MACS2, or SEACR

    MACS2: 
    macs2 callpeak -t ip.bam -c input.bam \
      -g hs -f BAMPE -n macs2_peak_q0.1 \
      --outdir outdir -q 0.1 
      --keep-dup all 
      2>log.txt
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        default_args = {
            'method': 'macs2',
            'ip': None,
            'input': None,
            'genome': None,
            'genome_size': None,
            'genome_size_file': None,
            'prefix': None,
            'outdir': None,
            'overwrite': False,
            'hiseq_type': 'macs2_call_peak',
            'keep_tmp': False
            }
        self = update_obj(self, default_args, force=False)

        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
            log.warning('set pwd as outdir: {}'.format(self.outdir))
        self.outdir = file_abspath(self.outdir)
        check_path(self.outdir)

        # bam files
        if not file_exists(self.ip):
            raise ValueError('ip, required')

        self.ip = file_abspath(self.ip)
        self.input = file_abspath(self.input)

        # file name
        self.ip_name = os.path.splitext(os.path.basename(self.ip))[0]
        if self.input is None:
            self.input_name = ''
        else:
            self.input_name = os.path.splitext(os.path.basename(self.input))[0]
        
        # prefix
        if not isinstance(self.prefix, str):
            self.prefix = self.ip_name

        self.prefix_top001 = self.prefix + '.top0.01'


    def init_genome(self):
        """
        Determine the genome_size, genome_size_file
        """
        # genome
        g = {
          'dm6': 'dm',
          'dm3': 'dm',
          'mm9': 'mm',
          'mm10': 'mm',
          'hg19': 'hs',
          'hg38': 'hs'
        }

        if isinstance(self.genome, str):
            self.genome_size = g.get(self.genome, 0)
            self.genome_size_file = Genome(genome=self.genome).get_fasize()

            if self.genome_size is None:
                raise ValueError('unknown genome: {}'.format(self.genome))

        elif isinstance(self.genome_size, int):
            if self.genome_size < 1:
                raise ValueError('genome_size, failed, {}'.format(
                    self.genome_size))
            elif not isinstance(self.genome_size_file, str):
                raise ValueError('genome_size_file, failed, {}'.format(
                    self.genome_size_file))
            else:
                pass
        
        else:
            raise ValueError('--genome, --genome-size, required')


    def init_files(self):
        # files
        default_files = {
            'config_pickle': self.outdir + '/config.pickle',
            'macs2_peak': self.outdir + '/' + self.prefix + '_peaks.narrowPeak',
            'macs2_log': self.outdir + '/' + self.prefix + '.callpeak.log',
            'macs2_peak_xls': self.outdir + '/' + self.prefix + '_peaks.xls',
            'ip_bed': self.outdir + '/' + self.ip_name + '.frag.bed',
            'ip_bg': self.outdir + '/' + self.ip_name + '.frag.bg',
            'input_bed': self.outdir + '/' + self.input_name + '.frag.bed',
            'input_bg': self.outdir + '/' + self.input_name + '.frag.bg',
            'seacr_peak': self.outdir + '/' + self.prefix + '.stringent.bed',
            'seacr_peak_top': self.outdir + '/' + self.prefix_top001 + '.stringent.bed'
        }
        self = update_obj(self, default_files, force=True) # key

        ## genome_size
        # self.genome_size_file = Genome(genome=self.genome).get_fasize()
        

    def run_macs2(self):
        input_arg = '' if self.input is None else '-c {}'.format(self.input)

        cmd = ' '.join([
            '{} callpeak'.format(shutil.which('macs2')),
            '-t {} {}'.format(self.ip, input_arg),
            '-g {} -f BAMPE'.format(self.genome_size),
            '-n {}'.format(self.outdir + '/' + self.prefix),
            '-q 0.1 --keep-dup all',
            '--nomodel --shift -100 --extsize 200',
            '2> {}'.format(self.macs2_log)
            ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.macs2.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(self.macs2_peak) and not self.overwrite:
            log.info('run_macs2() skipped, file exists:{}'.format(
                self.macs2_peak))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_macs2() failed, check {}'.format(
                    self.macs2_log))


    def bampe_to_bg(self, bam, bg):
        """
        Convert bam to bed
        bedpe, sort -n (by name)
        """
        bam_sorted = os.path.splitext(bg)[0] + '.sorted_by_name.bam'
        bed = os.path.splitext(bg)[0] + '.bed'

        cmd = ' '.join([
            'samtools sort -@ 4 -n -o {} {}'.format(bam_sorted, bam),
            '&& bedtools bamtobed -bedpe -i {}'.format(bam_sorted),
            "| awk '$1==$4 && $6-$2 < 1000 {print $0}'",
            '| cut -f 1,2,6',
            '| sort -k1,1 -k2,2n > {}'.format(bed),
            '&& bedtools genomecov -bg -i {}'.format(bed),
            '-g {}'.format(self.genome_size_file),
            '> {}'.format(bg)
        ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.bam2bg.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(bg) and not self.overwrite:
            log.info('bampe_to_bg() skipped, file exists:{}'.format(bg))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('bampe_to_bg() failed, check {}'.format(bg))

        # temp files
        del_list = [bam_sorted, bed]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


    def run_seacr(self):
        """
        bash $seacr ip.bedgraph input.bedgraph \
        non stringent output
    
        see: https://github.com/FredHutch/SEACR
        1. bam-to-bed
        bedtools bamtobed -bedpe -i in.bam > out.bed
        awk '$1==$4 && $6-$2 < 1000 {print $0}' out.bed > out.clean.bed
        cut -f1,2,6 out.clean.bed | sort -k1,1 -k2,2n -k3,3n > out.fragments.bed
        bedtools genomecov -bg -i out.fragments.bed -g genome > out.fragments.bg
        
        2. call peak
        bash seacr out.fragments.bg IgG.bg non stringent output
        bash seacr out.fragments.bg 0.01 non 0.01 non stringent top0.01.peaks
        """
        self.bampe_to_bg(self.ip, self.ip_bg)
        cmd = ' '.join([
            'bash {}'.format(shutil.which('SEACR_1.3.sh')),
            '{} 0.01'.format(self.ip_bg),
            'non stringent {}'.format(self.outdir + '/' + self.prefix_top001)
        ])

        if not self.input is None:
            self.bampe_to_bg(self.input, self.input_bg)

            cmd = ' '.join([
                '{} &&'.format(cmd),
                'bash {}'.format(shutil.which('SEACR_1.3.sh')),
                '{} {}'.format(self.ip_bg, self.input_bg),
                'non stringent {}'.format(self.outdir + '/' + self.prefix)            
            ])

        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.SEACR.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')

        # run
        if file_exists(self.seacr_peak_top) and not self.overwrite:
            log.info('run_seacr() skipped, file exists:{}'.format(
                self.seacr_peak))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_seacr() failed, check {}'.format(
                    self.seacr_peak))


    def run(self):
        # save config
        dict2pickle(self.__dict__, self.config_pickle)

        # generate cmd
        if self.method.lower() == 'macs2':
            self.run_macs2()
        elif self.method.lower() == 'seacr':
            self.run_seacr()
        else:
            raise ValueError('unknown method: {macs2|seacr}, got {}'.format(
                self.method))

        # remove files
        del_list = []


def dict2pickle(d, p, update=False):
    """Save dict as pickle
    d is dict
    p is pickle file
    """
    if isinstance(d, dict) and isinstance(p, str):
        if file_exists(p):
            with open(p, 'rb') as r:
                d_old = pickle.load(r)
        else:
            d_old = {}

        # update
        d_old.update(d)
        d_new = d_old if update else d

        # save to file
        with open(p, 'wb') as w:
            pickle.dump(d_new, w, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        log.error('dict2pickle() failed, dict, str expect')


def pickle2dict(p):
    """
    Read pickle, save as dict
    """
    d = {}
    if file_exists(p):
        with open(p, 'rb') as r:
            d = pickle.load(r)
    else:
        log.error('pickle file illegal')

    return d

