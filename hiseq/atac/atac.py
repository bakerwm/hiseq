"""
ATAC-seq pipeline

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
"""

import os
import sys
import re
import pandas as pd
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trimmer
from hiseq.align.alignment import Alignment
from hiseq.peak.call_peak import Macs2
from hiseq.utils.rep_cor import *
import hiseq


class AtacConfig(object):

    def __init__(self, fq1, fq2, genome, outdir, **kwargs):
        """
        Config for atac-seq pipeline
        directories
        files
        ...
        """
        # positional
        self.fq1 = fq1
        self.fq2 = fq2
        self.genome = genome
        self.outdir = outdir
        self.kwargs = kwargs

        # update
        self.args = self.args_init()
        self.check_status = self.check() # check


    def check(self):
        """
        Check arguments, target file exists
        """
        args = self.args.copy()

        ## check args:
        config_dir = os.path.join(self.outdir, 'config')
        assert is_path(config_dir)
        args_pickle = os.path.join(config_dir, 'arguments.pickle') # out_prefix
        args_file = os.path.join(config_dir, 'arguments.txt')
        args_json = os.path.join(config_dir, 'arguments.json')

        ## check files:
        chk1 = args_checker(args, args_pickle)
        chk2 = args['overwrite'] is False

        args_logger(args, args_file, True) # update arguments.txt
        Json(args).writer(args_json) # save to json

        return all([chk1, chk2]) # arguments.pickle, changed 


    def args_init(self):
        """
        Init arguments for ATAC-seq pipeline
        ##
        Create directories, filenames
        raw_data
        clean_data
        align
        bam_files
        bw_files
        peak
        motif
        report
        ...
        files:
        raw_fastq 
        clean fastq
        align.bam 
        rmdup bam
        peak (narrowPeak)
        motif (to-do)
        report
        """
        args = self.kwargs.copy()

        ## sample name
        args['smp_name'] = args.get('smp_name', None)
        fqname = file_prefix(self.fq1)[0]
        fqname = re.sub('[._][rR]?1$', '', fqname)
        if args['smp_name']:
            fqname = args['smp_name']

        self.fqname = fqname
        ## outdir
        self.rawdir = os.path.join(self.outdir, 'raw_data')
        self.cleandir = os.path.join(self.outdir, 'clean_data')
        self.aligndir = os.path.join(self.outdir, 'align')
        self.bamdir = os.path.join(self.outdir, 'bam_files')
        self.bwdir = os.path.join(self.outdir, 'bw_files')
        self.peakdir = os.path.join(self.outdir, 'peak')
        self.motifdir = os.path.join(self.outdir, 'motif')
        self.qcdir = os.path.join(self.outdir, 'qc')
        self.reportdir = os.path.join(self.outdir, 'report')
        self.out_prefix = os.path.join(self.outdir, fqname)

        ## names
        fqnames = list(map(os.path.basename, [self.fq1, self.fq2]))
        self.raw_fq_list = [os.path.join(self.rawdir, i) for i in fqnames]
        self.clean_fq_list = [os.path.join(self.cleandir, i) for i in fqnames]
        self.trim_stat = os.path.join(self.cleandir, fqname + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, fqname, '2.*', fqname + '.bam')
        self.align_stat = os.path.join(self.aligndir, fqname + '.align.txt')
        self.bam_rmdup = os.path.join(self.bamdir, fqname + '.rmdup.bam')
        self.bam_proper_pair = os.path.join(self.bamdir, fqname + '.proper_pair.bam')
        self.peak = os.path.join(self.peakdir, fqname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, fqname + '.bigWig')
        ## qc files
        self.lendist_txt = os.path.join(self.qcdir, 'length_distribution.txt')
        self.lendist_pdf = os.path.join(self.qcdir, 'length_distribution.pdf')
        self.frip_txt = os.path.join(self.qcdir, 'FRiP.txt')

        ## optional
        args['fq1'] = self.fq1
        args['fq2'] = self.fq2
        args['genome'] = self.genome
        args['outdir'] = self.outdir
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)
        args['fqname'] = self.fqname

        ## trimming, cutadapt
        adapter3 = hiseq.utils.args.Adapter('nextera').adapters[0] #
        args['len_min']  = args.get('len_min', 20)
        args['adapter3'] = args.get('adapter3', adapter3) # nextera
        args['AD3'] = args.get('AD3', adapter3) # nextera

        ## alignment
        args['aligner'] = args.get('aligner', 'bowtie2') # bowtie alignment
        args['n_map'] = args.get('n_map', 2)
        args['align_to_chrM'] = True

        ## update args
        args['rawdir'] = self.rawdir
        args['cleandir'] = self.cleandir
        args['aligndir'] = self.aligndir
        args['bamdir'] = self.bamdir
        args['bwdir'] = self.bwdir
        args['peakdir'] = self.peakdir
        args['motifdir'] = self.motifdir
        args['reportdir'] = self.reportdir
        args['out_prefix'] = self.out_prefix 
        args['raw_fq_list'] = self.raw_fq_list
        args['clean_fq_list'] = self.clean_fq_list
        args['trim_stat'] = self.trim_stat
        args['bam_raw'] = self.bam_raw
        args['align_stat'] = self.align_stat
        args['bam_rmdup'] = self.bam_rmdup
        args['bam_proper_pair'] = self.bam_proper_pair
        args['peak'] = self.peak
        args['bw'] = self.bw

        ## create directories
        tmp = [is_path(d) for d in [self.rawdir, self.cleandir, self.aligndir, 
               self.bamdir, self.bwdir, self.peakdir, self.motifdir, 
               self.qcdir, self.reportdir]]

        return args


class Atac(object):
    """
    ATAC-seq analysis pipeline
    replicates: 1,2,3,...
    paired-end:
    """
    def __init__(self, design, **kwargs):
        """
        Parsing design file: Json
        fq1, fq2, genome, outdir
        """
        self.design = design
        self.kwargs = kwargs
        args_design = self.atac_init()
        # required
        self.fq1 = args_design.pop('fq1')
        self.fq2 = args_design.pop('fq2')
        self.genome = args_design.pop('genome')
        self.outdir = args_design.pop('outdir')
        self.outdir = os.path.abspath(self.outdir)
        # final/global config, args
        self.config = AtacConfig(
            self.fq1, 
            self.fq2, 
            self.genome, 
            self.outdir, 
            **args_design)


    def atac_init(self):
        """
        Initiate directories, config 
        """
        args = self.read_design()
        # update args
        args.update(self.kwargs)
        # check required arguments
        # fq1, fq2, genome, outdir
        args_required = ['fq1', 'fq2', 'genome', 'outdir']
        chk = [i for i in args_required if not i in args] # missing args
        if len(chk) > 0:
            raise Exception('check design file for: {}'.format(' '.join(chk)))

        return args


    def read_design(self):
        """
        default parameters for ATAC-seq
        Json -> dict
        """
        return Json(self.design).dict


    def symlink(self, src, dest, absolute_path=True):
        """
        Create symlinks within output dir
        ../src
        """
        if absolute_path:
            # support: ~, $HOME,
            srcname = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        else:
            srcname = os.path.join('..', os.path.basename(src))

        if not os.path.exists(dest):
            os.symlink(srcname, dest)


    #######################################
    ## main pipeline
    def prep_raw(self, copy=False):
        """
        Copy raw data to dest dir
        if not: create a symlink
        
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.config.raw_fq_list

        # copy
        if copy is True:
            shutil.copy(self.fq1, raw_fq1)
            shutil.copy(self.fq2, raw_fq2)
        else:
            self.symlink(self.fq1, raw_fq1, absolute_path=True)
            self.symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.trim.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        # args = self.args.copy()
        fq1, fq2 = self.config.raw_fq_list

        # update args
        args = self.config.args # all
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.cleandir
        if not trimmed is True:
            Trimmer(**args).run()
        else:
            # create links
            for s, d in zip(self.config.raw_fq_list, self.config.clean_fq_list):
                self.symlink(s, d)


    def align(self):
        """
        Alignment

        -X 2000, ...
        """
        fq1, fq2 = self.config.clean_fq_list

        # update arguments
        args = self.config.args
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.aligndir
        #
        args['extra_para'] = '-X 2000'
        Alignment(**args).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bamdir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        # bamdir = os.path.join(self.config.aligndir, '2.*')
        bamdir = self.config.align_stat.rstrip('.align.txt')
        # print(bamdir)
        bamlist = listfiles2('*.bam', bamdir, recursive=True)
        # [chrM, genome]
        return(bamlist[1]) # genome


    def bam_rmdup(self, rmdup=True):
        """
        Remove PCR dup from BAM file using sambamfa/Picard 
        save only proper paired PE reads
        """
        bam_raw = self.get_raw_bam()

        # rmdup
        if rmdup is True:
            Bam(bam_raw).rmdup(self.config.bam_rmdup)
        else:
            shutil.copy(bam_raw, self.config.bam_rmdup)

        # save proper pair reads
        Bam(self.config.bam_rmdup).proper_pair(self.config.bam_proper_pair)
        # index
        Bam(self.config.bam_proper_pair).index()


    def callpeak(self):
        """
        Call peaks using MACS2
        """
        args = self.config.args.copy()
        bam = self.config.bam_proper_pair
        genome = args.pop('genome', None)
        output = args.pop('peakdir', None)
        prefix = args.pop('fqname', None)
        m = Macs2(bam, genome, output, prefix, atac=True, **args).callpeak()


    def qc_lendist(self):
        """
        Create length distribution plot
        """
        # compute length distribution
        x = frag_length(self.config.bam_proper_pair, self.config.lendist_txt)

        # create plot
        pkg_dir = os.path.dirname(hiseq.__file__)
        lendistR = os.path.join(pkg_dir, 'bin', 'qc.lendist.R')
        cmd = 'Rscript {} {} {}'.format(
            lendistR,
            self.config.lendist_txt,
            self.config.lendist_pdf)
        run_shell_cmd(cmd)


    def qc_frip(self):
        """
        Compute frip
        """
        frip, n, total = cal_FRiP(self.config.peak, self.config.bam_proper_pair)

        # head
        hd = ['FRiP', "peak_reads", "total_reads", "id"]
        n = list(map(str, [frip, n, total]))
        n.append('self.config.fqname')
        with open(self.config.frip_txt, 'wt') as w:
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
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report_single.R')
        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.config.outdir,
            self.config.reportdir)
        run_shell_cmd(cmd)


    def run(self):
        """
        Run all steps for ATACseq pipeline
        """
        # init dir
        args = self.config.args.copy()

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
        self.bam_rmdup(rmdup)

        # 5. peak 
        self.callpeak()
        
        # 6. motif

        # 7. qc
        self.qc()

        # 8.report
        self.report()


class AtacBatch(object):
    """
    option-1
    required: fq1, fq2, genome, outdir
    option: len_min, n_map, overwrite, smp_name, threads, adapter3, AD3

    fq1 (list), fq2 (list) => multiple samples

    option-2
    config: table
    fq1, fq2, genome, outdir, ...

    one sample per line
    """

    def __init__(self, **kwargs):
        """
        option-1:
        --config: (one sample per line)

        option-2:
        fq1, fq2, genome, outdir, ...
        """
        args = kwargs # udpate

        # required args
        self.fq1_list = args.pop('fq1', None)
        self.fq2_list = args.pop('fq2', None)
        self.genome = args.pop('genome', None)
        self.outdir = args.pop('outdir', None)
        self.config = args.pop('config', None)
        self.args = args # global

        # json_list
        if not self.config is None:
            json_list = self.config_to_json()
        else:
            json_list = self.args_to_json()

        if len(json_list) == 0:
            raise ValueError('required: --config or --fq1, --fq2, --genome, --outdir')

        # update: self.json_list


    def config_to_json(self):
        """
        convert config to multiple *.json files
        """
        args = self.args.copy()

        if self.config is None:
            return None

        # update arguments
        json_list = []
        nline = 0
        with open(self.config) as r:
            for line in r:
                nline += 1
                if line.startswith('#') or line.strip() == '':
                    continue # comment, blank rows
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    logging.warning('line-{}, arguments not enough: \n{}'.format(nline, line))
                    continue
                # required args
                fq1, fq2, genome, outdir = fields[:4]
                assert is_path(outdir)
                # prefix
                fq_prefix = file_prefix(fq1)[0] # name
                # suboutdir
                suboutdir = os.path.join(outdir, fq_prefix)
                config_json = suboutdir + '.config.json' # file
                config_pickle = suboutdir + '.config.pickle' # file
                d = {'fq1': fq1, 'fq2': fq2, 'genome': genome, 'outdir': suboutdir}
                d.update(args) # global config
                Json(d).writer(config_json)
                self.dict_to_pickle(d, config_pickle)
                json_list.append(config_json)

        self.json_list = json_list
        return json_list


    def args_to_json(self):
        """
        fq1, fq2 list, fastq files
        genome, str
        outdir, str
        ...
        """
        args = self.args.copy()

        # dir
        assert is_path(self.outdir)

        json_list = []
        if isinstance(self.fq1_list, list) and isinstance(self.fq2_list, list):
            # require genome, outdir
            assert isinstance(self.genome, str)
            assert isinstance(self.outdir, str)
            if len(self.fq1_list) == len(self.fq2_list):
                for fq1, fq2 in zip(self.fq1_list, self.fq2_list):
                    # prefix
                    fq_prefix = file_prefix(fq1)[0] # name
                    # suboutdir
                    suboutdir = os.path.join(self.outdir, fq_prefix)
                    config_json = suboutdir + '.config.json' # file
                    config_pickle = suboutdir + '.config.pickle' # file
                    d = {'fq1': fq1, 'fq2': fq2, 'genome': self.genome, 'outdir': suboutdir}
                    d.update(args) # global config
                    Json(d).writer(config_json)
                    self.dict_to_pickle(d, config_pickle)
                    json_list.append(config_json)
            else:
                raise Exception('Error, the --fq1 and --fq2 are not same in length')
        else:
            return None
        
        self.json_list = json_list
        return json_list


    def dict_to_log(self, d, x, overwrite=False):
        """
        Convert dict to log style
            key | value
        """
        assert isinstance(d, dict)
        logout = ['%30s |    %-40s' % (k, d[k]) for k in sorted(d.keys())]
        if overwrite is True or not os.path.exists(x): 
            with open(x, 'wt') as w:
                w.write('\n'.join(logout) + '\n')

        return '\n'.join(logout)


    def dict_to_pickle(self, d, x):
        """
        Convert dict to pickle
        """
        assert isinstance(d, dict)
        with open(x, 'wb') as w:
            pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)


    def pickle_to_dict(self, x):
        """
        Convert pickle file to dict
        """
        with open(x, 'rb') as r:
            return pickle.load(r)


    def run(self):
        """
        Run multiple fastq files
        with input: confi
        """
        args = self.args.copy()

        # iterator json_list
        for config in self.json_list:
            Atac(config, **args).run()


class AtacConfig2(object):
    """
    Config for replicates

    save: merge to the same directories
    merge: report + html 
    """
    def __init__(self, rep_list, **kwargs):
        self.smpname = kwargs.get('smpname', None)
        self.outdir = kwargs.get('outdir', None)
        self.rep_list = rep_list

        # rep_list: config
        self.rep_config = [self.rep_parser(i) for i in rep_list]

        # global vars
        self.align_txt = [c.align_stat for c in self.rep_config]
        self.bam_list = [c.bam_proper_pair for c in self.rep_config]
        self.peak_list = [c.peak for c in self.rep_config]
        self.frip_list = [c.frip_txt for c in self.rep_config]

        # outdir
        self.args_init()
        kwargs.update(self.config_to_dict())
        self.args = kwargs
        self.config_save()


    def args_init(self):
        """
        Init arguments for ATAC-seq pipeline, merge

        Create directories, filenames
        config/
        bam_files/
        bw_files/
        peak/
        motif/
        qc/
        report/
        ...

        proper_pair bam
        bigWig files
        peak (narrowPeak)
        motif (to-do)
        qc/
        report/html
        ...

        qc/
        align.txt
        peaks
        FRiP
        bam_cor
        peak_overlap
        IDR
        """
        # args = self.kwargs.copy()

        ## sample name
        if self.smpname is None:
            self.smpname = self.merge_name()

        ## outdir
        if self.outdir is None:
            # the first one
            self.outdir = os.path.join(
                os.path.dirname(self.rep_config[0].outdir), self.smpname)

        ## outdir
        self.configdir = os.path.join(self.outdir, 'config')
        self.aligndir = os.path.join(self.outdir, 'align')
        self.bamdir = os.path.join(self.outdir, 'bam_files')
        self.bwdir = os.path.join(self.outdir, 'bw_files')
        self.peakdir = os.path.join(self.outdir, 'peak')
        self.motifdir = os.path.join(self.outdir, 'motif')
        self.qcdir = os.path.join(self.outdir, 'qc')
        self.reportdir = os.path.join(self.outdir, 'report')

        ## names
        self.bam = os.path.join(self.bamdir, self.smpname + '.bam')
        self.peak = os.path.join(self.peakdir, self.smpname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, self.smpname + '.bigWig')
        self.lendist_txt = os.path.join(self.qcdir, 'length_distribution.txt')
        self.lendist_pdf = os.path.join(self.qcdir, 'length_distribution.pdf')
        self.frip_txt = os.path.join(self.qcdir, 'FRiP.txt')
        self.cor_npz = os.path.join(self.qcdir, 'cor.bam.npz')
        self.cor_counts = os.path.join(self.qcdir, 'cor.bam.counts.tab')
        self.idr_txt = os.path.join(self.qcdir, 'idr.txt')
        self.idr_log = os.path.join(self.qcdir, 'idr.log')
        self.idr_png = self.idr_txt + '.png'
        self.peak_overlap = os.path.join(self.qcdir, 'peak_overlap.pdf')

        ## create directories
        tmp = [is_path(d) for d in [self.configdir, self.aligndir, self.bamdir, 
            self.bwdir, self.peakdir, self.motifdir, self.qcdir, 
            self.reportdir]]


    def merge_name(self):
        """
        Get the name of replictes
        common in left-most
        """
        name_list = [os.path.basename(i) for i in self.rep_list]
        name_list = [re.sub('.rep[0-9]', '', i) for i in name_list]
        return list(set(name_list))[0]


    def is_abspath(self, x):
        """
        x is the absolute path
        """
        return x == os.path.abspath(x)


    def rep_parser(self, x):
        """
        Parse the config, pickle, for each replicate
        return:
        BAM,
        ...
        """
        design = x + '.config.json'
        args_pickle = x + '.config.pickle'
        args = self.pickle_to_dict(args_pickle)
        return Atac(design, **args).config


    def config_to_dict(self):
        """
        Save config to dict
        """
        args = {
            'rep_list': self.rep_list,
            'bam': self.bam_list,
            'peak': self.peak_list, 
            'idr_txt': self.idr_txt,
            'cor_counts': self.cor_counts,
            'peak_overlap': self.peak_overlap,
        }

        return args


    def config_save(self):
        """
        Save config to json file
        """
        config_json = os.path.join(self.configdir, 'arguments.json')
        config_pickle = os.path.join(self.configdir, 'arguments.pickle')
        config_txt = os.path.join(self.configdir, 'arguments.txt')
        replist_txt = os.path.join(self.configdir, 'sample_list.txt')
        Json(self.args).writer(config_json)
        self.dict_to_pickle(self.args, config_pickle)
        args_logger(self.args, config_txt, True) # update arguments.txt
        # save replist
        with open(replist_txt, 'wt') as w:
            w.write('\n'.join(self.rep_list) + '\n')

                
    def dict_to_pickle(self, d, x):
        """
        Convert dict to pickle
        """
        assert isinstance(d, dict)
        with open(x, 'wb') as w:
            pickle.dump(d, w, protocol=pickle.HIGHEST_PROTOCOL)


    def pickle_to_dict(self, x):
        """
        Convert pickle file to dict
        """
        with open(x, 'rb') as r:
            return pickle.load(r)


class Atac2(object):
    """
    For replicates:
    1. mito%
    2. frag size
    3. TSS enrich
    4. Correlation (pearson)
    5. IDR
    6. Peaks overlap


    directory:
    qc
    report

    report:
    bam correlation
    peaks overlap
    IDR
    """
    def __init__(self, rep_list, **kwargs):
        """
        A list of directories of replicates:
        """
        self.config = AtacConfig2(rep_list, **kwargs)
        self.bam_list = self.config.bam_list # [c.bam_proper_pair for c in self.config.rep_config]
        self.peak_list = self.config.peak_list  # [c.peak for c in self.config.rep_config]


    def mergebam(self):
        """
        Merge replicates, BAM
        """

        cmd = 'samtools merge {} {} && \
            samtools sort -o {} {} && \
            samtools index {}'.format(
            self.config.bam + '.tmp',
            ' '.join(self.bam_list),
            self.config.bam,
            self.config.bam + '.tmp',
            self.config.bam)

        if os.path.exists(self.config.bam):
            logging.info('file exists: {}'.format(self.config.bam))
        else:
            run_shell_cmd(cmd)


    def callpeak(self):
        """
        Call peaks using MACS2
        """        
        bam = self.config.bam
        genome = self.config.rep_config[0].genome
        output = self.config.peakdir
        prefix = self.config.smpname
        m = Macs2(bam, genome, output, prefix, atac=True).callpeak()


    def bam_cor(self, window = 500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp
        
        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz --outRawCounts *counts.tab -b bam
        """
        multiBamSummary = shutil.which('multiBamSummary')

        cmd = '{} bins --binSize {} --smartLabels -o {} --outRawCounts {} -b {}'.format(
            multiBamSummary, 
            window,
            self.config.cor_npz,
            self.config.cor_counts,
            ' '.join(self.config.bam_list))

        if os.path.exists(self.config.cor_counts):
            logging.info('file.exsits: {}'.format(self.config.cor_counts))
        else:
            run_shell_cmd(cmd)


    def peak_overlap(self):
        """
        Compute the overlaps between overlaps
        """
        pkg_dir    = os.path.dirname(hiseq.__file__)
        peak_overlapR = os.path.join(pkg_dir, 'bin', 'atac_peak_overlap.R')

        cmd = 'Rscript {} {} {}'.format(
            peak_overlapR,
            self.config.qcdir,
            ' '.join(self.peak_list))

        if os.path.exists(self.config.peak_overlap):
            logging.info('file exists: {}'.format(self.config.peak_overlap))
        else:
            run_shell_cmd(cmd)


    def rep_idr(self):
        """
        Calculate IDR for replicates
        1 vs 1
        peak files
        """
        idr = shutil.which('idr') # command
        # sort peakfile by p.value
        cmd = 'sort -k8,8nr -o {} {} && \
            sort -k8,8nr -o {} {} && \
            {} --input-file-type narrowPeak --rank p.value --plot --output-file {} \
            --log-output-file {} --samples {}'.format(
            self.peak_list[0],
            self.peak_list[0],
            self.peak_list[1],
            self.peak_list[1],
            idr,
            self.config.idr_txt,
            self.config.idr_log,
            ' '.join(self.peak_list))

        if os.path.exists(self.config.idr_txt):
            logging.info('file exists: {}'.format(self.config.idr_txt))
        else:
            run_shell_cmd(cmd)


    def get_frip(self):
        """
        Save all FRiP.txt file to one
        """
        # print(self.config.frip_txt)
        with open(self.config.frip_txt, 'wt') as w:
            for f in self.config.frip_list:
                with open(f) as r:
                    for line in r:
                        w.write(line)


    def get_align_txt(self):
        """
        Copy align.txt files to align/
        """
        pass


    def tss_enrich(self):
        """
        Calculate the TSS enrichment
        """
        pass


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'atac_report_single.R')
        cmd = 'Rscript {} {} {}'.format(
            qc_reportR,
            self.config.outdir,
            self.config.reportdir)
        run_shell_cmd(cmd)



    def run(self):
        ## question:
        ## merge bam or merge peaks?!
        # merge
        self.mergebam()
        self.callpeak()
        # qc
        self.bam_cor()
        self.peak_overlap()
        self.rep_idr()
        self.get_frip()
        # self.get_align_txt()
        # self.tss_enrich()
        self.report()


class AtacBatch2(object):
    """
    Run merge for all samples
    """
    def __init__(self, sample_dirs, **kwargs):
        """
        option-1:
        --config: list of config json files
                  the same directory with outdir
        """
        self.sample_dirs = [i for i in sample_dirs if self.is_atac_single(i)]
        self.sample_dirs = sorted(self.sample_dirs)
        self.outdir = kwargs.pop('outdir', None)
        self.name_list = self.unique_names(self.sample_dirs)
        self.name_list = sorted(self.name_list)
        self.args = kwargs # global

        if len(sample_dirs) < 2:
            raise ValueError('required: --sample-dirs, >=2 dirs')

        
    def unique_names(self, name_list):
        """
        Get the name of replictes
        common in left-most
        """
        name_list = [os.path.basename(i) for i in name_list]
        name_list = [re.sub('.rep[0-9].*', '', i) for i in name_list]
        return list(set(name_list))


    def is_atac_single(self, x):
        """
        Check input directory is ATAC_single outdir
        required directories
        required files
        
        # config file is 
        x.config.json
        x.config.pickle
        """
        x_config = self.rep_parser(x)

        if x_config is None:
            return None
        else:
            # required directories
            atac_required = [
                x_config.rawdir,
                x_config.cleandir,
                x_config.aligndir,
                x_config.bamdir,
                x_config.peakdir,
                x_config.qcdir,
                x_config.reportdir,
                x_config.frip_txt,
                x_config.lendist_txt,
                x_config.align_stat,
                x_config.peak,
                x_config.bam_proper_pair]
            return all([os.path.exists(i) for i in atac_required])


    def rep_parser(self, x):
        """
        Parse the config, pickle, for each replicate
        return:
        BAM,
        ...
        """
        design = x + '.config.json'
        args_pickle = x + '.config.pickle'
        if os.path.exists(args_pickle) and os.path.exists(design):
            args = self.pickle_to_dict(args_pickle)
            return Atac(design, **args).config
        else:
            return None


    def pickle_to_dict(self, x):
        """
        Convert pickle file to dict
        """
        with open(x, 'rb') as r:
            return pickle.load(r)


    def get_samples(self, x):
        """
        samples from list
        """
        return [i for i in self.sample_dirs if x in os.path.basename(i)]


    def run(self):
        """
        Iterator all sample names
        """
        for x in self.name_list:
            x_list = self.get_samples(x)
            # check outdir
            if self.outdir is None:
                x_outdir = os.path.join(
                    os.path.dirname(x_list[0]),
                    x)
            Atac2(x_list, outdir=x_outdir).run()
            # print(x, x_list)
            # print(x_outdir)
