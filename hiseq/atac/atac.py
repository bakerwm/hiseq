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
from hiseq.utils.helper import *
from hiseq.qc.trimmer import Trimmer
from hiseq.align.alignment import Alignment
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
        self.reportdir = os.path.join(self.outdir, 'report')
        self.out_prefix = os.path.join(self.outdir, fqname)

        ## names
        fqnames = list(map(os.path.basename, [self.fq1, self.fq2]))
        self.raw_fq_list = [os.path.join(self.rawdir, i) for i in fqnames]
        self.clean_fq_list = [os.path.join(self.cleandir, i) for i in fqnames]
        self.trim_stat = os.path.join(self.cleandir, fqname + '.qc.stat')
        self.bam_raw = os.path.join(self.aligndir, fqname, '2.*', fqname + '.bam')
        self.align_stat = os.path.join(self.aligndir, fqname + '.align.txt')
        self.bam_rmdup = os.path.join(self.bamdir, fqname + '.bam')
        self.peak = os.path.join(self.peakdir, fqname + '_peaks.narrowPeak')
        self.bw = os.path.join(self.bwdir, fqname + '.bigWig')

        ## optional
        args['fq1'] = self.fq1
        args['fq2'] = self.fq2
        args['genome'] = self.genome
        args['outdir'] = self.outdir
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)

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
        args['peak'] = self.peak
        args['bw'] = self.bw

        ## create directories
        tmp = [is_path(d) for d in [self.rawdir, self.cleandir, self.aligndir, 
               self.bamdir, self.bwdir, self.peakdir, self.motifdir, 
               self.reportdir]]

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


    def symlink(self, src, dest):
        """
        Create symlinks within output dir

        ../src
        """
        srcname = os.path.join('..', os.path.basename(src))
        os.symlink(srcname, dest)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        # args = self.args.copy()
        fq1, fq2 = self.config.raw_fq_list
        outdir = self.config.cleandir
        args = self.config.args # all

        if trimmed is True:
            # update clean fastq
            Trimmer(fq1, outdir, fq2, **args).run()
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
        # arguments
        args = self.config.args
        args['fq1'] = args['fq'] = fq1
        args['fq2'] = fq2
        args['outdir'] = self.config.aligndir
        Alignment(**args).run()


    def run(self):
        """
        Run all steps for ATACseq pipeline
        """
        args = self.args.copy() # copy

        # init dir

        # trim
        if args['trimmed'] is True:
            clean_fq1, clean_fq2 = (self.atac[''])
        else:
            args_trim = args_atac(args, trim=True)
            clean_fq1, clean_fq2 = trimmer.Trimmer(
                args['fq1'], args[''])
            # hiseq.qc.trimmer.Trimmer(fq1, outdir).run()

        # align

        # rmdup 

        # peak 

        # motif 

        # qc

        # report







    # def config(self):
    #     """
    #     Create directories, filenames
    #     raw_data
    #     clean_data
    #     align
    #     bam_files
    #     bw_files
    #     peak
    #     motif
    #     report
    #     ...
    #     files:
    #     raw_fastq 
    #     clean fastq
    #     align.bam 
    #     rmdup bam
    #     peak (narrowPeak)
    #     motif (to-do)
    #     report
    #     """
    #     args = self.kwargs.copy()

    #     assert os.path.exists(self.fq1)
    #     assert os.path.exists(self.fq2)
    #     assert is_path(self.outdir)

    #     ## sample name
    #     fqname = file_prefix(self.fq1)[0]
    #     fqname = re.sub('[._][rR]?1$', '', fqname)
    #     if args['smp_name']:
    #         fqname = args['smp_name']

    #     ## outdir
    #     self.rawdir = os.path.join(self.outdir, 'raw_data')
    #     self.cleandir = os.path.join(self.outdir, 'clean_data')
    #     self.aligndir = os.path.join(self.outdir, 'align')
    #     self.bamdir = os.path.join(self.outdir, 'bam_files')
    #     self.bwdir = os.path.join(self.outdir, 'bw_files')
    #     self.peakdir = os.path.join(self.outdir, 'peak')
    #     self.motifdir = os.path.join(self.outdir, 'motif')
    #     self.reportdir = os.path.join(self.outdir, 'report')

    #     ## global variables:
    #     #############################
    #     ## prefix
    #     self.out_prefix = os.path.join(self.outdir, fqname)
    #     ## raw data
    #     fq_names = list(map(os.path.basename, [self.fq1, self.fq2]))
    #     self.raw_fq_list = [os.path.join(self.rawdir, i) for i in fq_names]
    #     ## clean data
    #     self.clean_fq_list = [os.path.join(self.cleandir, i) for i in fq_names]
    #     self.trim_stat = os.path.join(self.cleandir, fqname + '.qc.stat')
    #     ## align data
    #     self.bam_raw = os.path.join(self.aligndir, fqname, '2.*', fqname + '.bam')
    #     self.align_stat = os.path.join(self.aligndir, fqname + '.align.txt')
    #     ## rmdup
    #     self.bam_rmdup = os.path.join(self.bamdir, fqname + '.bam')
    #     ## peak
    #     self.peak = os.path.join(self.peakdir, fqname + '_peaks.narrowPeak')
    #     ## bw file
    #     self.bw = os.path.join(self.bwdir, fqname + '.bigWig')

    #     ## update args
    #     args['rawdir'] = self.rawdir
    #     args['cleandir'] = self.cleandir
    #     args['aligndir'] = self.aligndir
    #     args['bamdir'] = self.bamdir
    #     args['bwdir'] = self.bwdir
    #     args['peakdir'] = self.peakdir
    #     args['motifdir'] = self.motifdir
    #     args['reportdir'] = self.reportdir
    #     args['out_prefix'] = self.out_prefix 
    #     args['raw_fq_list'] = self.raw_fq_list
    #     args['clean_fq_list'] = self.clean_fq_list
    #     args['trim_stat'] = self.trim_stat
    #     args['bam_raw'] = self.bam_raw
    #     args['align_stat'] = self.align_stat
    #     args['bam_rmdup'] = self.bam_rmdup
    #     args['peak'] = self.peak
    #     args['bw'] = self.bw

    #     ## create directories
    #     tmp = map(is_path, [self.rawdir, self.cleandir, self.aligndir, self.bamdir,
    #                   self.bwdir, self.peakdir, self.motifdir, self.reportdir])

    #     return args