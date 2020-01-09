
# -*- coding: utf-8 -*-

"""
Default arguments for subcommands
"""

import os
import re
import pathlib
import logging

from hiseq.utils.helper import *

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)



################################################
## BUG: this function could not be imported
## from hiseq.utils.helper import file_prefix
################################################
def file_prefix(fn, with_path=False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz2'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]



class Adapter(object):
    """
    The default adapters, TrueSeq, Nextera, smallRNA,
    TruSeq: nsr, rnaseq, chipseq, 
    Nextera: atacseq, dnaseq,
    smallRNA: smallRNAseq

    ## options for trimming
    nsr: cut 7, -7, -m 20
    chipseq: -m 20 
    iclip: -m 15, rmdup, cut 9, discard-untrimmed, 
    eclip: -m 15, rmdup, cut 10, discard-untrimmed
    atacseq: -m 20
    smrna: -m 18

    """

    def __init__(self, lib='truseq'):
        libtype = {
            'truseq': ['AGATCGGAAGAGC', 'AGATCGGAAGAGC'],
            'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'],
            'smrna': ['TGGAATTCTCGGGTGCCAAGG', '']}

        self.lib = lib
        self.adapters = libtype.get(lib, ['AGATCGGAAGAGC', 'AGATCGGAAGAGC'])
        self.libtype = libtype


    def adapter(self):
        return self.adapters


    def for_trimmer(self):
        ## minimum
        args_mini = {
            'adapter': self.libtype['truseq'],
            'len_min': 20,
            'rmdup': False,
            'cut_after_trim': '0'
        }


        ## default parameters
        args_trimmer = {
            'rnaseq' : {
                'adapter': self.libtype['truseq'],
                'len_min': 20,
                'rmdup': False,
                'cut_after_trim': '7,-7',
            },
            'chipseq': {
                'adapter': self.libtype['truseq'],
                'len_min': 20,
                'rmdup': False,
                'cut_after_trim': '0'
            },
            'iclip': {
                'adapter': self.libtype['truseq'],
                'len_min': 15,
                'rmdup': True,
                'cut_after_rmdup': '9',
                'adapter_sliding': True,
                'trim_times': 4,
                'rm_untrim': True
            },
            'eclip' : {
                'adapter': self.libtype['truseq'],
                'len_min': 15,
                'rmdup': True,
                'cut_after_trim': '-7',
                'cut_after_rmdup': '10',
                'adapter_sliding': True,
                'trim_times': 4,
                'rm_untrim': True
            },
            'clipnsr' : {
                'adapter': self.libtype['truseq'],
                'len_min': 15,
                'rmdup': True,
                'cut_after_trim': '7,-7',
                'adapter_sliding': True,
                'trim_times': 4,
                'rm_untrim': True
            },
            'atacseq': {
                'adapter': self.libtype['truseq'],
                'len_min': 20,
                'rmdup': False,
                'cut_after_trim': '0',
                'adapter_sliding': True,
            },
            'smrna': {
                'adapter': self.libtype['smrna'],
                'len_min': 18,
                'rmdup': True,
                'cut_after_trim': '7,-7',
            }
        }

        ## library-type
        return args_trimmer.get(self.lib, args_mini)



class ArgumentsInit(object):
    """
    Initiate the arguments for hiseq commands: demx, qc, align, ...
    """
    def __init__(self, *args, **kwargs):
        """
        Input: original kwargs; commands
        return the updated one
        # keep keys in kwargs
        """
        self.args_input = args[0]
        self.cmd_input = kwargs

        # save all kwargs to attributes
        for k, v in kwargs.items():
            setattr(self, k, v)

        self.dict = self.config()


    def args_global(self):
        """
        Global arguments for all commands
        fq1, fq2(optional), outdir, genome(optional), ...
        """
        args = self.args_input.copy()

        ## required, to absolute path
        self.fq = args.get('fq', None) # deprecated
        self.fq1 = args.get('fq1', None)
        self.fq2 = args.get('fq2', None)
        self.outdir = args.get('outdir', None)
        
        ## common args
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())

        ## convert
        if isinstance(self.fq1, str):
            self.fq1 = os.path.abspath(self.fq1)
        elif isinstance(self.fq2, list):
            self.fq1 = [os.path.abspath(i) for i in self.fq1]
        else:
            pass
            # log.warning('fq1=[], str and list expect, get {}'.format(
            #     type(self.fq1)))

        if isinstance(self.fq2, str):
            self.fq2 = os.path.abspath(self.fq2)
        elif isinstance(self.fq2, list):
            self.fq2 = [os.path.abspath(i) for i in self.fq2]
        else:
            pass
            # log.warning('fq2=[], str and list expect, get {}'.format(
            #     type(self.fq2)))

        self.outdir = os.path.abspath(self.outdir)

        ## optional
        self.genome_path = args.get('genome_path', 
                os.path.join(str(pathlib.Path.home()), 'data', 'genome'))            
        self.overwrite = args.get('overwrite', False)
        self.threads = args.get('threads', 8)
        self.smp_name = args.get('smp_name', None)


    def args_demx(self):
        """
        For demultiplex Illumina one lane
        using inline-barcode + p7-index
        """
        args = self.args_input.copy()
        self.demx_type = args.get('demx_type', 'p7') # p7, barcode, both
        self.n_mismatch = args.get('n_mismatch', 0)
        self.bc_n_left = args.get('bc_n_left', 0)
        self.bc_n_right = args.get('bc_n_right', 0)
        self.bc_in_read = args.get('bc_in_read', 1) # 1, 2
        self.cut = args.get('cut', False) # True, False


    def args_trim(self):
        args = self.args_input.copy()
        self.len_min = args.get('len_min', 15)
        self.adapter3 = args.get('adapter3', Adapter().adapters[0]) # TruSeq
        self.adapter5 = args.get('adapter5', None)
        self.qual_min = args.get('qual_min', 20)
        self.error_rate = args.get('error_rate', 0.1)
        self.overlap = args.get('overlap', 3)
        self.percent = args.get('percent', 80)

        self.keep_name = args.get('keep_name', True)
        self.save_temp = args.get('save_temp', False)
        self.gzip = args.get('gzip', True)
        self.trim_times = args.get('trim_times', 1)
        self.rm_untrim = args.get('rm_untrim', False)
        self.save_untrim = args.get('save_untrim', False)
        self.save_too_short = args.get('save_too_short', False)
        self.save_too_long = args.get('save_too_long', False)

        self.adapter_sliding = args.get('adapter_sliding', False)
        self.cut_before_trim = args.get('cut_before_trim', '0')
        self.cut_after_trim = args.get('cut_after_trim', '0')
        self.cut_after_rmdup = args.get('cut_after_rmdup', '0')
        self.cut_to_length = args.get('cut_to_length', 0)
        self.rmdup = args.get('rmdup', False)
        self.not_trim_adapter = args.get('not_trim_adapter', False)
        self.keep_temp_files = args.get('keep_temp_files', False)

        ## PE reads
        self.AD3 = args.get('AD3', Adapter().adapters[0])
        self.AD5 = args.get('AD5', None)


    def args_align(self):
        """
        For alignment
        """
        args = self.args_input.copy()

        self.aligner = args.get('aligner', 'bowtie')
        self.genome = args.get('genome', None)
        self.spikein = args.get('spikein', None)
        self.index = args.get('index', None)
        self.index_list = args.get('index_list', None)
        self.index_name = args.get('index_name', None)
        self.extra_index = args.get('extra_index', None)
        self.unique_only = args.get('unique_oinly', True)        
        self.align_by_order = args.get('align_by_order', True)
        self.n_map = args.get('n_map', 0)

        self.search_index = args.get('search_index', False)

        self.align_to_chrM = args.get('align_to_chrM', False)
        self.align_to_rRNA = args.get('align_to_rRNA', False)
        self.align_to_MT_trRNA = args.get('align_to_MT_trRNA', False)
        self.repeat_masked_genome = args.get('repeat_masked_genome', False)

        self.small_genome = args.get('small_genome', False)
        self.simplify_name = args.get('simplify_name', True)
        self.index_parallel = args.get('index_parallel', False)
        self.extra_para = args.get('extra_para', None)

        # check spikein
        if self.spikein == self.genome:
            self.spikein = None


    def args_peak(self):
        """
        Arguments for peak calling
        """
        args = self.args_input.copy()

        self.peak_caller = args.get('peak_caller', 'Macs2')
        self.threshold = args.get('threshold', 1)
        self.intersect = args.get('intersect', 0)


    def args_bam2bw(self):
        """
        Arguments for Bam to BigWig conversion
        """
        args = self.args_input.copy()

        self.filterRNAstrand = args.get('filterRNAstrand', None)
        self.samFlagExclude = args.get('samFlagExclude', None)
        self.samFlagInclude = args.get('samFlagInclude', None)
        self.binsize = args.get('binsize', 10)


    def check_cmd(self):
        """
        Which command to use
        """
        cmd = self.cmd_input

        if cmd.get('demx', False):
            self.args_demx()

        if cmd.get('trim', False):
            self.args_trim()

        if cmd.get('align', False):
            self.args_align()

        if cmd.get('peak', False):
            self.args_peak()

        if cmd.get('bam2bw', False):
            self.args_bam2bw()


    def config(self):
        """
        Config for all
        """
        ## update global
        self.args_global()

        ## update commands
        self.check_cmd()

        ## save all variable to "self" object
        return self



def args_init(x, **kwargs):
    """
    Alternative wrapper for arguments, (temp)
    replace: hiseq.utils.args.args_init()
    """
    assert isinstance(x, dict)
    args = ArgumentsInit(x, **kwargs).dict.__dict__
    args.pop('args_input', None)
    args.pop('cmd_input', None)
    args.pop('dict', None)
    return args



# ## deprecated: 2019-12-11
# ## not use dict()
# def args_init(kwargs={}, demx=False, trim=False, align=False, call_peak=False, 
#     bam2bw=False):
#         """
#         Inititate the arguments, assign the default values to arg
#         positional arg: smp, genome
#         for 1 read fastq:
#         """
#         args = kwargs.copy() # do not change original dict
#         if not isinstance(args, dict):
#             raise Exception('unknown argument: args=%s' % args)

#         # required 
#         args['fq1'] = args.get('fq1', None)
#         args['fq2'] = args.get('fq2', None)
#         args['outdir'] = args.get('outdir', str(pathlib.Path.cwd()))
#         assert isinstance(args['fq1'], str)

#         ## common args
#         if args['outdir'] is None:
#             args['outdir'] = str(pathlib.Path.cwd())

#         # absolute path
#         args['fq1'] = os.path.abspath(args['fq1'])
#         if isinstance(args['fq2'], str):
#             args['fq2'] = os.path.abspath(args['fq2'])
#         args['outdir'] = os.path.abspath(args['outdir'])

#         ## optional
#         args['genome_path'] = args.get('genome_path', 
#                 os.path.join(str(pathlib.Path.home()), 'data', 'genome'))            
#         args['overwrite'] = args.get('overwrite', False)
#         args['threads'] = args.get('threads', 8)
#         args['smp_name'] = args.get('smp_name', None)

#         ## demx
#         if demx:
#             args['demx_type'] = args.get('demx_type', 'p7') # p7, barcode, both
#             args['n_mismatch'] = args.get('n_mismatch', 0)
#             args['bc_n_left'] = args.get('bc_n_left', 3)
#             args['bc_n_right'] = args.get('bc_n_right', 2)
#             args['bc_in_read'] = args.get('bc_in_read', 1)
#             args['cut'] = args.get('cut', False)

#         ## trimming, cutadapt
#         if trim:
#             # output filename
#             fqname = file_prefix(args['fq1'])[0]
#             if not args['fq2'] is None:
#                 fqname = re.sub('_[rR]?1$', '', fqname)
#             args['fqname'] = fqname
#             args['fq_out_prefix'] = os.path.join(args['outdir'], fqname)

#             args['len_min'] = args.get('len_min', 15)
#             args['adapter3'] = Adapter().adapters[0] # TruSeq
#             args['adapter5'] = args.get('adapter5', None)
#             args['keep_name'] = args.get('keep_name', True)
#             args['gzip'] = args.get('gzip', True) # output fastq
#             args['save_temp'] = args.get('save_temp', False) # save temp files

#             args['qual_min'] = args.get('qual_min', 20)
#             args['error_rate'] = args.get('error_rate', 0.1)
#             args['overlap'] = args.get('overlap', 3)
#             args['percent'] = args.get('percent', 80)
#             args['trim_times'] = args.get('trim_times', 1)

#             args['rm_untrim'] = args.get('rm_untrim', False)
#             args['save_untrim'] = args.get('save_untrim', False)
#             args['save_too_short'] = args.get('save_too_short', False)
#             args['save_too_long'] = args.get('save_too_long', False)

#             args['adapter_sliding'] = args.get('adapter_sliding', False)
#             args['cut_before_trim'] = args.get('cut_before_trim', '0') # NSR
#             args['cut_after_trim'] = args.get('cut_after_trim', '0') # NSR
#             args['rmdup'] = args.get('rmdup', False)
#             args['cut_after_rmdup'] = args.get('cut_after_rmdup', '0')
#             args['cut_to_length'] = args.get('cut_to_length', 0)
            

#             ## PE trimming options
#             args['AD3'] = args.get('AD3', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
#             args['AD5'] = args.get('AD5', None)
#             args['not_trim_adapter'] = args.get('not_trim_adapter', False)
#             args['keep_temp_files'] = args.get('keep_temp_files', False)

#             ## deprecated
#             # args['rm_dup'] = args.get('rm_dup', True) # Deprecated since v0.3
#             # args['cut_before_trim'] = args.get('cut_before_trim', '0') # Deprecated since v0.3
#             # args['gzipped'] = args.get('gzipped', True) # Deprecated since v0.3
#             # args['double_trim'] = args.get('double_trim', False) # Deprecated since v0.3

#         ## alignment
#         if align:
#             args['genome'] = args.get('genome', None)
#             args['spikein'] = args.get('spikein', None)
#             #args['index_ext'] = args.get('index_ext', None)
#             args['extra_index'] = args.get('extra_index', None)
#             args['unique_only'] = args.get('unique_only', True) # unique map
#             args['aligner'] = args.get('aligner', 'bowtie') # bowtie alignment
#             args['te_index'] = args.get('te_index', None) #
#             # args['align_to_te'] = args.get('align_to_te', False) # deprecated
#             args['align_by_order'] = args.get('align_by_order', True) # align reads to multiple index by order
#             args['n_map'] = args.get('n_map', 0)

#             args['align_to_chrM'] = args.get('align_to_chrM', False)
#             args['align_to_rRNA'] = args.get('align_to_rRNA', False)
#             args['align_to_MT_trRNA'] = args.get('align_to_MT_trRNA', False)

#             args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
#             args['merge_rep'] = args.get('merge_rep', True)
#             args['small_genome'] = args.get('small_genome', False)
#             args['simplify_name'] = args.get('simplify_name', True)
#             args['simple_name'] = args.get('simple_name', False) # deprecated, instead simplify_name
#             args['index_parallel'] = args.get('index_parallel', False) # for multiple index
#             args['extra_para'] = args.get('extra_para', None)

#             # check-point
#             if args['spikein'] == args['genome']:
#                 args['spikein'] = None

#         ## peak-calling
#         if call_peak:
#             args['peak_caller'] = args.get('peak_caller', 'pyicoclip')

#             ## rtstop-calling
#             args['threshold'] = args.get('threshold', 1)
#             args['intersect'] = args.get('intersect', 0)
#             args['threads'] = args.get('threads', 8)

#         ## bam2bw
#         if bam2bw:
#             args['filterRNAstrand'] = args.get('filterRNAstrand', None)
#             args['samFlagExclude'] = args.get('samFlagExclude', None)
#             args['samFlagInclude'] = args.get('samFlagInclude', None)
#             args['binsize'] = args.get('binsize', 10)

#         return args


