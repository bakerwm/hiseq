
# -*- coding: utf-8 -*-

"""
Default arguments for subcommands
"""

import os
import pathlib


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


def args_init(args={}, demx=False, trim=False, align=False, call_peak=False, 
    bam2bw=False):
        """
        Inititate the arguments, assign the default values to arg
        positional arg: smp, genome
        """
        if not isinstance(args, dict):
            raise Exception('unknown argument: args=%s' % args)

        ## common args
        args['fq'] = args['fq1'] = args.get('fq', None)
        # args['fq'] = args['fq'] # for read1
        args['fq2'] = args.get('fq2', None)
        args['outdir'] = args.get('outdir', str(pathlib.Path.cwd()))
        if args['outdir'] is None:
            args['outdir'] = str(pathlib.Path.cwd())

        ## optional
        args['genome_path'] = args.get('genome_path', 
                os.path.join(str(pathlib.Path.home()), 'data', 'genome'))            
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)
        args['smp_name'] = args.get('smp_name', None)

        ## demx
        if demx:
            args['demx_type'] = args.get('demx_type', 'p7') # p7, barcode, both
            args['n_mismatch'] = args.get('n_mismatch', 0)
            args['bc_n_left'] = args.get('bc_n_left', 3)
            args['bc_n_right'] = args.get('bc_n_right', 2)
            args['bc_in_read'] = args.get('bc_in_read', 1)
            args['cut'] = args.get('cut', False)

        ## trimming, cutadapt
        if trim:
            args['len_min']  = args.get('len_min', 15)
            args['adapter3'] = Adapter().adapters[0] # TruSeq 
            args['keep_name'] = args.get('keep_name', True)
            args['gzip'] = args.get('gzip', True) # output fastq
            args['save_temp'] = args.get('save_temp', False) # save temp files

            args['qual_min'] = args.get('qual_min', 20)
            args['error_rate'] = args.get('error_rate', 0.1)
            args['overlap'] = args.get('overlap', 3)
            args['percent'] = args.get('percent', 80)
            args['trim_times'] = args.get('trim_times', 1)

            args['rm_untrim'] = args.get('rm_untrim', False)
            args['save_untrim'] = args.get('save_untrim', False)
            args['save_too_short'] = args.get('save_too_short', False)
            args['save_too_long'] = args.get('save_too_long', False)

            args['adapter_sliding'] = args.get('adapter_sliding', False)
            args['cut_before_trim'] = args.get('cut_before_trim', '0') # NSR
            args['cut_after_trim'] = args.get('cut_after_trim', '0') # NSR
            args['rmdup'] = args.get('rmdup', False)
            args['cut_after_rmdup'] = args.get('cut_after_rmdup', '0')
            args['cut_to_length'] = args.get('cut_to_length', 0)
            args['adapter5'] = args.get('adapter5', None)

            ## PE trimming options
            args['AD3'] = args.get('AD3', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
            args['AD5'] = args.get('AD5', None)
            args['not_trim_adapter'] = args.get('not_trim_adapter', False)
            args['keep_temp_files'] = args.get('keep_temp_files', False)

            ## deprecated
            # args['rm_dup'] = args.get('rm_dup', True) # Deprecated since v0.3
            # args['cut_before_trim'] = args.get('cut_before_trim', '0') # Deprecated since v0.3
            # args['gzipped'] = args.get('gzipped', True) # Deprecated since v0.3
            # args['double_trim'] = args.get('double_trim', False) # Deprecated since v0.3


        ## alignment
        if align:
            args['genome'] = args.get('genome', None)
            args['spikein'] = args.get('spikein', None)
            #args['index_ext'] = args.get('index_ext', None)
            args['extra_index'] = args.get('extra_index', None)
            args['unique_only'] = args.get('unique_only', True) # unique map
            args['aligner'] = args.get('aligner', 'bowtie') # bowtie alignment
            args['te_index'] = args.get('te_index', None) #
            # args['align_to_te'] = args.get('align_to_te', False) # deprecated
            args['align_by_order'] = args.get('align_by_order', True) # align reads to multiple index by order
            args['n_map'] = args.get('n_map', 0)
            args['align_to_rRNA'] = args.get('align_to_rRNA', True)
            args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
            args['merge_rep'] = args.get('merge_rep', True)
            args['small_genome'] = args.get('small_genome', False)
            args['simplify_name'] = args.get('simplify_name', True)
            args['simple_name'] = args.get('simple_name', False) # deprecated, instead simplify_name

            # check-point
            if args['spikein'] == args['genome']:
                args['spikein'] = None

        ## peak-calling
        if call_peak:
            args['peak_caller'] = args.get('peak_caller', 'pyicoclip')

            ## rtstop-calling
            args['threshold'] = args.get('threshold', 1)
            args['intersect'] = args.get('intersect', 0)
            args['threads'] = args.get('threads', 8)

        ## bam2bw
        if bam2bw:
            args['filterRNAstrand'] = args.get('filterRNAstrand', None)
            args['samFlagExclude'] = args.get('samFlagExclude', None)
            args['samFlagInclude'] = args.get('samFlagInclude', None)
            args['binsize'] = args.get('binsize', 10)

        return args


