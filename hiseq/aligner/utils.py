#!/usr/bin/env python3

"""Functions for aligner
check fastq/index, arguments, ...

file existence 
arguments 

"""

import pyfastx
from hiseq.utils.helper import *












def check_index_args(**kwargs):
    """Check the index for aligner
    exists
    valid
    aligner
    name
    ...

    Parameters
    ----------
    genome
    genome_index
    spikein
    spikein_index
    index_list (ignore all)
    index_extra (append)
    to_rRNA
    to_MT
    to_chrM
    to_MT_trRNA
    """
    # default arguments
    args = {
        'genome': None,
        'genome_index': None,
        'spikein': None,
        'spikein_index': None,
        'index_list': None,
        'index_extra': None,
        'to_rRNA': False,
        'to_MT': False,
        'to_chrM': False,
        'to_MT_trRNA': False
    }
    args.update(kwargs)
    # level-1
    if isinstance(args['index_list'], list):
        pass
    else:
        # spikein
        if isinstance(args['spikein_index'], str)


























def check_fx_args(fq1, fq2=None, **kwargs):
    """Check the fastx in arguments
    
    Parameters
    ----------
    fq1 : str or list
    fq2 : None, str or list

    check:
    1. file exists
    2. file type
    3. fq paired
    4. check_empty
    """
    if isinstance(fq1, str):
        out = check_fx(fq1, **kwargs)
        if isinstance(fq2, str):
            k2 = check_fx(fq2, **kwargs)
            p1 = check_fx_paired(fq1, fq2, **kwargs)
            out = all([out, k2, p1])
    elif isinstance(fq1, list):
        out = all([check_fx(i, **kwargs) for i in fq1])
        if isinstance(fq2, list):
            k2 = all([check_fx(i, **kwargs) for i in fq2])
            out = out and k2
            if len(fq1) == len(fq2):
                p2 = all([check_fx_paired(a, b) for a,b in zip(fq1, fq2)])
                out = out and p2
    else:
        log.error('fq1 expect str or list, got {}'.format(type(fq1).__name__))
        out = False
    return out


# fq files
def check_fx(fx, **kwargs):
    """Check the fastq/a files
    1. file exist
    2. fq1 required
    3. fq type: fasta/q
    """
    kwargs['check_empty'] = kwargs.get('check_empty', True)
    kwargs['show_error'] = kwargs.get('show_error', False)
    if isinstance(fx, str):
        if check_file(fx, **kwargs):
            try:
                fx_type = Fastx(fx).format # fasta/q
                out = fx_type in ['fasta', 'fastq']
            except ValueError as err:
                log.info('Failed to read file, with error: {}'.format(err))
                out = False
            finally:
                log.info('check_fx() ...')            
        else:
            if kwargs['show_error']:
                log.error('fx failed, {}'.format(fx))
            out = False
    else:
        out = False
    return out


def check_fx_paired(fq1, fq2, **kwargs):
    """Check the fq1 and fq2 are paired or not
    the name of fq1, fq2 are the same

    fq1: @the-name/1
    fq2: @the-name/2

    scan the first read from the two files
    """
    if isinstance(fq1, str) and isinstance(fq2, str):
        if check_fx(fq1, fq2):
            # check fq name:
            fx1 = pyfastx.Fastx(fq1)
            fx2 = pyfastx.Fastx(fq2)
            for a,b in zip(fx1, fx2):
                out = a[0][:-1] == b[0][:-1]
                break
        else:
            out = False
    else:
        out = False
    return out


def check_file(x, **kwargs):
    """Check the x file
    1. file exists
    """
    show_error = kwargs.get('show_error', False)
    check_empty = kwargs.get('check_empty', False)
    if isinstance(x, str):
        if file_exists(x):
            x_size = os.stat(x).st_size
            # empty gzipped file, size=20
            out = x_size > 20 if check_empty else True
        else:
            if show_error:
                log.error('file not exists: {}'.format(x))
            out = False # failed
    elif isinstance(x, list):
        out = all([check_file(i, **kwargs) for i in x])
    else:
        log.error('x expect str or list, got {}'.format(type(x).__name__))
        out = False
    return out





class AlignerConfig(object):
    """
    Config for single alignment, 1 fq, 1 index
  
    [1, 1] : 1 fastq, 1 index
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        default arguments
        required:
          - fq1
          - index
          - outdir
        
        optional:
          - aligner
        """
        args_default = {
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index': None,            
            'index_name': None,
            'smp_name': None,
            'extra_para': None,
            'smp_name': None,
            'threads': 1,
            'overwrite': False,
            'n_map': 1,
            'unique_only': False,
            'genomeLoad': 'NoSharedMemory',
            'genome_size': 0,
            'genome_size_file': None,
            'keep_tmp': False
        }
        self = update_obj(self, args_default, force=False)

        # outdir
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # smp_name
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        ## update
        self.init_fq()
        self.init_index()
        init_cpu()
        self.init_files()
        self.is_paired = True if self.fq2 else False


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (optional)
        """
        # fq1
        if not isinstance(self.fq1, str):
            raise ValueError('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))

        # abs
        self.fq1 = file_abspath(self.fq1)

        # file exists
        if not file_exists(self.fq1):
            raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

        # fq2
        if not self.fq2 is None:
            if not isinstance(self.fq2, str):
                raise ValueError('--fq2, None or str expected, got {}'.format(
                    type(self.fq2).__name__))

            # abs
            self.fq2 = file_abspath(self.fq2)

            if not file_exists(self.fq2):
                raise ValueError('--fq2, file not exists: {}'.format(self.fq2))

            # paired
            if not fq_paired(self.fq1, self.fq2):
                raise ValueError('--fq1, --fq2, file not paired')

        ## fq type (fa or fq)
        self.fq_format = self.fq_type = Fastx(self.fq1).format # to-do !!!


    def init_index(self):
        """
        alignment index = {genome_index}
        genome, extra_index, genome_index

        output: genome_index
        """
        if isinstance(self.index, str):
            # get name
            if not isinstance(self.index_name, str):
                self.index_name = AlignIndex(index=self.index).index_name()

            ai = AlignIndex(aligner=self.aligner, index=self.index)

            if not ai.is_index():
                raise ValueError('index: {} is not for {}'.format(
                    self.index, self.aligner))

            if self.genome_size < 1:
                self.genome_size = ai.index_size()

            if not isinstance(self.genome_size_file, str): 
                self.genome_size_file = ai.index_size(return_file=True)

        else:
            raise ValueError('index failed: {}'.format(self.index))


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
            'config_toml': os.path.join(subdir, 'config.toml'),
            'cmd_shell': os.path.join(subdir, 'cmd.sh'),
            'bam': prefix + '.bam',
            'sam': prefix + '.sam',
            'unmap': prefix + '.unmap.fastq',
            'unmap1': prefix + '.unmap.1.fastq',
            'unmap2': prefix + '.unmap.2.fastq',
            'align_log': prefix + '.align.log',
            'align_stat': prefix + '.align.stat',
            'align_toml': prefix + '.align.toml',
            'align_flagstat': prefix + '.flagstat'
        }
        self = update_obj(self, default_files, force=True)