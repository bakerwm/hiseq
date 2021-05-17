#!/usr/bin/env python3 

"""For aligner config
"""

import pyfastx
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import *
from align_index import *
from utils import *



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
        
        
        