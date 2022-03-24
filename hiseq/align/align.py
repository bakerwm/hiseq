#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
!!! switch from TOML to YAML, because of None type in python not saved in TOML format !!!

The structure of the classes:

Align()
AlignRx() : fq=n, index=n, multiple fastq to multiple index
AlignRn() : fq=1, index=n, single fastq to multiple index
AlignR1() : fq=1, index=1, single fastq to single index
AlignRp() : report()

AlignConfig()
AlignRxConfig()
AlignRnConfig()
AlignR1Config()
AlignRpConfig()

Working model:
Align {AlignConfig()}
  |- AlignRx()
  |- AlignRn()
  |- AlignR1()


AlignR1()
  |- {bowtie2, bowtie, STAR, bwa, ...}
"""

# import sys
# import re
# import shutil
# import logging
# import tempfile
# import collections
# import pandas as pd
# from Levenshtein import distance

import os
from pathlib import Path # cwd
import argparse
from multiprocessing import Pool
from hiseq.utils.seq import Fastx
from hiseq.align.align_rp import AlignRp
from hiseq.utils.utils import (
    update_obj, Config, get_date, read_hiseq, is_supported
)
from hiseq.utils.file import (
    check_path, file_abspath, file_prefix, file_exists, symlink_file,
    list_dir, check_fx_args, check_fx_paired, fx_name
)
from hiseq.align.bowtie import Bowtie
from hiseq.align.bowtie2 import Bowtie2
from hiseq.align.hisat2 import Hisat2
from hiseq.align.star import Star
from hiseq.align.salmon import Salmon
# from hiseq.align.utils import check_fx_args
from hiseq.align.align_index import AlignIndex, check_index_args


class Align(object):
    """The main port: fx: N; index: N
    support: multiple fx, multiple index

    force to list, even if `None`
    fq2: [None, None, ...]
    smp_name: [str, str, ...]
    index_name: [str, str, ...]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = AlignConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        Config().dump(self.__dict__, self.config_yaml)


    def run(self):
        msg = '\n'.join([
            '-'*80,
            '{:>14s} : {}'.format('Date', get_date()),
            '{:>14s} : {}'.format('program', 'Align'),
            '{:>14s} : {}'.format('config', self.config_yaml),
            '{:>14s} : {}'.format('aligner', self.aligner),
            '{:>14s} : {}'.format('genome', self.genome),
            '{:>14s} : {}'.format('index_list', self.index_list),
            '{:>14s} : {}'.format('extra_index', self.extra_index),
            '{:>14s} : {}'.format('fq1', self.fq1),
            '{:>14s} : {}'.format('fq2', self.fq2),
            '{:>14s} : {}'.format('outdir', self.outdir),
            '{:>14s} : {}'.format('unique_only', self.unique_only),
            '{:>14s} : {}'.format('extra_para', self.extra_para),
            '{:>14s} : {}'.format('threads', self.threads),
            '{:>14s} : {}'.format('parallel_jobs', self.parallel_jobs),
            '{:>14s} : {}'.format('to_rRNA', self.to_rRNA),
            '{:>14s} : {}'.format('to_chrM', self.to_chrM),
            '{:>14s} : {}'.format('to_MT_trRNA', self.to_MT_trRNA),
            '-'*80,
        ])
        print(msg)
        AlignRx(**self.__dict__).run()


class AlignRx(object):
    """
    Alignment: fx: N; index: N (same as: Align)
    Align to each index, one-by-one; 
    Align fq files in parallel, if specified;
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        
        
    def run_single_fx(self, i):
        """Run Alignment for single fx file
        parallel
        """
        args_local = self.__dict__.copy()
        # update fq1, fq2, smp_name
        args_local.update({
            'fq1': self.fq1[i],
            'fq2': self.fq2[i] if isinstance(self.fq2, list) else None,
            'smp_name': self.smp_name[i],
            'rep_list': None,
            'is_paired': None,
            })
        AlignRn(**args_local).run()


    def run(self):
        # run in parallel
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fx, range(len(self.fq1)))
        else:
            [self.run_single_fx(i) for i in range(len(self.fq1))]
        # add project_dir
        self.project_dir = self.outdir
        AlignRp(**self.__dict__).run()


class AlignRn(object):
    """Alignment: fx: 1, index: N
    Align to each index, one-by-one, not need: parallel
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        global config
        """
        args_local = AlignRnConfig(**self.__dict__) # update
        self = update_obj(self, args_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)


    def wrap_stat(self):
        """
        Save alignment to one file
        all files in project_dir
        """
        # save the last index-bam as output
        last_bam = None
        last_stat = None
        last_flagstat = None
        df = {} 
        r1_dirs = list_dir(self.project_dir, include_dir=True)
        for r1 in r1_dirs:
            a = read_hiseq(r1, 'r1')
            if a.is_hiseq:
                s1 = Config().load(getattr(a, 'align_json'))
                df.update(s1)
                last_bam = a.bam
                last_stat = a.align_stat
                last_json = a.align_json
                last_flagstat = a.align_flagstat
                last_unmap = a.unmap
                last_unmap1 = a.unmap1
                last_unmap2 = a.unmap2
        symlink_file(last_bam, self.bam) # save bam
        symlink_file(last_stat, self.align_stat) # save flagstat
        symlink_file(last_json, self.align_json) # save flagstat
        symlink_file(last_flagstat, self.align_flagstat) # save flagstat
        symlink_file(last_unmap, self.unmap)
        symlink_file(last_unmap1, self.unmap1)
        symlink_file(last_unmap2, self.unmap2)
#         Config().dump(df, self.align_json) # save stat


    def run(self):
        args_local = self.__dict__.copy()
        for index,index_name in zip(self.index_list, self.index_name):
            args_local.update({
                'index': index,
                'index_name': index_name,
                'index_list': None,
                'rep_list': None,
                'is_paired': None,
            })
            bam, unmap1, unmap2 = AlignR1(**args_local).run()
            # update the unmap files for next round
            args_local.update({
                'fq1': unmap1,
                'fq2': unmap2,
            })
        # combine stat
        self.wrap_stat()
        AlignRp(**self.__dict__).run()


class AlignR1(object):
    """Alignment: fx: 1; index: 1
    return: bam, unmap-1, unmap-2
    save: log, stat, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = AlignR1Config(**self.__dict__) # update
        self = update_obj(self, args_local.__dict__, force=True)


    def run(self):
        """Run specific aligner
        bowtie
        bowtie2
        STAR
        ...
        """
        aligner = {
            'bowtie': Bowtie,
            'bowtie2': Bowtie2,
            'star': Star,
            'salmon': Salmon,
            'hisat2': Hisat2,
#             'bwa': BWA,
#             'kallisto': Kallisto,
        }
        port = aligner.get(self.aligner.lower(), None)
        if port is None:
            raise ValueError('unknown aligner: {}'.format(self.aligner))
        args_local = self.__dict__.copy()
        align = port(**args_local)
        return align.run() # bam, unmap-1, unmap-2


class AlignConfig(object):
    """The main port: fx: N; index: N
    force, even if `None`
    fq1: [fq1]
    fq2: None or list
    smp_name: [str, str, ...]
    index_name: [str, str, ...]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        if not is_supported(self.aligner, key='supported_aligner'):
            raise ValueError('-a {}, unknown'.format(self.aligner))
        if isinstance(self.genome, str) and not is_supported(self.genome, key='supported_genome'):
            raise ValueError('-g {}, unknown'.format(self.genome))


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        see: check_index_args for default arguments
        """
        args_init = {
            'aligner': 'bowtie',
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'smp_name': None,
            'unique_only': False,
            'index_list': None,
            'extra_index': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'to_rRNA': False,
            'rRNA_index': None,
            'to_chrM': False,
            'to_MT_trRNA': False,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'keep_tmp': False,
            'genome_size': 0,
            'genome_size_file': None,
            'extra_para': None,
            'verbose': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'alignment_rx' # force ?@!
        if self.outdir is None:
            self.outdir = str(Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.init_fq()
        self.init_index()
        self.init_files()


    def init_fq(self):
        """Force fx to list
        fq1 required. to_list
        fq2 optional, None
        smp_name to_list
        """
        # force fq1 to list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        if not check_fx_args(self.fq1, self.fq2, check_empty=True):
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        snames = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)
        if isinstance(self.smp_name, str):
            self.smp_name = [self.smp_name]
        if isinstance(self.smp_name, list):
            if len(self.smp_name) == len(self.fq1):
                if all([isinstance(i, str) for i in self.smp_name]):
                    snames = self.smp_name
        self.smp_name = snames
        # update rep_list
        self.rep_list = [os.path.join(self.outdir, i) for i in self.smp_name]
    

    def init_index(self):
        self.index_list = check_index_args(**self.__dict__)
        if len(self.index_list) == 0:
            raise ValueError('no index found')
        self.index_name = [AlignIndex(i).index_name() for i in self.index_list]
        # auto index_name, 01, 02, ...
        if len(self.index_name) > 1:
            self.index_name = ['{:02d}_{}'.format(i, n) for i,n in zip(
                range(len(self.index_name)), self.index_name)]


    def init_files(self):
        """
        Config only
        config_yaml
        """
        self.config_dir = os.path.join(self.outdir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        check_path(self.config_dir, create_dirs=True)


class AlignRnConfig(object):
    """The main port: fx: 1; index: N
    index_list
    index_name

    force, even if `None`
    *fq1: str
    *fq2: None or str
    smp_name: str
    index_name: [str, str, ...]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'fq1': None,
            'fq2': None,
            'smp_name': None,
            'index_list': None,
            'index_name': None,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'alignment_rn'
        self.init_fq()
        self.init_index()
        self.init_files()


    def init_fq(self):
        """
        fq1 : str required. 
        fq2 : None or str, optional, None
        smp_name ： str
        index_list : list
        index_name : list
        """
        if isinstance(self.fq1, str):
            flag = check_fx_args(self.fq1, self.fq2)
        else:
            flag = False
        if not flag:
            raise ValueError('fq1, fq2 not valid, str expect, got {}'.format(
            type(self.fq1).__name__))
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.fx_format = Fastx(self.fq1).format # fasta/q
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)
        # update rep_list
        self.rep_list = os.path.join(self.outdir, self.smp_name)

            
    def init_index(self):
        """Force: 
        index_list: list
        index_name: list
        """
        flag_err = True
        # update index_name
        if self.index_name is None:
            self.index_name = [AlignIndex(i, self.aligner).index_name() \
                               for i in self.index_list]
        if isinstance(self.index_list, list) and \
            isinstance(self.index_name, list):
            flag_err = not len(self.index_list) == len(self.index_name)
        if flag_err:
            raise ValueError('index_list, index_name, not valid')


    def init_files(self):
        """Config only
        config_yaml
        """
        self.project_dir = os.path.join(self.outdir, self.smp_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        prefix = os.path.join(self.project_dir, self.smp_name)
        default_args = {
            'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
            'bam': prefix + '.bam',
            'unmap': prefix + '.unmap.' + self.fx_format,
            'unmap1': prefix + '.unmap.1.' + self.fx_format, # 
            'unmap2': prefix + '.unmap.2.' + self.fx_format, # 
            'align_stat': prefix + '.align.stat',
            'align_json': prefix + '.align.json',
            'align_flagstat': prefix + '.align.flagstat',
        }
        self = update_obj(self, default_args, force=True)
        check_path(self.config_dir, create_dirs=True)


class AlignR1Config(object):
    """The main port: fx: 1; index: 1
    check arguments by Aligner()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'fq1': None,
            'fq2': None,
            'smp_name': None,
            'aligner': None,
            'index': None,
            'index_name': None,
            'unique_only': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'alignment_r1'
        self.init_fq()
        self.init_index()


    def init_fq(self):
        """Force fx to list
        fq1 : str required. 
        fq2 : None or str, optional, None
        smp_name ： str
        index_list : list
        index_name : list
        """
        if isinstance(self.fq1, str):
            flag = check_fx_args(self.fq1, self.fq2)
        else:
            flag = False
        if not flag:
            raise ValueError('fq1, fq2 not valid, str expect, got {}'.format(
            type(self.fq1).__name__))
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.fx_format = Fastx(self.fq1).format # fasta/q
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)
        # update rep_list
        self.rep_list = os.path.join(self.outdir, self.smp_name)


    def init_index(self):
        """Force:
        index: str
        index_name: str
        :for single fq:
        """
        flag_err = True
        # update index_name
        if isinstance(self.index, str):
            flag_err = not AlignIndex(self.index, self.aligner).is_valid()
            if self.index_name is None:
                self.index_name = AlignIndex(self.index).index_name()
        else:
            raise ValueError('index, not valid: {}'.format(self.index))


def get_args():
    """Parsing arguments for Align
    The main port, aligner
    """
    example = '\n'.join([
        'Examples:',
        '$ python align.py -1 f1.fq -x genome -o output',
        '# add extra para',
        '$ python align.py -1 f1.fq -2 f2.fq -x genome -o output -X "-X 2000"',
        '# unique reads, update index_name',
        '$ python align.py -1 f1.fq -x genome -o output -u -in 01.genome',
    ])    
    parser = argparse.ArgumentParser(
        prog='align',
        description='run algner {bowtie|bowtie2|STAR}',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--aligner', default='bowtie2', type=str,
        help='The aligner for alignment, default: [bowtie2]')
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='Fasta/q file, read1 of PE, or SE read')
    parser.add_argument('-2', '--fq2', nargs='+', required=False, default=None,
        help='Fasta/q file, read2 of PE, or SE read, optional')
    parser.add_argument('-o', '--outdir', default=None,
        help='Directory saving results, default: [cwd]')
    parser.add_argument('-n', '--smp-name', nargs='+', dest='smp_name', default=None, 
        help='The name of the sample')

    # index
    parser.add_argument('--index-list', dest="index_list", nargs='+',
        help='A list of index')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('-g', '--genome', default=None, 
        help='The name of the genome, [dm6, hg38, mm10]')
    parser.add_argument('--genome-index', dest="genome_index", default=None,
        help='align index of genome')
    parser.add_argument('-k', '--spikein', default=None,
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--spikein-index', dest="spikein_index", default=None,
        help='align index of spikein')    
    parser.add_argument('--to-rRNA', dest='to_rRNA', action='store_true',
        help='Align to rRNA')
    parser.add_argument('--rRNA-index', dest="rRNA_index", default=None,
        help='align index of rRNA')   
    parser.add_argument('--to-chrM', action='store_true', dest='to_chrM',
        help='Align reads to mitochromosome first')
    parser.add_argument('--to-MT-trRNA', action='store_true', dest='to_MT_trRNA',
        help='Align reads to chrM, tRNA and rRNAs first')    
    
    # extra para
    parser.add_argument('-p', '--threads', default=4, type=int,
        help='Number of threads, default: [4]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        help='Number of jobs to run in parallel, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite the exist files')
    parser.add_argument('-u', '--unique-only', action='store_true',
        dest='unique_only', help='Report unique mapped reads only')
    parser.add_argument('--n-map', dest='n_map', type=int, default=0,
        help='Number of hits per read, default [0]')
    parser.add_argument('-l', '--largs-insert', action='store_true',
        dest='large_insert',
        help='For large insert, use: -X 1000 --chunkmbs 128')
    parser.add_argument('--keep-tmp', dest='keep_tmp', action='store_true',
                        help='save temp files')
    parser.add_argument('--rm-unmap', dest='keep_unmap', action='store_false',
                        help='remove unmap fastq files')
    parser.add_argument('-X', '--extra-para', dest='extra_para', default=None,
        help='Add extra parameters, eg: "-X 2000"')
    parser.add_argument('--verbose', action='store_true', 
        help='Show message in details')
    return parser


def main():
    args = vars(get_args().parse_args())
    # update: keep_tmp, keep_unmap
    # args['keep_unmap'] = args['keep_tmp']
    Align(**args).run()


if __name__ == '__main__':
    main()
    
    
    
"""The port for alignment
support: multi fx files, multi indexes

Align()
AlignRx() : fq=n, index=n, multiple fastq to multiple index
AlignRn() : fq=1, index=n, single fastq to multiple index
AlignR1() : fq=1, index=1, single fastq to single index
AlignRp() : report()

AlignConfig()
AlignRxConfig()
AlignRnConfig()
AlignR1Config()
AlignRpConfig()

Working model:
Align {AlignConfig()}
  |- AlignRx()
  |- AlignRn()
  |- AlignR1()


AlignR1()
  |- {bowtie2, bowtie, STAR, bwa, ...}

## requirements
- unique mapper
- multiple mapper
- config (yaml)
- saving unmapped ?!

## index
genome + {spikein, rRNA, MT}
extra_index {list}

top-level for alignment

input:
  - pickle: str or None
  - fq1: list (required)
  - fq2: list or None
  - aligner: str (required)
  - genome: str or None
  - spikein: str or None
  - align-to-rRNA: bool
  - align-to-chrM: bool
  - smp_name: list or None (auto)
  - index_list: list or None
  - index_name: list or None (auto)

  - index_list_equal: bool

  - n_map
  - unique_only
  - extra_para
  - parallel_jobs
  - threads

output:
  - bam: str
  - unmap1: str
  - unmap2: str or None 

priority:
  - 1. pickle
  - 2. genome, spikein, align-to-chrM, ...
  - x. extra_index

return bam, unmap1, unmap2


changelog

## update: 2021-03-08
1. run the codes independently
2. force 'bowtie2' to align rRNAs/chrM/...

## update: 2020-01-08
1. uniform code style: self ->  dict -> arguments.txt + pickle

## update: 2020-04-24
1. rewrite the script, frame updated

to-do
1. force rRNA/chrM, unique + multiple

"""