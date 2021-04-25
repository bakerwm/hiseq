#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
###############
1. unmap files:
  - bowtie  : _1.fastq, _2.fastq
  - bowtie2 : .1.fastq, .2.fastq
  - STAR    : .unmap.1.fastq (default: Unmapped.out.mate1/2)
###############

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
- config (toml)
- saving unmapped ?!

## index
genome + {spikein, rRNA, MT}
extra_index {list}



#######################################################
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
#######################################################


## update: 2020-01-08
1. uniform code style: self ->  dict -> arguments.txt + pickle

## update: 2020-04-24
1. rewrite the script, frame updated

to-do
1. force rRNA/chrM, unique + multiple

"""

import os
import sys
import re
import pathlib
import shutil
import logging
import tempfile
import collections
import pandas as pd
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions


def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


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


class Align(object):
    """
    Main port for Alignment

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()

    def init_args(self):
        args_local = AlignConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update


    def run(self):
        args_local = self.__dict__
        if self.alignment_type == 'alignment_rx':
            AlignRx(**args_local).run()
        elif self.alignment_type == 'alignment_rn':
            AlignRn(**args_local).run()
        elif self.alignment_type == 'alignment_r1':
            AlignR1(**args_local).run()
        else:
            raise ValueError('unknown alignment_type, check fq1')


class AlignRx(object):
    """
    Multiple index for multiple fastq

    {AlignRn, AlignR1}

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = AlignRxConfig(**self.__dict__) # udpate
        self = update_obj(self, args_local, force=True) # update
        Toml(self.__dict__).to_toml(self.config_toml)


    def run_single_fq(self, i):
        """Run Alignment for single fastq file
        could be run in parallel

        Pool
        """
        args_local = self.__dict__.copy()
        # update fq1, fq2, smp_name
        args_local['fq1'] = self.fq1[i]
        args_local['fq2'] = self.fq2[i] if isinstance(self.fq2, list) else None
        args_local['smp_name'] = self.smp_name[i] \
            if isinstance(self.smp_name, list) else None

        AlignRn(**args_local).run()


    def run(self):
        i_list = list(range(len(self.fq1)))
        # run in parallel
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fq, i_list)
        else:
            # map(self.run_single_fq, i_list)
            [self.run_single_fq(i) for i in i_list]


class AlignRn(object):
    """
    Multiple index for single fastq

    fq=1, index=n


    {AlignR1}

    Example:
    1. 

    2. 

    3. 
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
        Toml(self.__dict__).to_toml(self.config_toml)


    def run_single_index(self, i):
        """Run for single index
        Iterate index_list
        """
        args_local = self.__dict__.copy()
        args_local['index'] = self.index_list[i]
        args_local['index_name'] = self.index_name[i]
        args_local.pop('index_list', None) # remove index_list
        # args_local['keep_tmp'] = True # force
        ai = AlignR1(**args_local)
        # align
        ai.run()

        # update unmap
        self.fq1 = ai.unmap1
        self.fq2 = ai.unmap2 if file_exists(ai.unmap2) else None


    def wrap_stat(self):
        """Save alignment to one file

        to-do
        """
        pass


    def run(self):
        i_list = list(range(len(self.index_list)))
        if len(self.index_list) > 1:
            self.keep_tmp = True # force, for multiple index, alignment
        [self.run_single_index(i) for i in i_list]
        # for index in self.index_list:
        #     self.run_single_index(index)

        self.wrap_stat()


class AlignR1(object):
    """
    Single index for single fastq

    {Bowtie2, ...}

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        global config
        """
        args_local = AlignR1Config(**self.__dict__) # update
        self = update_obj(self, args_local.__dict__, force=True)
        Toml(self.__dict__).to_toml(self.config_toml)


    def run(self):
        """
        Pick the aligner porter
        """
        aligner = {
            'bowtie': Bowtie,
            'bowtie2': Bowtie2,
            'star': STAR,
            'bwa': BWA,
            'hisat2': Hisat2,
            'kallisto': Kallisto,
            'salmon': Salmon
        }
        port = aligner.get(self.aligner.lower(), None)

        if port is None:
            raise ValueError('unknown aligner: {}'.format(self.aligner))

        args_local = self.__dict__.copy()
        align = port(**args_local)
        align.run()


class AlignRp(object):
    """
    Report for Alignment

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)


class AlignConfig(object):
    """
    Config files for Align

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': None,
            'index_list': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'extra_index': None,
            'align_to_rRNA': False,
            'align_to_MT_trRNA': False,
            'align_to_chrM': False,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'keep_tmp': False,
            'genome_size': 0,
            'genome_size_file': None
        }
        self = update_obj(self, args_init, force=False)
        self.alignment_type = 'alignment_rx'

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # default files
        self.init_files()

        # update mission
        self.init_mission()


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (optional)
        """
        # convert str to list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]

        if isinstance(self.fq1, list):
            # fq1
            self.fq1 = file_abspath(self.fq1)

            # file exists
            if not file_exists(self.fq1):
                raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

            # fq2
            if isinstance(self.fq2, list):
                self.fq2 = file_abspath(self.fq2)

                if not file_exists(self.fq2):
                    raise ValueError('--fq2, file not exists: {}'.format(
                        self.fq2))

                # paired
                if not fq_paired(self.fq1, self.fq2):
                    raise ValueError('--fq1, --fq2, file not paired')
        else:
            raise ValueError('fq1, fq2 required')


    def init_index(self):
        """
        alignment index
        genome, extra_index, genome_index

        output: index_list
        """
        print('!CCCC-1', self.genome_index)
        index_all = {
            'spikein': None,
            'tags': None, # rRNA, tRNA, chrM, ...
            'genome': None,
            'extra_index': None
        }

        # index_list
        if isinstance(self.index_list, str):
            self.index_list = [self.index_list]
        elif isinstance(self.index_list, list):
            pass
        else:
            self.spikein_index = None # init
            self.tag_index = None # init
            # self.genome_index = None # init
            ####################
            # index_list       #
            ####################
            if isinstance(self.extra_index, str):
                self.extra_index = [self.extra_index]
            elif isinstance(self.extra_index, list):
                self.extra_index = self.extra_index
            else:
                # genome, spikein, rRNA, MT, ...
                ####################
                # spikein index    #
                ####################
                if isinstance(self.spikein_index, str):
                    ai = AlignIndex(
                        index=self.spikein_index, 
                        aligner=self.aligner)
                    if not ai.is_index():
                        raise ValueError('spikein_index failed, {}'.format(
                            self.spikein_index))
                elif isinstance(self.spikein, str):
                    self.spikein_index = AlignIndex(
                        aligner=self.aligner).search(
                        genome=self.spikein, group='genome')
                else:
                    pass

                ####################
                # tag index        #
                ####################
                if self.align_to_MT_trRNA is True:
                    tag = 'MT_trRNA'            
                elif self.align_to_rRNA is True:
                    tag = 'rRNA'
                elif self.align_to_chrM is True:
                    tag = 'chrM'
                else:
                    tag = None

                if tag and isinstance(self.genome, str):
                    tag_index = AlignIndex(aligner=self.aligner).search(
                        genome=self.genome, group=tag)
                else:
                    tag_index = None

                if AlignIndex(aligner=self.aligner, index=tag_index).is_index():
                    self.tag_index = tag_index
                else:
                    self.tag_index = None

                ####################
                # genome index     #
                ####################
                if isinstance(self.genome_index, str):
                    ai = AlignIndex(index=self.genome_index)
                    # valid
                    if not ai.is_index():
                        raise ValueError('genome_index failed, {}'.format(
                            self.genome_index))
                    # chrsizes
                    if self.genome_size < 1:
                        self.genome_size = ai.index_size()
                    if not isinstance(self.genome_size_file, str): 
                        self.genome_size_file = ai.index_size(return_file=True)

                elif isinstance(self.genome, str):
                    self.genome_index = AlignIndex(
                        aligner=self.aligner).search(
                        genome=self.genome, group='genome')
                    self.genome_size_file = Genome(
                        genome=self.genome).get_fasize()
                    if self.genome_size < 1:
                        with open(self.genome_size_file, 'rt') as r:
                            s = [i.strip().split('\t')[1] \
                                for i in r.readlines()]
                        self.genome_size = sum(map(int, s))
                
                else:
                    pass
                    # self.genome_index = None

            # construct all index
            self.index_list = [] # int

            # spikein
            if self.spikein_index:
                self.index_list.append(self.spikein_index)
            if self.tag_index:
                self.index_list.append(self.tag_index)
            if self.genome_index:
                self.index_list.append(self.genome_index)
            if self.extra_index:
                self.index_list.extend(self.extra_index)

        if len(self.index_list) < 1:
            raise ValueError('index not found, check, genome, extra_index')


    def init_files(self):
        """
        default files for output
        bam, sam, log, unmap, ...
        """
        self.config_toml = self.outdir + '/config.toml'


    def init_mission(self):
        """Determine the type of Alignment
        Rx: fq=n, index=n
        Rn: fq=1, index=n
        R1: fq=1, index=1
        Rp: generate report()
        """
        if len(self.fq1) > 1: # Rx, Rn
            self.alignment_type = 'alignment_rx'
        elif len(self.fq1) == 1: # Rn, R1
            self.fq1 = self.fq1.pop()
            self.fq2 = self.fq2.pop() if isinstance(self.fq2, list) else None

            if isinstance(self.index_list, str):
                self.index = self.index_list
                self.alignment_type = 'alignment_r1'
            else:
                self.alignment_type = 'alignment_rn'
        else:
            raise ValueError('unknown Alignment type, check fq1')


class AlignRxConfig(object):
    """
    Config files for AlignRx

    fq1=n, fq2=n, index_list={list}

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': None,
            'index_list': None,
            'threads': 1,
            'overwrite': False,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)
        self.alignment_type = 'alignment_rx'

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # default files
        self.init_files()


    def init_fq(self):
        """
        Support fastq files
        fq1 (required)
        fq2 (optional)
        """
        # convert str to list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]

        if isinstance(self.fq1, list):
            # fq1
            self.fq1 = file_abspath(self.fq1)

            # file exists
            if not file_exists(self.fq1):
                raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

            # fq2
            if isinstance(self.fq2, list):
                self.fq2 = file_abspath(self.fq2)

                if not file_exists(self.fq2):
                    raise ValueError('--fq2, file not exists: {}'.format(
                        self.fq2))

                # paired
                if not fq_paired(self.fq1, self.fq2):
                    raise ValueError('--fq1, --fq2, file not paired')
        else:
            raise ValueError('fq1, fq2 required')


    def init_index(self):
        """
        alignment index_list = {genome_index}
        
        The index_list are all valid, for the aligner
        """
        if isinstance(self.index_list, list):
            chk = [AlignIndex(aligner=self.aligner).is_index(i) \
                for i in self.index_list]
            if not all(chk):
                msg = '\n'.join([
                    '{:>6s} : {}'.format(str(j), k) for j, k in zip(chk, 
                        self.index_list)])
                print(msg)
                raise ValueError('index_list failed')
        else:
            raise ValueError('index_list failed: {}'.format(self.index_list))


    def init_files(self):
        """
        default files for output
        bam, sam, log, unmap, ...
        """
        self.config_toml = self.outdir + '/config.toml'
        check_path(self.outdir)


class AlignRnConfig(object):
    """
    Config files for AlignRn

    fq1=1, fq2=1, index_list=list

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': None,
            'index_list': None,
            'index_name': None,
            'threads': 1,
            'overwrite': False,
            'keep_tmp': False
        }
        self = update_obj(self, args_init, force=False)
        self.alignment_type = 'alignment_rn'

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # default files
        self.init_files()


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


    def init_index(self):
        """
        alignment index_list = {genome_index}
        
        The index_list are all valid, for the aligner
        """
        if isinstance(self.index_list, list):
            chk = [AlignIndex(aligner=self.aligner).is_index(i) \
                for i in self.index_list]
            if not all(chk):
                msg = '\n'.join([
                    '{:>6s} : {}'.format(str(j), k) for j, k in zip(chk, 
                        self.index_list)])
                print(msg)
                raise ValueError('index_list failed')

            # update index name
            self.index_name = [AlignIndex(index=i).index_name() \
                for i in self.index_list]
            self.index_name = ['{:02d}_{}'.format(i+1, n) for i, n in \
                enumerate(self.index_name)]
        else:
            raise ValueError('index_list failed: {}'.format(self.index_list))


    def init_files(self):
        """
        default files for output
        bam, sam, log, unmap, ...
        """
        # subdir
        subdir = os.path.join(self.outdir, self.smp_name)
        check_path(subdir)

        self.config_toml = subdir + '/config.toml'


class AlignR1Config(object):
    """
    Config files for AlignR1
    
    force: bowtie2, for rRNA

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """Required arguments, files, for Alignment()
        default values
        """
        args_init = {
            'smp_name': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': None,
            'index': None,
            'index_name': None,
            'threads': 1,
            'overwrite': False,
            'keep_tmp': False,
            'genome_size': 0,
            'genome_size_file': None,
        }
        self = update_obj(self, args_init, force=False)
        self.alignment_type = 'alignment_r1'

        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = getattr(self, 'smp_name', None)
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one


        # fq files
        self.init_fq()

        # alignment index
        self.init_index()

        # default files
        self.init_files()


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

            ai = AlignIndex(index=self.index)

            if self.genome_size < 1:
                self.genome_size = ai.index_size()

            if not isinstance(self.genome_size_file, str): 
                self.genome_size_file = ai.index_size(return_file=True)

        else:
            raise ValueError('index failed: {}'.format(self.index))

        if not AlignIndex(aligner=self.aligner, index=self.index).is_index():
            raise ValueError('not a {} index: {}'.format(
                self.aligner, self.index))


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
            'cmd_shell': os.path.join(subdir, 'cmd.txt'),
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


class AlignRpConfig(object):
    """
    Config files for AlignRp

    Example:
    1. 

    2. 

    3. 
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)


class AlignReader(object):
    """
    Read config file
    """
    def __init__(self, x):
        self.x = x
        self.read()
        self.alignment_type = self.args.get('alignment_type', None)

        self.is_alignment_r1 = self.alignment_type == 'alignment_r1'
        self.is_alignment_rn = self.alignment_type == 'alignment_rn'
        self.is_alignment_rx = self.alignment_type == 'alignment_rx'
        self.is_alignment_rt = self.alignment_type == 'alignment_rt'
        self.is_alignment_rd = self.alignment_type == 'build_design'


        if isinstance(self.alignment_type, str):
            self.is_hiseq = self.alignment_type.startswith('alignment_')
        else:
            self.is_hiseq = False


    def read(self):
        """
        # hiseq
        # add support for TOML
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
        p3x = os.path.join(self.x, 'config', 'config.toml')
        p4x = os.path.join(self.x, '*', '*', 'config.toml')
        p1 = glob.glob(p1x)
        p2 = glob.glob(p2x)
        p3 = glob.glob(p3x)
        p4 = glob.glob(p4x)

        # read config
        if len(p1) == 1:
            self.args = pickle2dict(p1[0])
        elif len(p2) == 1:
            self.args = pickle2dict(p2[0])
        elif len(p3) == 1:
            self.args = pickle2dict(p3[0])
        elif len(p4) == 1:
            self.args = pickle2dict(p4[0])        
        else:
            self.args = {}


############################################################
## aligner                                                ##
## base-level: aligner
## align, report in json
## 
## input:
##   - fq1: str (required)
##   - fq2: str or None
##   - smp_name: str or None (auto)
##   - index: str (required)
##   - index_name: str or None (auto)
##
##   - n_map
##   - unique_only
##   - extra_para
##   - parallel_jobs
##   - threads
##
## output:
##   - bam: str
##   - unmap1: str
##   - unmap2: str or None 
##
## return bam, unmap1, unmap2
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


class Bowtie(object):
    """
    Run bowtie for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.prep_cmd()


    def init_args(self):
        """
        check
        """
        args_local = AlignerConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'bowtie' # force changed.
        Toml(self.__dict__).to_toml(self.config_toml)        


    def prep_cmd(self):
        """
        unique
        n_map
        extra_para
        """
        # unique/multiple
        if self.extra_para is None:
            self.extra_para = ''

        if self.n_map < 1:
            self.n_map = 1

        if self.unique_only:
            arg_unique = '-m 1'
        else:
            arg_unique = '-v 2 -k {}'.format(self.n_map)

        # file type
        arg_fx = '-f' if self.fq_format == 'fasta' else '-q'

        # input file
        arg_input = '{}'.format(self.fq1) if not self.is_paired \
            else '-1 {} -2 {}'.format(self.fq1, self.fq2)

        self.cmd = ' '.join([
            '{}'.format(shutil.which('bowtie')),
            '--mm --best --sam --no-unal',
            arg_unique, 
            arg_fx,
            self.extra_para,
            '--un {}'.format(self.unmap),
            '{}'.format(self.index),
            arg_input,
            '2> {}'.format(self.align_log),
            '&& samtools view -@ {}'.format(self.threads),
            '-Sub -F 0x4 -',
            '| samtools sort -@ {}'.format(self.threads),
            '-o {} -'.format(self.bam),
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(
                self.bam, 
                self.align_flagstat),
        ])

        # save cmd
        with open(self.cmd_shell, 'wt') as w:
            w.write(self.cmd + '\n')


    def parse_align(self, to_toml=True):
        """
        Wrapper bowtie log

        Bowtie:
        # reads processed: 10000
        # reads with at least one reported alignment: 3332 (33.32%)
        # reads that failed to align: 457 (4.57%)
        # reads with alignments suppressed due to -m: 6211 (62.11%)

        or:

        # reads processed: 10000
        # reads with at least one reported alignment: 9543 (95.43%)
        # reads that failed to align: 457 (4.57%)

        unique, multiple, unmap, map, total

        skip: Warning, ...
        """
        dd = {}
        with open(self.align_log) as r:
            for line in r:
                # if not ':' in line or line.startswith('Warning'):
                #     continue
                if not line.startswith('#'):
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = eval(value)
                if 'reads processed' in line:
                    dd['total'] = value
                elif 'at least one reported alignment' in line:
                    dd['map'] = value
                elif 'failed to align' in line:
                    dd['unmap'] = value
                elif 'alignments suppressed due to -m' in line:
                    dd['multiple'] = value
                else:
                    pass

        # unique_only
        dd['unique'] = dd['map']
        dd['multiple'] = dd.get('multiple', 0) # default 0
        
        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']

        # save fqname, indexname,
        dd['fqname'] = self.smp_name
        dd['index_name'] = self.index_name

        # # sort by keys
        self.log_dict = dd

        # save dict to plaintext file
        with open(self.align_stat, 'wt') as w:
            ## version-1
            # for k, v in sorted(dd.items()):
            #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')

            # ## version-2
            # w.write('#') # header line
            # w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            # w.write('\t'.join(list(map(str, dd.values()))) + '\n')

            ## version-3
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 
                'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        if to_toml:
            Toml(dd).to_toml(self.align_toml)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))

        # rename unmap files
        unmap1 = self.subdir + '/' + self.smp_name + '.unmap_1.fastq'
        unmap2 = self.subdir + '/' + self.smp_name + '.unmap_2.fastq'

        if self.is_paired:
            # file_symlink(unmap1, self.unmap1)
            # file_symlink(unmap2, self.unmap2)
            file_copy(unmap1, self.unmap1)
            file_copy(unmap2, self.unmap2)
            file_remove([unmap1, unmap2], ask=False)
        else:
            self.unmap1 = self.unmap
            self.unmap2 = None

        # parse log file
        if file_exists(self.align_log):
            self.parse_align(to_toml=True)

        # temp files
        del_list = [self.sam, self.unmap1, self.unmap2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)

        return (self.bam, self.unmap1, self.unmap2)


class Bowtie2(object):
    """
    Run bowtie2 for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.prep_cmd()


    def init_args(self):
        """
        check
        """
        args_local = AlignerConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'bowtie2'
        Toml(self.__dict__).to_toml(self.config_toml)        


    def prep_cmd(self):
        """
        unique： samtools view -q 30

        Bowtie2 unique reads
        filt by tag: YT:Z:CP
        
        YT:Z: String representing alignment type
        CP: Concordant; DP: Discordant; UP: Unpaired Mate; UU: Unpaired.
        
        see1: https://www.biostars.org/p/19283/#19292
        see2: https://www.biostars.org/p/133094/#133127
         
        -F 2048: suppress Supplementary alignments
        """ 
        # nmap
        if self.n_map < 1: 
            self.n_map = 1 # 

        # extra_para
        if self.extra_para is None:
            self.extra_para = ''

        # nmap
        arg_nmap = '-k {}'.format(self.n_map)

        # file type
        arg_fx = '-f' if self.fq_format == 'fasta' else '-q'

        # input
        if self.is_paired:
            arg_io = '--un-conc {} -1 {} -2 {}'.format(
                self.unmap,
                self.fq1,
                self.fq2
            )
        else:
            arg_io = '--un {} -U {}'.format(
                self.unmap,
                self.fq1
            )

        # Bowtie2 unique reads
        # filt by tag: YT:Z:CP
        # 
        # YT:Z: String representing alignment type
        # CP: Concordant; DP: Discordant; UP: Unpaired Mate; UU: Unpaired.
        #
        # see1: https://www.biostars.org/p/19283/#19292
        # see2: https://www.biostars.org/p/133094/#133127
        #  
        # -F 2048: suppress Supplementary alignments
        if self.unique_only:
            if self.is_paired:
                arg_unique = ' '.join([
                    '&& samtools view -Sb',
                    '<(samtools view -H {} ;'.format(self.sam),
                    'samtools view -F 2048 {}'.format(self.sam),
                    "| grep 'YT:Z:CP')",
                    '| samtools sort -o {} -'.format(self.bam)])
            else:
                arg_unique = ' '.join([
                    '&& samtools view -Sb -F 2048 {}'.format(self.sam),
                    '| samtools sort -o {} -'.format(self.bam)])
        else:
            arg_unique = ' '.join([
                    '&& samtools view -Sb -F 2048 {}'.format(self.sam),
                    '| samtools sort -o {} -'.format(self.bam)])

        # command-line
        self.cmd = ' '.join([
            '{}'.format(shutil.which('bowtie2')),
            '--mm -p {}'.format(self.threads),
            '--local --very-sensitive --no-unal --no-mixed --no-discordant',
            arg_nmap,
            arg_io,
            '-x {}'.format(self.index),            
            '1> {} 2> {}'.format(self.sam, self.align_log),
            arg_unique,
            '&& samtools index {}'.format(self.bam),
            '&& samtools flagstat {} > {}'.format(self.bam, self.align_flagstat)
        ])

        # save cmd
        with open(self.cmd_shell, 'wt') as w:
            w.write(self.cmd + '\n')


    def parse_align(self, to_toml=True):
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

        # save dict
        # sort by keys
        # dd = dict(sorted(dd.items(), key=lambda kv: kv[1], reverse=True))
        self.log_dict = dd

        # save dict to plaintext file
        with open(self.align_stat, 'wt') as w:
            ## version-1
            # for k, v in sorted(dd.items()):
            #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')

            # ## version-2
            # w.write('#') # header line
            # w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            # w.write('\t'.join(list(map(str, dd.values()))) + '\n')

            ## version-3
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 
                'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        ## save to json
        if to_toml:
            Toml(dd).to_toml(self.align_toml)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))

        # parse log file
        if file_exists(self.align_log):
            self.parse_align(to_toml=True)

        # unmap files
        if not self.is_paired:
            self.unmap1 = self.unmap
            self.unmap2 = None

        # temp SAM
        file_remove(self.sam, ask=False)
        # temp fastq files
        del_list = [self.unmap1, self.unmap2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)

        return (self.bam, self.unmap1, self.unmap2)


## caution:
## --genomeLoad: LoadAndRemove (NoSharedMemory)
class STAR(object):
    """
    Run bowtie2 for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.prep_cmd()


    def init_args(self):
        """
        check
        """
        args_local = AlignerConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True) # update
        self.aligner = 'STAR'
        Toml(self.__dict__).to_toml(self.config_toml)

        ## prefix for files
        self.align_prefix = os.path.join(self.subdir, self.smp_name)

        ## check genome size (chrLength.txt)
        self.small_genome = self.genome_size < 10000000 # 10M

        ## genomeLoad
        genomeLoad = ['NoSharedMemory', 'LoadAndKeep', 'LoadAndRemove', 
            'LoadAndExit', 'Remove', 'NoSharedMemory']
        if not self.genomeLoad in genomeLoad:
            msg = '\n'.join([
                'unknown --genomeLoad: {}'.format(self.genomeLoad),
                'expected: {}'.format(' '.join(gl)),
                'auto switch to: NoSharedMemory'
                ])
            log.warning(msg)
            self.genomeLoad = 'NoSharedMemory'


    def prep_cmd(self):
        """
        ## for small genome
        mapping takes too long,
        99% of the reads not mapped
        change --seedPerWindowNmax
        max number of seeds per window
        https://github.com/alexdobin/STAR/issues/329
        https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ
        """
        if self.small_genome:
            log.warning('STAR on small genome (<10 Mb): {}'.format(
                self.index))
            seed_max = 5 # even smaller
        else:
            seed_max = 50 # default

        ## for unique map
        ## --outFilterMultimapNmax
        ## maximum number of loci the read is allowed to map to. [default: 10]
        ##
        ## or filt unique reads by: samtools view -q 255
        # n_map = 1 if self.unique_only else n_map
        if self.unique_only: 
            self.n_map = 1
        arg_unique = ' '.join([
            '--outFilterMultimapNmax {}'.format(self.n_map),
            '--seedPerWindowNmax {}'.format(seed_max)])

        ## STAR, .gz file input
        arg_reader = 'zcat' if self.fq1.endswith('.gz') else '-'

        ## fq2
        if self.fq2 is None: self.fq2 = '' # empty

        ## extra para
        arg_extra_para = '' if self.extra_para is None else self.extra_para

        ## For sharing memory in STAR
        ## by Devon Ryan: https://www.biostars.org/p/260069/#260077 
        ## by Dobin: https://github.com/alexdobin/STAR/pull/26
        ##
        ##  --genomeLoad
        ## 
        ##  NoSharedMemory: each job use its own copy
        ##  LoadAndExit:   load genome to memory, do not run alignment
        ##  Remove:        remove genome from memory, do not run alignment
        ##  LoadAndRemove: load genome to memory, remove genome after alignment
        ##  LoadAndKeep:   load genome to memory, keep it after run
        ##
        ##  pratice: for general usage
        ##  1. NoSharedMemory: (?) what if other jobs using it?
        ##
        ##  pratice: for general usage
        ##  1. LoadAndRemove: (?) 
        ##
        ##  pratice: for multiple alignment together.
        ##  1. LoadAndExit
        ##  2. (loop over samples): LoadAndKeep
        ##  3. Remove
        ##    
        self.cmd = ' '.join([
            '{}'.format(shutil.which('STAR')),
            '--genomeLoad {}'.format(self.genomeLoad), #
            '--runMode alignReads',
            '--genomeDir {}'.format(self.index),
            '--readFilesIn {} {}'.format(self.fq1, self.fq2),
            '--readFilesCommand {}'.format(arg_reader),
            '--outFileNamePrefix {}'.format(self.align_prefix),
            '--runThreadN {}'.format(self.threads),
            '--limitBAMsortRAM 10000000000',
            '--outSAMtype BAM SortedByCoordinate',
            '--outFilterMismatchNoverLmax 0.07',
            '--seedSearchStartLmax 20',
            '--outReadsUnmapped Fastx', # self.unmap1,
            arg_unique,
            arg_extra_para])

        # save cmd
        with open(self.cmd_shell, 'wt') as w:
            w.write(self.cmd + '\n')


    def update_files(self):
        """
        Update the filenames of STAR output
        bam: *Aligned.sortedByCoord.out.bam -> *.bam # mv
        log: *Log.final.out -> *.log # copy
        log: *Log.out -> *.out
        unmap: *Unmapped.out.mate1 -> *.unmap.1.fastq
               *Unmapped.out.mate1 -> *.unmap.1.fastq
        """
        ## default output of STAR
        ## bam, log, unmap files
        bam_from = self.align_prefix + 'Aligned.sortedByCoord.out.bam'
        log_from = self.align_prefix + 'Log.final.out'
        unmap1 = self.align_prefix + 'Unmapped.out.mate1'
        unmap2 = self.align_prefix + 'Unmapped.out.mate2'

        # new files
        file_symlink(bam_from, self.bam)
        file_symlink(log_from, self.align_log)
        if self.is_paired:
            file_symlink(unmap1, self.unmap1)
            file_symlink(unmap2, self.unmap2)
        else:
            file_symlink(unmap1, self.unmap1)
            self.unmap2 = None

        # # remove old files
        # del_list = [bam_from, unmap1, unmap2]
        # file_remove(del_list, ask=False)


    def parse_align(self, to_toml=True):
        """
        Wrapper

        STAR:
        *final.Log.out, (changed to *.log, in this script)

                                     Started job on |       Sep 12 11:08:57
                                 Started mapping on |       Sep 12 11:11:27
                                        Finished on |       Sep 12 11:11:29
           Mapping speed, Million of reads per hour |       18.00

                              Number of input reads |       10000
                          Average input read length |       73
                                        UNIQUE READS:
                       Uniquely mapped reads number |       47
                            Uniquely mapped reads % |       0.47%
                              Average mapped length |       51.66
                           Number of splices: Total |       5
                Number of splices: Annotated (sjdb) |       0
                           Number of splices: GT/AG |       3
                           Number of splices: GC/AG |       0
                           Number of splices: AT/AC |       0
                   Number of splices: Non-canonical |       2
                          Mismatch rate per base, % |       2.14%
                             Deletion rate per base |       0.04%
                            Deletion average length |       1.00
                            Insertion rate per base |       0.00%
                           Insertion average length |       0.00
                                 MULTI-MAPPING READS:
            Number of reads mapped to multiple loci |       83
                 % of reads mapped to multiple loci |       0.83%
            Number of reads mapped to too many loci |       19
                 % of reads mapped to too many loci |       0.19%
                                      UNMAPPED READS:
           % of reads unmapped: too many mismatches |       0.02%
                     % of reads unmapped: too short |       98.31%
                         % of reads unmapped: other |       0.18%
                                      CHIMERIC READS:
                           Number of chimeric reads |       0
                                % of chimeric reads |       0.00%

        unique, multiple, unmap, map, total
        total unique multiple map unmap fqname index
        """
        dd = {}
        with open(self.align_log, 'rt') as r:
            for line in r:
                value = line.strip().split('|')
                if not len(value) == 2:
                    continue
                value = value[1].strip()
                if 'Number of input reads' in line:
                    dd['total'] = int(value)
                elif 'Uniquely mapped reads number' in line:
                    dd['unique'] = int(value)
                elif 'Number of reads mapped to multiple loci' in line:
                    dd['multiple'] = int(value)
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

        # sort by keys
        # dd = dict(sorted(d.items(), key=lambda kv: kv[1], reverse=True))
        self.log_dict = dd

        # save dict to plaintext file
        # fixed order
        # total unique multiple map unmap fqname index
        with open(self.align_stat, 'wt') as w:
            ## version-1
            # for k, v in sorted(dd.items()):
            #     w.write('\t'.join([self.config.fqname, self.config.index_name, k, str(v)]) + '\n')
            
            # ## version-2
            # w.write('#') # header line
            # w.write('\t'.join(list(map(str, dd.keys()))) + '\n')
            # w.write('\t'.join(list(map(str, dd.values()))) + '\n')
            
            ## version-3
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 
                'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        ## save to json
        if to_toml:
            Toml(dd).to_toml(self.align_toml)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        # run cmd
        if file_exists(self.bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(self.bam))
        else:
            try:
                run_shell_cmd(self.cmd)
            except:
                log.error('align() failed, check {}'.format(self.align_log))

        # rename files
        self.update_files()

        # parse log file
        if file_exists(self.align_log):
            self.parse_align(to_toml=True)

        # unmap files
        if not self.is_paired:
            self.unmap1 = self.unmap
            self.unmap2 = None

        # # temp fastq files
        # del_list = [self.unmap1, self.unmap2]
        # if not self.keep_tmp:
        #    file_remove(del_list, ask=False)

        return (self.bam, self.unmap1, self.unmap2)


## to-do
class BWA(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        pass


    def init_args(self):
        pass
        

    def run(self):
        pass


class Hisat2(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        pass

    def init_args(self):
        pass
        

    def run(self):
        pass


class Kallisto(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        pass


    def init_args(self):
        pass
        

    def run(self):
        pass


class Salmon(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        pass

    def init_args(self):
        pass
        

    def run(self):
        pass


## for index
## to-do
##   - build index (not recommended)
class AlignIndex(object):

    def __init__(self, **kwargs):
        """
        Two keywords: index, aligner
        Required args:
          - aligner
          - index (optional)
          - genome
          - group : genome, rRNA, transposon, piRNA_cluster, ...
          - genome_path
        """
        self.update(kwargs)
        self.init_args()
        # self.name = self.index_name()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        force: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        force exists attr
        """
        # fresh start
        if remove is True:
            for k in self.__dict__:
                # self.__delattr__(k)
                delattr(self, k)
        # add attributes
        if isinstance(d, dict):
            for k, v in d.items():
                if not hasattr(self, k) or force:
                    setattr(self, k, v)


    def init_args(self):
        args_default = {
            'index': None,
            'aligner': None
        }
        self.update(args_default, force=False)
        # update: remove `genome` from object
        if hasattr(self, 'genome'):
            delattr(self, 'genome')


    def get_aligner(self, index=None):
        """
        Search the available index for aligner:
        bowtie, [*.[1234].ebwt,  *.rev.[12].ebwt]
        bowtie2, [*.[1234].bt2, *.rev.[12].bt2]
        STAR,
        bwa,
        hisat2,
        """
        # unknown
        if index is None:
            index = self.index

        if index is None: # required
            log.warning('AlignIndex(index=), required for guessing the aligner')
            return None

        # check
        bowtie_files = ['{}.{}.ebwt'.format(index, i) for i in range(1, 5)]
        bowtie2_files = ['{}.{}.bt2'.format(index, i) for i in range(1, 5)]
        hisat2_files = ['{}.{}.ht2'.format(index, i) for i in range(1, 4)]
        bwa_files = ['{}.{}'.format(index, i) for i in ['sa', 'amb', 'ann', 'pac', 'bwt']]
        star_files = [os.path.join(index, i) for i in [
            'SAindex',
            'Genome',
            'SA',
            'chrLength.txt',
            'chrNameLength.txt',
            'chrName.txt',
            'chrStart.txt',
            'genomeParameters.txt']]

        ## check
        chk0 = all(file_exists(bowtie_files))
        chk1 = all(file_exists(bowtie2_files))
        chk2 = all(file_exists(hisat2_files))
        chk3 = all(file_exists(bwa_files))
        chk4 = all(file_exists(star_files))

        ## check file exists
        if chk0:
            aligner = 'bowtie'
        elif chk1:
            aligner = 'bowtie2'
        elif chk2:
            aligner = 'hisat2'
        elif chk3:
            aligner = 'bwa'
        elif chk4:
            aligner = 'star' # STAR
        else:
            aligner = None

        return aligner


    def is_index(self, index=None):
        """
        guesses the aligner, from index
        """
        if index is None:
            index = self.index
        
        ## return the aligner, from index
        if self.aligner is None:
            chk0 = not self.get_aligner(index=index) is None # 
        else:
            chk0 = self.aligner.lower() == self.get_aligner(index=index)

        return chk0


    def search(self, **kwargs):
        """
        Search the index for aligner: 
        STAR, bowtie, bowtie2, bwa, hisat2
        para:

        *genome*    The ucsc name of the genome, dm3, dm6, mm9, mm10, hg19, hg38, ...
        *group*      Choose from: genome, rRNA, transposon, piRNA_cluster, ...

        structure of genome_path:
        default: {HOME}/data/genome/{genome_version}/{aligner}/


        ## bowtie/bowtie2/hisat2/...
        path-to-genome/
            |- Bowtie_index /
                |- genome
                |- rRNA
                |- MT_trRNA
                |- transposon
                |- piRNA_cluster

        ## STAR
        path-to-genome/
            |- Bowtie_index /
                |- genome/
                |- rRNA/
                |- MT_trRNA/
                |- transposon/
                |- piRNA_cluster/
        """
        self.update(kwargs, force=True) # input args

        args_default = {
            'genome': None,
            'group': None,
            'genome_path': os.path.join(str(pathlib.Path.home()), 'data', 'genome'),
        }
        self.update(args_default, force=False) # assign default values

        ## required arguments: aligner
        aligner_supported = ['bowtie', 'bowtie2', 'STAR', 'hisat2', 'bwa', 
                             'kallisto', 'salmon']
        if not self.aligner in aligner_supported:
            log.error('AlignIndex(aligner=) required, candidate: {}'.format(aligner_supported))
            return None

        ## required arguments: genome
        if self.genome is None:
            log.error('AlignIndex().search(), require, genome=.')
            return None

        ## required arguments: group
        group_list = ['genome', 'genome_rm', 'MT_trRNA', 'rRNA', 'chrM', 
                      'structural_RNA', 'transposon', 'te', 'piRNA_cluster', 
                      'miRNA', 'miRNA_hairpin']
        if not self.group in group_list:
            log.error('AlignIndex().search(group={}) unknown, expect {}'.format(self.group, group_list))
            return None

        ## create index path
        p0 = os.path.join(self.genome_path, self.genome, self.aligner + '_index') # [case sensitive] STAR bowtie
        # p1 = [os.path.join(p0, i) for i in self.group_list]
        p1 = os.path.join(p0, self.group)

        if self.is_index(index=p1) and self.get_aligner(index=p1) == self.aligner.lower():
            return p1
        else:
            log.warning('index not exists: {}'.format(p1))
            return None


    def index_name(self, index=None):
        """
        Get the name of index
        basename: bowtie, bowtie2, hisqt2, bwa
        folder: STAR
        """
        if index is None:
            index = self.index
        
        ## check
        if index is None:
            log.warning('AlignIndex(index=) or AlignIndex().index_name(index=) required')
            return None

        if not self.is_index(index=index):
            log_msg = '\n'.join([
                'index not exists, or not match the aligner:',
                '{:>30s}: {}'.format('Index', index),
                '{:>30s}: {}'.format('Aligner expected', self.get_aligner(index=index)),
                '{:>30s}: {}'.format('Aligner get', self.aligner)])
            log.warning(log_msg)
            return None

        if file_exists(index):
            # STAR
            return os.path.basename(index)
        elif os.path.basename(index) == 'genome':
            # ~/data/genome/dm3/bowtie2_index/genome
            # bowtie, bowtie2, bwa, hisat2
            # iname = os.path.basename(index)
            return os.path.basename(os.path.dirname(os.path.dirname(index)))
        else:
            # other groups
            return os.path.basename(index)


    def _tmp(self, is_dir=False, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='tmp')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
                delete=False)
        return tmp.name


    def index_size(self, index=None, return_file=False):
        """
        chr size of index

        bowtie:  bowtie-inspect -s <index>
        bowtie2: bowtie-inspect -s <index>
        STAR:  chrNameLength.txt
        """
        if index is None:
            index = self.index

        ## check
        if index is None:
            log.warning('AlignIndex(index=) or AlignIndex().index_name(index=) required')
            return None

        if not self.is_index(index=index):
            log_msg = '\n'.join([
                'index not exists, or not match the aligner:',
                '{:>30s}: {}'.format('Index', index),
                '{:>30s}: {}'.format('Aligner expected', self.get_aligner(index=index)),
                '{:>30s}: {}'.format('Aligner get', self.aligner)])
            log.warning(log_msg)
            return None

        ## aligner
        gsize = self._tmp(suffix='.chrom.sizes')
        chrLength = 0
        aligner = self.get_aligner(index).lower()

        if aligner in ['bowtie', 'bowtie2', 'hisat2', 'star']:
            # get genome size
            if aligner.lower() == 'star':
                gsize = os.path.join(index, 'chrNameLength.txt')
            else:
                if aligner == 'bowtie':
                    x_inspect = shutil.which('bowtie-inspect')
                elif aligner == 'bowtie2':
                    x_inspect = shutil.which('bowtie2-inspect')
                elif aligner == 'hisat2':
                    x_inspect = shutil.which('hisat2-inspect')
                else:
                    pass

                # inspect
                cmd = ' '.join([
                    '{}'.format(x_inspect),
                    '-s {} |'.format(index),
                    'grep ^Sequence |',
                    "sed -E 's/^Sequence-[0-9]+\t//' > {}".format(gsize)])

                # run
                try:
                    os.system(cmd)
                except:
                    log.error('failed to run: {}'.format(x_inspect))

            # read size
            with open(gsize, 'rt') as r:
                s = [i.strip().split('\t')[-1] for i in r.readlines()]
            chrLength = sum(map(int, s))

        else:
            log.error('unknown aligner: {}'.format(aligner))

        ## output
        return gsize if return_file else chrLength



