# -*- coding:utf-8 -*-

"""
###############
1. unmap files:
  - bowtie  : _1.fastq, _2.fastq
  - bowtie2 : .1.fastq, .2.fastq
  - STAR    : .unmap.1.fastq
###############



AlignConfig() : init arguments, determine alignment type, tools, idnex, ...

AlignFqSingle() : 1 fastq to N index

AlignFqMultiple() : N fastq to N index

AlignFqSingleIndexSingle() : 1 fastq to 1 index

AlignFqSingleIndexMultiple() : 1 fastq to N index

Align() : wrap all


## Structure:
Align()
    |- AlignConfig()
        |- AlignFqMultiple()
                |- AlignFqSingle()
                        |- AlignFqSingleIndexMultiple()
                                |- AlignFqSingleIndexSingle ()
                                         |- Aligners()

## utils
AlignIndex()



Align fasta/q files to reference sequence

bowtie
bowtie2
STAR
bwa
hisat2
...


arguments:

fq
index
outdir

# optional

unique
mismatch
insertion (bowtie2 -X)

max-multi

# output
sam

unmap

# custom args
arg_str


## standard

input: fq(s), index, args
output: sam, unmap, log

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
import copy # to copy objects
import pathlib
import shutil
import logging
import json
import pandas as pd
from multiprocessing import Pool
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions


def print_df(d):
    if isinstance(d, dict):
        for k, v in sorted(d.items(), key=lambda item: item[0]):  # 0:key, 1:value
            print('{:>15} : {}'.format(k, v))
    else:
        print(d)


@Logger('INFO')
class Alignment(object):
    """
    The wrapper for alignment
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        local_config = AlignConfig(**self.__dict__)
        self.update(local_config.__dict__, force=True) # update
        self.pickle = None # remove


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


    def run(self):
        """
        Run all Alignment
        """
        args = self.__dict__.copy() # local

        if self.align_type == [2, 2]:
            AlignFqNIndexN(**args).run()
        elif self.align_type == [2, 1]:
            AlignFqNIndex1(**args).run()
        elif self.align_type == [1, 2]:
            AlignFq1IndexN(**args).run()
        elif self.align_type == [1, 1]:
            AlignFq1Index1(**args).run()
        else:
            pass # no except



############################################################
## Alignment                                              ##
## top-level for alignment
## 
## input:
##   - pickle: str or None
##   - fq1: list (required)
##   - fq2: list or None
##   - aligner: str (required)
##   - genome: str or None
##   - spikein: str or None
##   - align-to-rRNA: bool
##   - align-to-chrM: bool
##   - smp_name: list or None (auto)
##   - index_list: list or None
##   - index_name: list or None (auto)
##
##   - index_list_equal: bool
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
##
## priority:
##   - 1. pickle
##   - 2. genome, spikein, align-to-chrM, ...
##   - x. extra_index
##
##
## return bam, unmap1, unmap2

class AlignConfig(object):
    """
    The wrapper for alignment, init arguments, top-level:
    
    Determine the align_type:
    [2, 2] : N fastq, N index
    [1, 2] : 1 fastq, N index
    [1, 1] : 1 fastq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args() # all args
        self.align_type = self.mission()


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
        """
        default arguments
        required:
          - fq1
          - genome/ext-index
          - outdir
        
        optional:
          - aligner


        prepare index
          - align-to-rRNA
          - align-to-chrM
        """
        args_default = {
            'fq1': None,
            'fq2': None,
            'index_list': None,
            'outdir': str(pathlib.Path.cwd()),
            'aligner': 'bowtie',
            'genome': None,
            'spikein': None,
            'index_name': None,
            'extra_index': None,
            'index_list_equal': False,
            'align_parallel': False,
            'align_to_chrM': False,
            'align_to_rRNA': False,
            'align_to_MT_trRNA': False,
            'repeat_masked': False,
            'extra_para': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'pickle': None,
            'genomeLoad': 'NoSharedMemory',
            }
        self.update(args_default, force=False) # update missing attrs
        
        # 1st level: pickle / update all config
        # fresh start by pickle
        if not self.pickle is None:
            if os.path.exists(self.pickle) and self.pickle.endswith('.pickle'):
                args_pickle = pickle_to_dict(self.pickle)
                args_pickle.pop('pickle', None) # empty, for further round
                self.update(args_pickle, force=True) # fresh start

        ## outdir
        self.outdir = file_abspath(self.outdir)

        ## spikein
        if self.spikein == self.genome:
            self.spikein = None #

        ## args
        self.init_fq()
        self.init_index()
        self.init_cpu()

        ## args
        self.is_paired = self.check_paired() # is_paired ?


    def init_fq(self):
        """
        Make sure:
        fq1: str or list, convert to list
        fq2: None or str or list, convert to list
        exists
        """
        ####################################################
        ## 1. for fastq1 files                              #
        if isinstance(self.fq1, list):
            self.fq1 = file_abspath(self.fq1)
        elif isinstance(self.fq1, str):
            self.fq1 = [file_abspath(self.fq1)]
        else:
            log.error('failed fq1, str, list expected, {} found'.format(type(self.fq1).__name__))

        ####################################################
        ## 2. for fastq2 files                              #
        if self.fq2 is None:
            self.fq2 = None
        elif isinstance(self.fq2, list):
            self.fq2 = file_abspath(self.fq2)
        elif isinstance(self.fq2, str):
            self.fq2 = [file_abspath(self.fq2)]
        else:
            log.error('failed fq2, None, list expected, {} found'.format(type(self.fq2).__name__))

        ####################################################
        ## 3. for smp_name                                 #
        ## smp_name [from fq1, smp_path]
        if isinstance(self.smp_name, list):
            pass
        elif isinstance(self.smp_name, str):
            self.smp_name = [self.smp_name]
        elif self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True)
        else:
            raise Exception('failed, smp_name, None, str expected, {} found'.format(type(self.smp_name).__name__))

        ## fq type (fa or fq)
        tmp_fq1 = next(iter(self.fq1), None)
        self.fq_format = self.fq_type = Fastx(tmp_fq1).format # to-do !!!

        ## file exists
        chk0 = isinstance(self.fq1, list)
        chk1 = isinstance(self.fq2, list) or self.fq2 is None
        chk2 = isinstance(self.smp_name, list)
        chk3 = True # len(self.fq1) == len(self.fq2)
        chk4 = len(self.fq1) == len(self.smp_name)
        chk5 = True # all(file_exists(self.fq1))
        chk6 = True # all(file_exists(self.fq2, isfile=True)) or self.fq2 is None
        if not all([chk0, chk1, chk2, chk3, chk4, chk5, chk6]):
            raise Exception('failed, fq1, fq2: \nfq1: {}, \nfq2: {}, {}'.format(self.fq1, self.fq2, [chk0, chk1, chk2, chk3, chk4, chk5, chk6]))


    def init_index(self):
        """
        Make sure, index_list = list
        """
        self.index_list = self.get_index_list()

        ## check: is list
        if isinstance(self.index_list, list):
            pass
        else:
            raise Exception('index_list, list expected, {} found'.format(type(self.index_list).__name__))
        
        ## index correct
        chk0 = [AlignIndex(aligner=self.aligner, index=i).is_index() for i in self.index_list]
        log_msg0 = ''
        for a, b in zip(chk0, self.index_list):
            log_msg0 += '{}: {}\n'.format(a, b)
        if not all(chk0):
            raise Exception(log_msg0)

        ## index_name
        index_name_auto = [AlignIndex(index=i).index_name() for i in self.index_list]
        if self.index_name is None:
            self.index_name = index_name_auto
        elif isinstance(self.index_name, list):
            if not len(self.index_name) == len(self.index_list):
                log.warning('failed index_name, auto fetch from index_list')
                self.index_name = index_name_auto
        elif isinstance(self.index_name, str):
            log.warning('failed index_name, auto fetch from index_list')
            self.index_name = index_name_auto
        else:
            raise Exception('failed index_name, list expected, {} found'.format(type(self.index_name).__name__))

        ## check:
        chk1 = isinstance(self.index_list, list)
        chk2 = isinstance(self.index_name, list)
        chk3 = len(self.index_list) == len(self.index_name)
        ##
        log_msg1 = '\n'.join([
            'arguments failed: index_list',
            '{}: {:<20s}: {}'.format(chk0, 'list', self.index_list),
            '{}: {:<20s}: {}'.format(chk1, 'list', self.index_name),
            '{}: {:<20s}: {}, {}'.format(chk2, 'equal_number', len(self.index_list), len(self.index_name))])
        if not all([chk1, chk2, chk3]):
            raise Exception(log_msg1)


    def init_cpu(self):
        """
        threads, CPUs
        """
        ## check number of threads, parallel_jobs
        ## parallel jobs * threads
        n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

        max_jobs = int(n_cpu / 4.0)
        ## check parallel_jobs (max: 1/4 of n_cpus)
        if self.parallel_jobs > max_jobs: 
            log.warning('Too large, change parallel_jobs from {} to {}'.format(
                self.parallel_jobs, max_jobs))
            self.parallel_jobs = max_jobs

        ## check threads
        max_threads = int(0.8 * n_cpu / self.parallel_jobs)
        if self.threads * self.parallel_jobs > 0.8 * n_cpu:
            log.warning('Too large, change threads from {} to {}'.format(
                self.threads, max_threads))
            self.threads = max_threads        


    ## build index for alignment
    def get_index_list(self):
        """
        Create index list
        based on the arguments
        
        1st level: index_list (ordered, ext_index, TE/piRNA_cluster/consensus/repeat/...)
        2nd level: spikein (rRNA/chrM) | genome (rRNA/chrM) | extra_index
        """
        index_list = [] # init

        # 1st level: index_list, ...
        if isinstance(self.index_list, str):
            if AlignIndex(aligner=self.aligner).is_index(index=self.index_list):
                index_list += [self.index_list]
        elif isinstance(self.index_list, list):
            tmp1 = [AlignIndex(aligner=self.aligner).is_index(index=i) for i in self.index_list]
            index_list += self.index_list
        # 2nd level: [spikein | genome | extra_index] : auto generated, eaual-level
        else: 
            # 2nd level: spikein
            if isinstance(self.spikein, str):
                if self.spikein == self.genome:
                    self.spikein = None
                else:
                    # for spikein: choose tag
                    if self.align_to_MT_trRNA is True:
                        tag = 'MT_trRNA'            
                    elif self.align_to_rRNA is True:
                        tag = 'rRNA'
                    elif self.align_to_chrM is True:
                        tag = 'chrM'
                    else:
                        tag = None

                    # for group
                    index0 = None
                    if not tag is None:
                        index0 = AlignIndex(aligner=self.aligner).search(genome=self.spikein, group=tag)
                        if index0 is None:
                            log.warning('index {} for {}, not detected, skipped'.format(
                                tag, self.spikein))
                        else:
                            index_list.append(index0)

                    # for genome
                    index1 = AlignIndex(aligner=self.aligner).search(genome=self.spikein, group='genome')
                    if index1 is None:
                        log.warning('index {} for {}, not detected, skipped'.format(
                                tag, self.spikein))
                    else:
                        index_list.append(index1)

            # 2nd level: genome
            if isinstance(self.genome, str):
                # choose tag
                if self.align_to_MT_trRNA is True:
                    tag = 'MT_trRNA'            
                elif self.align_to_rRNA is True:
                    tag = 'rRNA'
                elif self.align_to_chrM is True:
                    tag = 'chrM'
                else:
                    tag = None

                # for group
                index0 = None
                if not tag is None:
                    index0 = AlignIndex(aligner=self.aligner).search(genome=self.genome, group=tag)
                    if index0 is None:
                        log.warning('index {} for {}, not detected, skipped'.format(
                            tag, self.spikein))
                    else:
                        index_list.append(index0)

                # for genome
                index1 = AlignIndex(aligner=self.aligner).search(genome=self.genome, group='genome')
                if index1 is None:
                    raise Exception('index {} for {}, not detected, skipped'.format(
                            tag, self.genome))
                else:
                    index_list.append(index1)

            # 2nd level: extra-index
            if isinstance(self.extra_index, list):
                tmp2 = [AlignIndex(aligner=self.aligner).is_index(i) for i in self.extra_index]
                index_list += self.extra_index

        ## check, index_list, valid
        chk0 = [AlignIndex(aligner=self.aligner, index=i).is_index() for i in index_list]
        ## message
        log_msg0 = 'Status of the index:\n'
        for i, j in zip(chk0, index_list):
            log_msg0 += '{}: {}\n'.format(i, j)

        if not all(chk0):
            raise Exception(log_msg0)

        ## output
        return index_list


    def check_paired(self):
        """
        Determine the type of Alignment:
        is_paired: bool
        is_multi_index: bool ?
        """
        # chk0 = len(self.fq1) == len(self.fq2)

        # if not chk0:
        #     raise Exception('fq1 and fq2 not match in numbers: \nfq1: {} \nfq2: {}'.format(
        #         len(self.fq1),
        #         len(self.fq2)))

        if self.fq2 is None:
            is_paired = False
        else:
            is_paired = True

        return is_paired

        # is_paired = []
        # check_log = 'fastq files status:\n'
        # chk0 = False # status, errors
        # for fq1, fq2 in zip(self.fq1, self.fq2):
        #     check_log += 'fq1: {}\tfq2: {}\n'.format(fq_name(fq1), fq_name(fq2))
        #     if fq2 is None:
        #         is_paired.append(False)
        #     elif os.path.exists(fq2):
        #         if fq_name(fq1, pe_fix=True) == fq_name(fq2, pe_fix=True):
        #             is_paired.append(True) # check file name
        #         else:
        #             chk0 = True
        
        # if chk0:
        #     raise Exception(check_log)

        # return all(is_paired)


    def mission(self):
        """
        Determine the align_type:
        [2, 2] : N fastq, N index
        [1, 2] : 1 fastq, N index
        [1, 1] : 1 fastq, 1 index
        """
        create_dirs = getattr(self, 'create_dirs', True)

        # fq1 number
        align_type = [2] if len(self.fq1) > 1 else [1] # single / multiple

        # index number
        align_type += [2] if len(self.index_list) > 1 else [1] # [2, 2], [2, 1], [1, 2], [1, 1]

        # choose sub-program
        if align_type == [2, 2]:
            self.init_fq_n_index_n(create_dirs)
        elif align_type == [1, 2]:
            self.init_fq_1_index_n(create_dirs)
        elif align_type == [1, 1]:
            self.init_fq_1_index_1(create_dirs)
        elif align_type == [2, 1]:
            self.init_fq_2_index_1(create_dirs)
        else:
            pass # no except

        return align_type


    def init_fq_n_index_n(self, create_dirs=True):
        """
        align_type: [2, 2] #fq, index
        for: multiple fq
        # parallel-jobs
        """
        chk0 = isinstance(self.fq1, list) and len(self.fq1) >= 1 # multiple
        chk1 = isinstance(self.smp_name, list) and len(self.smp_name) >= 1 # multiple
        ckk2 = isinstance(self.index_name, list) and len(self.index_name) >= 1 # multiple
        chk3 = isinstance(self.index_list, list) and len(self.index_list) >= 1 # multiple
        chk4 = len(self.fq1) == len(self.smp_name) # multiple
        chk5 = len(self.index_name) == len(self.index_list) # multiple

        #############################################################
        ## save file_path / default
        ## update these attributes, everytime
        ## auto-generated
        # project_dir = os.path.join(self.outdir, self.smp_name)
        project_dir = self.outdir
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': os.path.join(project_dir, 'config'),
            'report_dir': os.path.join(project_dir, 'report'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json')
        }
        self.update(auto_files, force=True) # key point

        ## create directories
        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


    ## deprecated, not used
    def init_fq_n_index_1(self, create_dirs=True):
        """
        ! question: skipped, replaced by [2, 2] ! to-do
        align_type: [2, 1] #fq, index
        for: multiple fq
        # parallel-jobs
        """
        chk0 = isinstance(self.fq1, list) and len(self.fq1) >= 1 # multiple
        chk1 = isinstance(self.smp_name, list) and len(self.smp_name) >= 1 # multiple
        chk2 = isinstance(self.index_list, list) and len(self.index_list) == 1 # multiple
        ckk3 = isinstance(self.index_name, list) and len(self.index_name) == 1 # multiple
        chk4 = len(self.fq1) == len(self.smp_name) # multiple
        chk5 = len(self.index_name) == len(self.index_list) # multiple

        index_name = next(iter(self.index_name), None)
        #############################################################
        ## save file_path / default
        ## update these attributes, everytime
        ## auto-generated
        project_dir = os.path.join(self.outdir, self.smp_name)
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'report_dir': os.path.join(project_dir, 'report'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json')
        }
        self.update(auto_files, force=True) # key point

        ## create directories
        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


    def init_fq_1_index_n(self, create_dirs=True):
        """
        align_type: [1, 2] #fq, index
        for: multiple index
        # parallel-jobs
        """
        chk0 = isinstance(self.fq1, list) and len(self.fq1) == 1 # single
        chk1 = isinstance(self.smp_name, list) and len(self.smp_name) == 1 # single
        ckk2 = isinstance(self.index_name, list) and len(self.index_name) >= 1 # multiple
        chk3 = isinstance(self.index_list, list) and len(self.index_list) >= 1 # multiple
        chk4 = len(self.index_name) == len(self.index_list) # multiple

        smp_name = next(iter(self.smp_name), None) # 1-st item
        #############################################################
        ## save file_path / default
        ## update these attributes, everytime
        ## auto-generated
        project_dir = os.path.join(self.outdir, smp_name)
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'report_dir': os.path.join(project_dir, 'report'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json'),
            'align_stat': os.path.join(project_dir + '.align.txt')
        }
        self.update(auto_files, force=True) # key point

        ## create directories
        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


    def init_fq_1_index_1(self, create_dirs=True):
        """
        align_type: [1, 1] #fq, index
        for: single index
        # single
        """
        chk0 = isinstance(self.fq1, list) and len(self.fq1) == 1 # single
        chk1 = isinstance(self.smp_name, list) and len(self.smp_name) == 1 # single
        chk2 = isinstance(self.index_name, list) and len(self.index_name) == 1 # single
        chk3 = isinstance(self.index_list, list) and len(self.index_list) == 1 # single
        if not all([chk0, chk1, chk2, chk3]):
            raise Exception('check, fq1, smp_name, index_list, index_name: single expected')

        smp_name = next(iter(self.smp_name), None) # 
        index_name = next(iter(self.index_name), None)
        #############################################################
        ## save file_path / default
        ## update these attributes, everytime
        ## auto-generated
        project_dir = os.path.join(self.outdir, smp_name, index_name)
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'data_dir': os.path.join(project_dir, 'data'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json'),
            'report_dir': os.path.join(project_dir, 'report'),
            'align_stat': os.path.join(self.outdir, index_name + '.align.txt' )
        }
        self.update(auto_files, force=True) # key point

        #############################################################
        ## default files
        self.align_log = os.path.join(project_dir, smp_name + '.log')
        self.align_json = os.path.join(project_dir, smp_name + '.json')
        self.align_sam = os.path.join(project_dir, smp_name + '.sam')
        self.align_bam = os.path.join(project_dir, smp_name + '.bam')
        if self.is_paired is True:
            self.unmap1 = os.path.join(project_dir, smp_name + '.unmap.1.' + self.fq_format)
            self.unmap2 = os.path.join(project_dir, smp_name + '.unmap.2.' + self.fq_format)
        else:
            self.unmap1 = os.path.join(project_dir, smp_name + '.unmap.' + self.fq_format)
            self.unmap2 = None

        ## create directories
        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])


class AlignFqNIndexN(object):
    """
    Align N fq to N index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args()


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
        """
        Initiate directories, config:
        save config files
        outdir/config/*json, *pickle, *txt
        """
        local_config = AlignConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        
        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def report(self):
        """
        Create Alignment report for N fq, N index
        """
        pass


    def run_fq_single(self, fq1):
        """
        run 1 fq on N index
        """
        obj_i = copy.copy(self)
        i = self.fq1.index(fq1) # index in fq1

        ## organize args for single run:
        fq1 = [fq1]
        fq2 = None if obj_i.fq2 is None else [obj_i.fq2[i]]
        smp_name = [obj_i.smp_name[i]]
        outdir = obj_i.outdir # os.path.join(obj_i.outdir, obj_i.smp_name[i])
        args_i = {
            'fq1': fq1,
            'fq2': fq2,
            'smp_name': smp_name, 
            'outdir': outdir,
        }

        # update: obj_i
        for k, v in args_i.items():
            setattr(obj_i, k, v) # force, update

        AlignFq1IndexN(**obj_i.__dict__).run()


    def run(self):
        """
        run in parallel
        """
        ## run in parallel !!!
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_fq_single, self.fq1)

        ## report
        self.report()

        # out
        return self.outdir


class AlignFqNIndex1(object):
    """
    Align 1 fq to N index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)


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
        pass
        

    def run(self):
        pass


class AlignFq1IndexN(object):
    """
    Align 1 fq to N index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args()


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
        local_config = AlignConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        
        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## add rank in index name
        self.index_name = ['{}.{}'.format(i + 1, v) for i, v in enumerate(self.index_name)]


    def save_stat(self):
        """
        Save alignment stat files of all index, into one file
        #
        in_files: smp_name.stat
        out_file: smp_name.align.txt

        #total  map     unmap   unique  multiple        fqname  index_name
        100000  2716    97284   2716    0       demo_control_rep1       1.rRNA
        """
        # save to
        outfile = self.align_stat 

        # prepare the target name
        smp_name = next(iter(self.smp_name), None) #
        stat_name = smp_name + '.stat' # target name in the directory
        stat_list = [] # init

        # sub-dirs:
        subdirs = listdir(self.project_dir, include_dir=True)
        stat_list = [os.path.join(i, stat_name) for i in subdirs if os.path.isdir(i)] # dir only
        stat_list = [i for i in stat_list if os.path.exists(i)] # file exists

        if len(stat_list) == 0:
            return None

        # read content
        with open(outfile, 'wt') as w:
            # save header
            f1 = next(iter(stat_list), None)
            with open(f1) as r:
                w.write(next(r)) # the header, 1-st line

            # save content
            for x in stat_list:
                with open(x) as r:
                    for line in r:
                        if line.startswith('#'): 
                            continue
                        w.write(line)

        ## out
        return outfile


    def run(self):
        """
        Run 1 fastq to N index
        """
        obj_a = copy.copy(self)

        # for idx in self.index_list:
        for i, idx in enumerate(self.index_list):
            obj_i = copy.copy(obj_a)

            ## update args for 1fq,1index
            outdir = obj_i.outdir # os.path.join(obj_i.outdir, obj_i.index_name[i])
            index_list = [idx]
            index_name = [obj_i.index_name[i]]
            args_i = {
                'index_list': index_list,
                'index_name': index_name,
                'outdir': outdir,
            }
            for k, v in args_i.items():
                setattr(obj_i, k, v) # force, updated

            _, unmap1, unmap2 = AlignFq1Index1(**obj_i.__dict__).run()
            # update unmapped reads, to fq1, fq2 in next round
            if not obj_a.index_list_equal:
                obj_a.fq1 = unmap1
                obj_a.fq2 = unmap2

        # organize stat
        self.save_stat()

        return self.outdir


class AlignFq1Index1(object):
    """
    Align 1 fq to 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args()


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
        local_config = AlignConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        
        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)


    # assign: Class
    def pick_ports(self):
        aligner = self.aligner.lower()

        ## check
        if aligner == 'bowtie':
            port = Bowtie
        elif aligner == 'bowtie2':
            port = Bowtie2
        elif aligner == 'star':
            port = STAR
        elif aligner == 'bwa':
            port = BWA
        elif aligner == 'hisat2':
            port = Hisat2
        elif aligner == 'kallisto':
            port = Kallisto
        elif aligner == 'salmon':
            port = Salmon
        else:
            raise Exception('{:>10} : aligner not supported {}'.format(args['aligner']))

        return port


    def run(self):
        port = self.pick_ports()
        port(**self.__dict__).run()
        return (self.align_bam, self.unmap1, self.unmap2)
        

############################################################
## aligner                                                ##
## base-level: aligner
## align, report in json
## 
## input:
##   - fq1: str (required)
##   - fq2: str or None
##   - smp_name: str or None (auto)
##   - index_list: str (required)
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
class AlignConfig2(object):
    """
    Config for single alignment, 1 fq, 1 index
  
    [1, 1] : 1 fastq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs, force=True)
        self.init_args()


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
        """
        default arguments
        required:
          - fq1
          - genome/ext-index
          - outdir
        
        optional:
          - aligner

        prepare index
          - align-to-rRNA
          - align-to-chrM
        """
        args_default = {
            'fq1': None,
            'fq2': None,
            'index_list': None,
            'outdir': str(pathlib.Path.cwd()),
            'index_name': None,
            'smp_name': None,
            'extra_para': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'pickle': None,
            'n_map': 1,
            'unique_only': False,
            'genomeLoad': 'NoSharedMemory',
            }
        self.update(args_default, force=False) # update missing attrs
        
        # 1st level: pickle / update all config
        # fresh start by pickle
        if not self.pickle is None:
            if os.path.exists(self.pickle) and self.pickle.endswith('.pickle'):
                args_pickle = pickle_to_dict(self.pickle)
                args_pickle.pop('pickle', None) # empty, for further round
                self.update(args_pickle, force=True) # fresh start

        ## outdir
        self.outdir = file_abspath(self.outdir)

        ## update
        self.init_fq()
        self.init_index()
        self.init_cpu()
        self.init_files()
        self.is_paired = self.check_paired() # is_paired ?


    def init_fq(self):
        """
        Check out fastq files
        make sure: 
        fq1: str 
        fq2: str or None
        """
        ####################################################
        ## 1. for fastq1 files                             #
        if isinstance(self.fq1, list):
            if len(self.fq1) > 1:
                log.warning('str expected, str found, pick the 1 fastq: {}'.format(
                    next(iter(self.fq1), None)))
            self.fq1 = file_abspath(next(iter(self.fq1), None))

        if isinstance(self.fq1, str):
            self.fq1 = file_abspath(self.fq1)
        else:
            raise Exception('failed fq1: str expected, {} found'.format(
                type(self.fq1).__name__))

        ####################################################
        ## 2. for fastq2 files                             #
        if isinstance(self.fq2, list):
            if len(self.fq2) > 1:
                log.warning('str expected, list found, pick the 1 fastq: {}'.format(
                    next(iter(self.fq2), None)))
            self.fq2 = file_abspath(next(iter(self.fq2), None))

        if self.fq2 is None:
            pass
        elif isinstance(self.fq2, str):
            self.fq2 = file_abspath(self.fq2)
        else:
            raise Exception('failed, fq2, None, str expected: {} found'.format(
                type(self.fq2).__name__))

        ####################################################
        ## 3. for smp_name                                 #
        ## smp_name [from fq1, smp_path]
        if isinstance(self.smp_name, list):
            if len(self.smp_name) > 1:
                log.warning('str expected, list found, pick the 1 smp_name: {}'.format(
                    next(iter(self.smp_name), None)))
            self.smp_name = next(iter(self.smp_name), None)

        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True)
        elif isinstance(self.smp_name, str):
            pass
        else:
            raise Exception('failed, smp_name, str expected, {} found'.format(
                type(self.smp_name).__name__))

        ## fq type (fa or fq)
        self.fq_format = self.fq_type = Fastx(self.fq1).format # to-do !!!

        ## file exists
        chk0 = isinstance(self.fq1, str)
        chk1 = isinstance(self.fq2, str) or self.fq2 is None
        chk2 = isinstance(self.smp_name, str)
        chk3 = file_exists(self.fq1, isfile=True)
        chk4 = file_exists(self.fq2, isfile=True) or self.fq2 is None
        log_msg0 = '\n'.format([
            'Failed on arguments:'
            '{:>30s}: {}'.format('fq1, str or list', self.fq1),
            '{:>30s}: {}'.format('fq2, str, list, None', self.fq2),
            '{:>30s}: {}'.format('smp_name, str or list', self.smp_name)])
        if not all([chk0, chk1, chk2, chk3, chk4]):
            # raise Exception('failed, fq1, fq2: \nfq1: {}, \nfq2: {}'.format(
            #     self.fq1, self.fq2))
            raise Exception(log_msg0)


    def init_index(self):
        """
        Checkout index_list, index_name
        index_list: str
        """
        ## index_list, genome, 
        if isinstance(self.index_list, list):
            if len(self.index_list) > 1:
                log.warning('str expected, str found, pick the 1st index: {}'.format(
                    next(iter(self.index_list), None)))
            self.index_list = next(iter(self.index_list), None)

        if isinstance(self.index_list, str):
            chk0 = AlignIndex(aligner=self.aligner, index=self.index_list).is_index()
            if not chk0:
                raise Exception('index_list, not a actual index: {} {}'.format(
                    self.aligner, self.index_list))
        else:
            raise Exception('failed index_list, str expected, {} found'.format(
                type(self.index_list).__name__))

        ## index_name: 1.rRNA, 2.genome, ...
        if isinstance(self.index_name, list):
            if len(self.index_name) > 1:
                log.warning('str expected, str found, pick the 1st index_name: {}'.format(
                    next(iter(self.index_name), None)))
            self.index_name = next(iter(self.index_name), None)

        if self.index_name is None:
            self.index_name = AlignIndex(aligner=self.aligner, index=self.index_list).index_name()
        elif isinstance(self.index_name, str):
            pass
        else:
            raise Exception('failed index_name, str expected, {} found'.format(
                type(self.index_name).__name__))


    def init_cpu(self):
        """
        threads, CPUs
        """
        ## check number of threads, parallel_jobs
        ## parallel jobs * threads
        n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

        max_jobs = int(n_cpu / 4.0)
        # check parallel_jobs (max: 1/4 of n_cpus)
        if self.parallel_jobs > max_jobs: 
            log.warning('Too large, change --parallel-jobs from {} to {}'.format(
                self.parallel_jobs, max_jobs))
            self.parallel_jobs = max_jobs

        # check threads
        max_threads = int(0.8 * n_cpu / self.parallel_jobs)
        if self.threads * self.parallel_jobs > 0.8 * n_cpu:
            log.warning('Too large, change --threads from {} to {}'.format(
                self.threads, max_threads))
            self.threads = max_threads        


    def init_files(self, create_dirs=True):
        """
        default files, paths 
        """
        create_dirs = True
        #############################################################
        ## save file_path / default
        ## update these attributes, everytime
        ## auto-generated
        project_dir = os.path.join(self.outdir, self.smp_name, self.index_name)
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'align_cmd_file': os.path.join(project_dir, 'cmd.sh'),
            'report_dir': os.path.join(project_dir, 'report'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json')
        }
        self.update(auto_files, force=True) # key point
        ## create directories
        if create_dirs is True:
            check_path([
                self.project_dir,
                self.config_dir,
                self.report_dir])

        smp_name = self.smp_name # 
        #############################################################
        ## default files
        self.align_log = os.path.join(project_dir, smp_name + '.log')
        self.align_json = os.path.join(project_dir, smp_name + '.json')
        self.align_stat = os.path.join(project_dir, smp_name + '.stat')
        self.align_sam = os.path.join(project_dir, smp_name + '.sam')
        self.align_bam = os.path.join(project_dir, smp_name + '.bam')
        self.unmap_prefix = os.path.join(project_dir, smp_name + '.unmap.' + self.fq_format)
        if self.is_paired is True:
            self.unmap1 = os.path.join(project_dir, smp_name + '.unmap.1.' + self.fq_format)
            self.unmap2 = os.path.join(project_dir, smp_name + '.unmap.2.' + self.fq_format)
        else:
            self.unmap1 = os.path.join(project_dir, smp_name + '.unmap.' + self.fq_format)
            self.unmap2 = None


    def check_paired(self):
        """
        Determine the type of Alignment:
        is_paired: bool
        is_multi_index: bool ?
        """
        is_paired = False
        if self.fq2 is None:
            is_paired = False
        else:
            # check file names
            if fq_name(self.fq1, pe_fix=True) == fq_name(self.fq2, pe_fix=True):
                is_paired = True

        return is_paired


class Bowtie(object):
    """
    Run bowtie for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.aligner = 'bowtie' # force changed.
        self.init_args()


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
        """
        check
        """
        local_config = AlignConfig2(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes

        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## check index (aligner)
        chk2 = AlignIndex(aligner=self.aligner, index=self.index_list).is_index()
        if not chk2:
            raise Exception('not a Bowtie index: {}'.format(self.index_list))

        

    def get_cmd(self):
        """
        unique
        n_map
        extra_para
        """
        aligner_exe = shutil.which('bowtie')
        
        # extra_para
        if self.extra_para is None:
            self.extra_para = ''

        # nmap
        if self.n_map < 1: 
            self.n_map = 1 # 
 
        # unique
        if self.unique_only is True: 
            c_unique = '-m 1'
        else:
            c_unique = '-v 2 -k {}'.format(self.n_map)

        # fq type
        c_fx = '-f' if self.fq_format == 'fasta' else '-q'

        # common
        c = '{} {} '.format(aligner_exe, self.index_list)
        c += '{} {} {} -p {} '.format(c_unique, c_fx, self.extra_para, self.threads)
        c += '--mm --best --sam --no-unal --un {} '.format(self.unmap_prefix)

        # SE or PE
        if self.is_paired:
            c += '-1 {} -2 {} 1> {} 2> {}'.format(
                self.fq1, 
                self.fq2, 
                self.align_sam, 
                self.align_log)
        else:
            c += '{} 1> {} 2> {}'.format(
            self.fq1, 
            self.align_sam, 
            self.align_log)

        # convert sam to bam
        c_samtools_para = ''
        c += '&& samtools view -Sub -F 0x4 {} {} | samtools sort -o {} -'.format(
            c_samtools_para,
            self.align_sam,
            self.align_bam)

        # rename the unmap files for PE reads
        # bowtie auto generated: _1.fastq, _2.fastq
        # expect output:         .1.fastq, .2.fastq
        if self.is_paired:
            unmap1 = fq_name(self.unmap_prefix, include_path=True, pe_fix=False) + '_1.fastq'
            unmap2 = fq_name(self.unmap_prefix, include_path=True, pe_fix=False) + '_2.fastq'
            c += '&& mv {} {} '.format(unmap1, self.unmap1)
            c += '&& mv {} {} '.format(unmap2, self.unmap2)

        return c


    def read_log(self, to_json=True):
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
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')



        ## save to json
        if to_json:
            Json(dd).writer(self.align_json)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        cmd = self.get_cmd()

        if file_exists(self.align_bam):
            log.warning('file exists, alignment skipped: {} {}'.format(self.smp_name, self.index_name))
        else:
            with open(self.align_cmd_file, 'wt') as w:
                w.write(cmd + '\n')
            run_shell_cmd(cmd)
            self.read_log(to_json=True)

        ## chek unmap file
        if self.is_paired:
            unmap1, unmap2 = (self.unmap1, self.unmap2)
        else:
            unmap1, unmap2 = (self.unmap_prefix, None)
        chk0 = os.path.exists(self.align_bam)
        chk1 = os.path.exists(unmap1)
        chk2 = os.path.exists(unmap2) if self.is_paired else True
        if not all([chk0, chk1, chk2]):
            raise Exception('Check the output files: {}'.format(self.project_dir))

        return (self.align_bam, unmap1, unmap2)


class Bowtie2(object):
    """
    Run bowtie2 for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.aligner = 'bowtie2' # force changed.
        self.init_args()


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
        """
        check
        """
        local_config = AlignConfig2(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes

        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## check index (aligner)
        chk2 = AlignIndex(aligner=self.aligner, index=self.index_list).is_index()
        if not chk2:
            raise Exception('not a Bowtie2 index: {}'.format(self.index_list))


    def get_cmd(self):
        """
        unique
        n_map
        extra_para
        """
        aligner_exe = shutil.which('bowtie2')
 
        # nmap
        if self.n_map < 1: 
            self.n_map = 1 # 

        # extra_para
        if self.extra_para is None:
            self.extra_para = ''

        # nmap
        c_multi = '-k {}'.format(self.n_map) if self.n_map > 0 else ''
 
        # fq type
        c_fx = '-f' if self.fq_format == 'fasta' else '-q'

        # common
        c = '{} -x {} '.format(aligner_exe, self.index_list)
        c += '{} {} -p {} '.format(c_fx, self.extra_para, self.threads)
        c += '--mm --sensitive-local --no-unal '

        # SE or PE
        if self.is_paired:
            c += '--un-conc {} -1 {} -2 {} 1> {} 2> {}'.format(
                self.unmap_prefix, # unmap.fastq 
                self.fq1, 
                self.fq2, 
                self.align_sam, 
                self.align_log)
        else:
            c += '--un {} -U {} 1> {} 2> {}'.format(
                self.unmap_prefix, # unmap.fastq 
                self.fq1, 
                self.align_sam, 
                self.align_log)

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
                c += ' '.join([
                    '&& samtools view -Sb',
                    '<(samtools view -H {} ;'.format(self.align_sam),
                    "samtools view -F 2048 {} | grep 'YT:Z:CP')".format(self.align_sam),
                    '| samtools sort -o {} -'.format(self.align_bam)])
            else:
                c += ' '.join([
                    '&& samtools view -Sb -F 2048 {}'.format(self.align_sam),
                    '| samtools sort -o {} -'.format(self.align_bam)])
        else:
            c += ' '.join([
                    '&& samtools view -Sb -F 2048 {}'.format(self.align_sam),
                    '| samtools sort -o {} -'.format(self.align_bam)])

        return c


    def read_log(self, to_json=True):
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
                if len(re.sub('[0-9]', '', value)) > 0:
                    continue
                if '%' in value:
                    continue
                if line.strip().startswith('----'):
                    continue
                value = eval(value)

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
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')


        ## save to json
        if to_json:
            Json(dd).writer(self.align_json)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        cmd = self.get_cmd()

        if file_exists(self.align_bam):
            log.warning('file exists, alignment skipped: {} {}'.format(self.smp_name, self.index_name))
        else:
            try:
                with open(self.align_cmd_file, 'wt') as w:
                    w.write(cmd + '\n')
                run_shell_cmd(cmd)
                self.read_log(to_json=True)
            except:
                log.error('Bowtie2() failed, outdir: {}'.format(
                    self.outdir))
        ## chek unmap file
        if self.is_paired:
            unmap1, unmap2 = (self.unmap1, self.unmap2)
        else:
            unmap1, unmap2 = (self.unmap_prefix, None)
        chk0 = os.path.exists(self.align_bam)
        chk1 = os.path.exists(unmap1)
        chk2 = os.path.exists(unmap2) if self.is_paired else True
        if not all([chk0, chk1, chk2]):
            raise Exception('Check the output files: {}'.format(self.project_dir))

        return (self.align_bam, unmap1, unmap2)


## caution:
## --genomeLoad: LoadAndRemove (NoSharedMemory)
class STAR(object):
    """
    Run bowtie2 for: 1 fq, 1 index
    [fq1|fq2], [index_list], [smp_name], [index_name]
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.aligner = 'STAR' # force changed.
        self.init_args()


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
        """
        check
        """
        local_config = AlignConfig2(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        
        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk1 = args_logger(self.__dict__, self.config_txt)

        ## prefix for files
        ## STAR auto generate files
        self.align_prefix = os.path.join(self.project_dir, self.smp_name)

        ## check genome size (chrLength.txt)
        chrSize = 0
        chrLength = os.path.join(self.index_list, 'chrLength.txt')
        if os.path.exists(chrLength):
            with open(chrLength) as r:
                for line in r:
                    chrSize += int(line.strip())
        self.small_genome = True if chrSize < 1000000 else False # 1M genome

        ## genomeLoad
        gl = ['NoSharedMemory', 'LoadAndKeep', 'LoadAndRemove', 'LoadAndExit', 'Remove', 'NoSharedMemory']
        if not self.genomeLoad in gl:
            log_msg0 = '\n'.join([
                'unknown --genomeLoad: {}'.format(self.genomeLoad),
                'expected: {}'.format(' '.join(gl)),
                'auto switch to: NoSharedMemory'
                ])
            log.warning(log_msg0)
            self.genomeLoad = 'NoSharedMemory'


    def get_cmd(self):
        """
        unique
        n_map
        extra_para
        """
        aligner_exe = shutil.which('STAR')
        
        ## for small genome
        ## mapping takes too long,
        ## 99% of the reads not mapped
        ## change --seedPerWindowNmax
        ## max number of seeds per window
        ## https://github.com/alexdobin/STAR/issues/329
        ## https://groups.google.com/d/msg/rna-star/hJL_DUtliCY/HtpiePlMBtYJ
        if self.small_genome:
            log.warning('STAR on small genome (<1 Mb): {}'.format(self.index_list))
            seed_max = 5 # even smaller
        else:
            seed_max = 50 # default

        ## for unique map
        ## --outFilterMultimapNmax
        ## maximum number of loci the read is allowed to map to. [default: 10]
        ##
        ## or filt unique reads by: samtools view -q 255
        # n_map = 1 if self.unique_only else n_map
        if self.unique_only: self.n_map = 1
        c_unique = ' '.join([
            '--outFilterMultimapNmax {}'.format(self.n_map),
            '--seedPerWindowNmax {}'.format(seed_max)])

        ## STAR, .gz file input
        c_reader = 'zcat' if self.fq1.endswith('.gz') else '-'

        ## fq2
        if self.fq2 is None: self.fq2 = '' # empty

        ## extra para
        c_extra_para = '' if self.extra_para is None else self.extra_para

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
        c = ' '.join([
            aligner_exe,
            '--genomeLoad {}'.format(self.genomeLoad), # default: NoSharedMemory
            '--runMode alignReads',
            '--genomeDir {}'.format(self.index_list),
            '--readFilesIn {} {}'.format(self.fq1, self.fq2),
            '--readFilesCommand {}'.format(c_reader),
            '--outFileNamePrefix {}'.format(self.align_prefix),
            '--runThreadN {}'.format(self.threads),
            '--limitBAMsortRAM 10000000000',
            '--outSAMtype BAM SortedByCoordinate',
            '--outFilterMismatchNoverLmax 0.07',
            '--seedSearchStartLmax 20',
            '--outReadsUnmapped Fastx', # self.unmap1,
            c_unique,
            c_extra_para])

        return c


    def update_names(self, keep_old=False):
        """
        Update the filenames of STAR output
        bam: *Aligned.sortedByCoord.out.bam -> *.bam # mv
        log: *Log.final.out -> *.log # copy
        log: *Log.out -> *.out
        unmap: *Unmapped.out.mate1 -> *.unmap.1.fastq
               *Unmapped.out.mate1 -> *.unmap.1.fastq
        """
        ## default output of STAR
        bam_from = self.align_prefix + 'Aligned.sortedByCoord.out.bam'
        log_from = self.align_prefix + 'Log.final.out'
        unmap1 = self.align_prefix + 'Unmapped.out.mate1'
        unmap2 = self.align_prefix + 'Unmapped.out.mate2'

        if keep_old:
            return bam_from
        else:
            ## rename files
            ## move BAM, unmap
            ## copy log
            flag = 0
            if file_exists(bam_from) and not file_exists(self.align_bam):
                shutil.move(bam_from, self.align_bam)
                flag += 1
            if file_exists(log_from):
                shutil.copy(log_from, self.align_log)
                flag += 1
            ## for unmap files
            if file_exists(unmap1):
                shutil.move(unmap1, self.unmap1)
                flag += 1
            if file_exists(unmap2):
                shutil.move(unmap2, self.unmap2)
                flag += 1


    def read_log(self, to_json=True):
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
            groups = ['total', 'map', 'unique', 'multiple', 'unmap', 'fqname', 'index_name']
            h = '\t'.join(groups)
            v = '\t'.join([str(dd.get(i, 0)) for i in groups])
            w.write('#' + h + '\n')
            w.write(v + '\n')

        ## save to json
        if to_json:
            Json(dd).writer(self.align_json)

        return dd['total'], dd['map'], dd['unique'], dd['multiple'], dd['unmap']


    def run(self):
        cmd = self.get_cmd()

        if file_exists(self.align_bam):
            log.warning('file exists, alignment skipped: {}'.format(self.align_bam))
        else:
            try:
                with open(self.align_cmd_file, 'wt') as w:
                    w.write(cmd + '\n')
                run_shell_cmd(cmd)
                self.update_names() # update
                self.read_log(to_json=True) # save to json, stat
            except:
                log.error('STAR().run() failed, outdir: {}'.format(
                    self.project_dir))

        ## chek unmap file
        if self.is_paired:
            unmap1, unmap2 = (self.unmap1, self.unmap2)
        else:
            unmap1, unmap2 = (self.unmap_prefix, None)
        chk0 = os.path.exists(self.align_bam)
        chk1 = os.path.exists(unmap1)
        chk2 = os.path.exists(unmap2) if self.is_paired else True
        if not all([chk0, chk1, chk2]):
            raise Exception('Check the output files: {}'.format(self.project_dir))

        return (self.align_bam, unmap1, unmap2)


## to-do
class BWA(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)


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
        pass
        

    def run(self):
        pass


class Hisat2(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)


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
        pass
        

    def run(self):
        pass


class Kallisto(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)


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
        pass
        

    def run(self):
        pass


class Salmon(object):
    """
    Run bowtie for: 1 fq, 1 index
    """
    def __init__(self, **kwargs):
        self.update(kwargs)


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

