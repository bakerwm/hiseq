#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Working mode:

1. build design

Create RNAseq_auto_design.txt for further analysis

2. pickle

Pass all arguments from *.pickle file

3. design

pass values from design.txt



## Design

input:

- control (rep1, rep2, ...)
- treatment (rep1, rep2, ...)
- genome

output:

  - single:
    - raw_data
    - clean_data
    - align
    - count
    - qc
    - report
    ...

  - pair:
    - count
    - deseq
    - qc
    - report


## figures/tables
1. scatter plot
2. MA plot
3. volcano
4. heatplot
...

## subfunctions:

class RNAseqSingle()
class RNAseqPair()
class RNAseqConfig()
...
object to dict, ...


DE analysis:

DEseq2, edgeR, ...



Alignment:

gene,
TE
piRNA cluster
...

2020-04-18
## 1. save arguments in self.
## 2. save default paths, files, in self.files
## 3. update all args/

"""

import os
import re
import hiseq
import pysam
import shutil
import tempfile
import pandas as pd
from multiprocessing import Pool
import copy # copy objects
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trimmer
from hiseq.align.alignment import Alignment
from hiseq.rnaseq.rnaseq_builder import RNAseqFqDesign
from hiseq.atac.atac_utils import Bam2bw


def print_df(x):
    # for k, v in x.items():
    #     print('{:30s} : {}'.format(k, v))
    for k in sorted(x.keys()):
        print('{:30s} : {}'.format(k, x[k]))


class RNAseq(object):
    def __init__(self, **kwargs):
        """
        options:
        global options: feature(gene/te/piRNAcluster/...), gtf,
        1. design: design.txt (single/multiple/deseq)
        2. single/multiple: --fq1, --fq2, --genome, --outdir
        3. deseq: --dirs-ctl, --dirs-exp, --outdir
        ...
        """
        self.update(kwargs) # fresh new
        local_config = RNAseqConfig(**self.__dict__)
        self.update(local_config.__dict__, force=True)
        self.pickle = None # terminate `pickle` option

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
        Run all RNAseq analysis
        """
        args = self.__dict__.copy()



        if self.rnaseq_type == 'build_design':
            # RNAseqBuildDesign(**args).run()
            print('!AAAA1')
            RNAseqFqDesign(**args)
        elif self.rnaseq_type == 'deseq_single':
            print('!AAAA2')
            RNAseqDeseqSingle(**args).run()
        elif self.rnaseq_type == 'deseq_multiple':
            print('!AAAA3')
            RNAseqDeseqMultiple(**args).run()
        elif self.rnaseq_type == 'rnaseq_single':
            print('!AAAA4')
            RNAseqSingle(**args).run()
        elif self.rnaseq_type == 'rnaseq_multiple':
            print('!AAAA5')
            # env: rnaseq_multiple
            smp_path = RNAseqMultiple(**args).run()
            # env: deseq_multiple
            args_i = self.__dict__.copy()
            args_i['rnaseq_type'] = 'rnaseq_single' # update
            args_i['smp_path'] = smp_path
            local_config = RNAseqConfig(**args_i)
            args_i.update(local_config.__dict__, force=True)
            args_i['smp_path'] = smp_path # covered by RNAseqConfig
            RNAseqDeseqMultiple(**args_i).run()
        else:
            pass


class RNAseqConfig(object):
    """
    UPDATED:
    directory type: single, pair
    config: init_args, single, merge, multiple
    list all files (objects)

    required args:
    fq1, fq2, genome, outdir, ...

    options:
    """
    def __init__(self, **kwargs):
        self.update(kwargs) # update from args
        self.init_args() # prepare args


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
        check arguments conflicts, defaults...
        """
        # default values
        ## critical args:
        ## PE mode get low pct unique mapped reads, but SE mode not.
        ## so force to SE mode
        args_default = {
            'align_to_rRNA': True, # default
            'aligner': 'STAR', # default
            'binsize': 10,
            'build_design': False, # boolean
            'ctl': None,
            'ctl_read2': None,
            'design': None, # json file
            'dirs_ctl': None,
            'dirs_exp': None,
            'exp': None,
            'exp_read2': None,
            'feature': 'gene',
            'fq1': None, 
            'fq2': None,
            'genome': 'dm6',
            'genomeLoad': 'LoadAndRemove',
            'genome_size': None,
            'group': None,
            'gtf': None,
            'outdir': str(pathlib.Path.cwd()),
            'overwrite': False,
            'parallel_jobs': 1,
            'read1_only': False,
            'smp_name': None,
            'smp_path': None,
            'threads': 1,
            'unique_only': True,
        }
        self.update(args_default, force=False) # update missing attrs
        # 1st level: pickle / update all config
        # # fresh start by pickle
        # if not self.pickle is None:
        #     if os.path.exists(self.pickle) and self.pickle.endswith('.pickle'):
        #         args_pickle = pickle_to_dict(self.pickle)
        #         args_pickle.pop('pickle', None) # empty, next round
        #         self.update(args_pickle, force=True) # fresh start

        # 2nd level: design (input)
        # fresh start by design (specific args)
        # if not self.design is None:
        #     # self.init_design_arg(create_dirs)
        #     args_design = DesignReader(self.design).to_dict()
        #     self.update(args_design, force=True) # update specific args

        ## 1st level: creat design.json
        if self.build_design:
            # require: ctl, exp
            if self.ctl is None or self.exp is None:
                log.warning('--ctl, --exp required for build-design')
        else:
            self.init_fq()
            self.init_group()

        ## outdir
        self.outdir = file_abspath(self.outdir)

        ## gtf
        if isinstance(self.gtf, list):
            self.gtf = next(iter(self.gtf), None)
            log.warning('gtf, list found, choose the first one: {}'.format(self.gtf))

        if self.gtf is None:
            self.gtf = Genome(genome=self.genome).gene_gtf('ensembl')
            log.warning('gtf, auto-geenrated: {}'.format(self.gtf))
        elif isinstance(self.gtf, str):
            if not os.path.exists(self.gtf):
                raise Exception('gtf, file not exists: {}'.format(self.gtf))
        else:
            raise Exception('failed gtf, str expected, {} found'.format(type(self.gtf).__name__))

        ## sub-cmd
        # self.init_fq()
        # self.init_group()
        self.init_cpu()
        self.init_mission()


    def init_fq(self):
        """
        check out input fastq files
        input: could be None, str, list
        output: None or list
        """
        ####################################################
        # 1st level: read1 only                            #
        if self.read1_only is True:
            self.fq2 = None

        ####################################################
        ## 1. for fastq1 files                             #
        if self.fq1 is None:
            pass
        elif isinstance(self.fq1, str):
            self.fq1 = file_abspath(self.fq1)
        elif isinstance(self.fq1, list):
            self.fq1 = file_abspath(self.fq1)
        else:            
            log.error('failed fq1, None, str, list expected, {} found'.format(type(self.fq1).__name__))

        ####################################################
        ## 2. for fastq2 files                              #
        if self.fq2 is None:
            pass
        elif isinstance(self.fq2, str):
            self.fq2 = file_abspath(self.fq2)
        elif isinstance(self.fq2, list):
            self.fq2 = file_abspath(self.fq2)
        else:
            log.error('failed fq2, None, str, list expected, {} found'.format(type(self.fq2).__name__))

        ####################################################
        ## 3. for smp_name                                 #
        ## smp_name [from fq1, smp_path]
        if self.smp_name is None:
            if isinstance(self.fq1, list):
                self.smp_name = fq_name(self.fq1, pe_fix=True)
            elif isinstance(self.smp_path, list):
                self.smp_name = fq_name(self.smp_path, pe_fix=True)
            else:
                log.warning('errors for smp_name')
                pass
        elif isinstance(self.smp_name, str):
            self.smp_name = self.smp_name
        elif isinstance(self.smp_name, list):
            pass
        else:
            log.error('failed, smp_name, None, str, list expected, {} found'.format(type(self.smp_name).__name__))


    def init_group(self):
        """
        Determine the group names of samples
        input:
        output: list or None
        """
        if self.group is None:
            if isinstance(self.smp_path, list):
                self.group = fq_name_rmrep(self.smp_path)
            elif isinstance(self.smp_name, list):
                self.group = fq_name_rmrep(self.smp_name)
            else:
                log.warning('group is None')
                pass # None
        elif isinstance(self.group, str):
            self.group = [self.group]
        elif isinstance(self.group, list):
            pass
        else:
            log.error('failed, group, None, list expected, {} found'.format(type(self.group).__name__))


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


    def init_mission(self):
        """
        Determine the purpose the RNAseq analysis
        1. build_design
        2. deseq_multiple
        3. deseq_single
        4. rnaseq_multiple
        5. rnaseq_single
        """
        create_dirs = getattr(self, 'create_dirs', True)

        # create design?!
        if self.build_design is True:
            self.rnaseq_type = 'build_design' # del 
            self.init_build_design() # fresh start 
        # read data from design.json
        elif file_exists(self.design):  
            self.rnaseq_type = 'deseq_multiple'
            self.init_deseq_multiple(create_dirs)
        # # determine the rnaseq type:
        # elif isinstance(self.smp_path, list):
        #     self.rnaseq_type = 'deseq_multiple'
        #     self.init_deseq_multiple(create_dirs)
        elif isinstance(self.dirs_ctl, list) and isinstance(self.dirs_exp, list):
            self.rnaseq_type = 'deseq_single'
            self.init_deseq_single(create_dirs)
        elif isinstance(self.fq1, list) and len(self.fq1) >= 1:
            self.rnaseq_type = 'rnaseq_multiple'
            self.init_rnaseq_multiple(create_dirs)
        elif isinstance(self.fq1, str):
            self.rnaseq_type = 'rnaseq_single'
            self.init_rnaseq_single(create_dirs)
        else:
            self.rnaseq_type = None # unknown
            log_msg0 = '\n'.join([
                'Failed to determine the type of RNAseq',
                '{:>30s} : {}'.format('RNAseq single', '[fq1], [fq2], genome'),
                '{:>30s} : {}'.format('RNAseq multiple', '[fq1], [fq2], genome'),
                '{:>30s} : {}'.format('DESeq single', '[dirs_ctl], [dirs_exp]'),
                '{:>30s} : {}'.format('DESeq multiple', '[smp_path]')])
            raise Exception(log_msg0)


    def init_build_design(self):
        """
        Create design for RNAseq
        ctl,
        exp, 
        (fq1, fq2)
        
        arguments:
        --ctl
        --exp
        --ctl-read2
        --exp-read2
        --output
        """
        if self.ctl is None or self.exp is None:
            raise Exception('--ctl, --exp required for build-design')
        pass


    def init_deseq_multiple(self, create_dirs=True):
        """
        Read ctl/exp from design.json
        remove the "design", after parsing
        """
        pass


    def init_deseq_single(self, create_dirs=True):
        """
        prepare args, files for DEseq single
        input: dirs_ctl, dirs_exp, genome
        output: deseq_single
        """
        ## smp_name
        smp_path = self.dirs_ctl + self.dirs_exp
        smp_path = [i.rstrip('/') for i in smp_path]
        if self.smp_name is None:
            self.smp_name = fq_name(smp_path)
        if self.group is None:
            self.group = fq_name_rmrep(smp_path)

        ## check
        chk0 = isinstance(self.dirs_ctl, list) and len(self.dirs_ctl) > 0
        chk1 = isinstance(self.dirs_exp, list) and len(self.dirs_exp) > 0
        chk2 = isinstance(self.genome, str)
        chk3 = all(RNAseqReader(i, feature=self.feature).is_rnaseq_single() for i in self.dirs_ctl + self.dirs_exp)
        chk4 = len(self.smp_name) == len(self.dirs_ctl + self.dirs_exp)
        if not all([chk0, chk1, chk2, chk3, chk4]):
            raise Exception('fq1, fq2, genome, failed')

        # absolute path
        self.dirs_ctl = file_abspath(self.dirs_ctl)
        self.dirs_exp = file_abspath(self.dirs_exp)
        self.prefix_ctl, self.prefix_exp = list_uniquer(self.group, sorted=False)
        self.project_name = '{}.vs.{}'.format(self.prefix_ctl, self.prefix_exp)
        self.project_dir = os.path.join(self.outdir, self.project_name, self.feature)
        self.config_dir = os.path.join(self.project_dir, 'config')
        #############################################################
        ## save file_path / default  !, not required
        ## update these attributes, everytime
        ## auto-generated
        auto_files = {
            'project_dir': self.project_dir,
            'config_dir': self.config_dir,
            'count_dir': os.path.join(self.project_dir, 'count'),
            'deseqdir': os.path.join(self.project_dir, 'deseq'),
            'enrichdir': os.path.join(self.project_dir, 'enrich'),
            'pdfdir': os.path.join(self.project_dir, 'pdf'),
            'report_dir': os.path.join(self.project_dir, 'report'),
            'config_txt': os.path.join(self.config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(self.config_dir, 'arguments.pickle'),
            'config_json': os.path.join(self.config_dir, 'arguments.json'),
            'deseq_design': os.path.join(self.config_dir, 'deseq_design.txt'),
            'count_ctl': [RNAseqReader(i, feature=self.feature).count_sens for i in self.dirs_ctl],
            'count_exp': [RNAseqReader(i, feature=self.feature).count_sens for i in self.dirs_exp],
        }
        self.update(auto_files, force=True) # !!!! key-point

        ## create directories
        if create_dirs is True:
            check_path([
                self.config_dir,
                self.count_dir,
                self.deseqdir,
                self.enrichdir,
                self.report_dir])


    def init_rnaseq_multiple(self, create_dirs=True):
        """
        prepare args, files for multiple
        input: fq1, fq2, genome, outdir
        output: rnaseq_single
        """
        ## smp_name
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True)
        ## check
        chk0 = isinstance(self.fq1, list) and len(self.fq1) > 0
        chk1 = isinstance(self.fq2, list) and len(self.fq2) > 0 or self.fq2 is None
        chk2 = all(file_exists(self.fq1))
        chk3 = file_exists(self.fq2) or self.fq2 is None
        chk4 = isinstance(self.genome, str)
        chk5 = isinstance(self.smp_name, list) and len(self.smp_name) > 0
        log_msg0 = '\n'.join([
            'failed for RNAseq multiple',
            '{:>30s} : {} {}'.format('fq1, list > 0', chk0, self.fq1),
            '{:>30s} : {} {}'.format('fq1 exists', chk2, self.fq1),
            '{:>30s} : {} {}'.format('fq2, None or list > 0', chk1, self.fq2),
            '{:>30s} : {} {}'.format('fq2 None or list', chk3, self.fq2),
            '{:>30s} : {} {}'.format('genome, str', chk4, self.genome),
            '{:>30s} : {} {}'.format('smp_name, list > 0', chk5, self.smp_name),
            ])
        if not all([chk0, chk1, chk2, chk3, chk4, chk5]):
            raise Exception(log_msg0)


        #############################################################
        ## save file_path / default  !, not required
        ## update these attributes, everytime
        ## auto-generated
        project_dir = os.path.join(self.outdir, 'config', self.feature)
        config_dir = os.path.join(project_dir, 'config')
        auto_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'report_dir': os.path.join(project_dir, 'report'),
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json'),
        }
        self.update(auto_files, force=True) # !!!! key-point
        #############################################################
 
        # create
        if create_dirs is True:
            check_path([
                self.config_dir, 
                self.report_dir])


    def init_rnaseq_single(self, create_dirs=True):
        """
        prepare args, files for single
        input: fq1, fq2, genome, outdir : str
        output: rnaseq_single
        """
        ## smp_name
        if self.smp_name is None:
            self.smp_name = fq_name(self.fq1, pe_fix=True)        
        ## check
        chk0 = isinstance(self.fq1, str)
        chk1 = isinstance(self.fq2, str) or self.fq2 is None
        chk2 = file_exists(self.fq1)
        chk3 = file_exists(self.fq2) or self.fq2 is None
        chk4 = isinstance(self.genome, str)
        chk5 = isinstance(self.smp_name, str)
        log_msg1 = '\n'.join([
            'arguments failed',
            '{:>30s}: {} {}'.format('fq1, str', chk0, self.fq1),
            '{:>30s}: {} {}'.format('fq1, exists', chk2, self.fq1),
            '{:>30s}: {} {}'.format('fq2, str or None', chk1, self.fq2),
            '{:>30s}: {} {}'.format('fq1, exists or None', chk3, self.fq2),
            '{:>30s}: {} {}'.format('genome, str', chk4, self.genome),
            '{:>30s}: {} {}'.format('smp_name, str', chk5, self.smp_name),
            ])
        if not all([chk0, chk1, chk2, chk3, chk4, chk5]):
            #raise Exception('fq1, fq2, genome, failed')
            raise Exception(log_msg1)

        #############################################################
        ## save file_path / default  !, not required
        ## update these attributes, everytime
        ## auto-generated
        # smp_name = next(iter(self.smp_name), None) # list -> str
        project_dir = os.path.join(self.outdir, self.feature)
        config_dir = os.path.join(project_dir, 'config')
        init_files = {
            'project_dir': project_dir,
            'config_dir': config_dir,
            'config_txt': os.path.join(config_dir, 'arguments.txt'),
            'config_pickle': os.path.join(config_dir, 'arguments.pickle'),
            'config_json': os.path.join(config_dir, 'arguments.json'),
            'raw_dir': os.path.join(project_dir, 'raw_data'),
            'clean_dir': os.path.join(project_dir, 'clean_data'),
            'align_dir': os.path.join(project_dir, 'align'),
            'bam_dir': os.path.join(project_dir, 'bam_files'),
            'bw_dir': os.path.join(project_dir, 'bw_files'),
            'count_dir': os.path.join(project_dir, 'count'),
            'report_dir': os.path.join(project_dir, 'report'),
            'out_prefix': os.path.join(project_dir, self.smp_name)
        }
        self.update(init_files, force=True) # !!!! key-point
        #############################################################
        ## local files: list -> str
        # fq1 = next(iter(self.fq1), None)
        # fq2 = None if self.fq2 is None else next(iter(self.fq2), None)
        # smp_name = next(iter(self.smp_name), None)
        ## raw data
        self.raw_fq_list = [os.path.join(self.raw_dir, os.path.basename(self.fq1))]
        fq2_raw = None if self.fq2 is None else os.path.join(self.raw_dir, os.path.basename(self.fq2))
        self.raw_fq_list.append(fq2_raw)

        ## clean data
        self.clean_fq_list = [os.path.join(self.clean_dir, fq_name(self.fq1, pe_fix=False) + '.fq.gz')]
        fq2_clean = None if self.fq2 is None else os.path.join(self.clean_dir, fq_name(self.fq2, pe_fix=False) + '.fq.gz')
        self.clean_fq_list.append(fq2_clean)

        ## files
        self.trim_stat = os.path.join(self.clean_dir, self.smp_name + '.qc.stat')
        self.bam_raw = os.path.join(self.align_dir, self.smp_name, '2.*', self.smp_name + '.bam')
        self.bam_out = os.path.join(self.bam_dir, self.smp_name + '.bam')
        self.align_stat = os.path.join(self.align_dir, self.smp_name + '.align.txt')
        self.bw_fwd = os.path.join(self.bw_dir, self.smp_name + '.fwd.bigWig')
        self.bw_rev = os.path.join(self.bw_dir, self.smp_name + '.rev.bigWig')
        self.count_sens = os.path.join(self.count_dir, 'count.sens.txt')
        self.count_anti = os.path.join(self.count_dir, 'count.anti.txt')
        self.strandness_status = os.path.join(self.count_dir, 'strandness_status.out')

        ## create directories
        if create_dirs is True:
            check_path([
                self.config_dir,
                self.raw_dir,
                self.clean_dir,
                self.align_dir,
                self.bam_dir,
                self.bw_dir,
                self.count_dir,
                self.report_dir])


class RNAseqMultiple(object):
    """
    Run RNAseq for multiple samples; RNAseqSingle()
    if groups > 1: run DEseq()    
    """
    def __init__(self, **kwargs):
        self.update(kwargs) # fresh new
        self.init_args() # update all variables: *.config, *.args


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
        local_config = RNAseqConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update unknown args
        assert self.rnaseq_type == 'rnaseq_multiple'

        # check arguments
        chk0 = args_checker(self.__dict__, self.config_pickle)
        chk1 = args_logger(self.__dict__, self.config_txt)


    def report(self):
        """
        Create alignment report for RNAseq multiple
        1. trim
        2. align
        3. quant
        ...
        """        
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'rnaseq_report.R')
        cmd_file = os.path.join(self.report_dir, 'cmd.sh')
        report_html = os.path.join(
            self.report_dir, 
            'rnaseq_report.html')

        cmd = 'Rscript {} {} {} {}'.format(
            qc_reportR,
            self.outdir,
            self.report_dir,
            self.feature)

        # save cmd
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')

        if check_file(report_html):
            log.info('RNAseq multiple report() skipped, file exists: {}'.format(
                report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run_rnaseq_single(self, fq1):
        """
        For parallel reason, 
        wrap analysis on single fastq files, into function
        """
        args_local = self.__dict__.copy()
        i = self.fq1.index(fq1)

        # fq1,fq2,smp_name,...
        fq2_list = args_local.get('fq2', None)
        fq2 = fq2_list[i] if isinstance(fq2_list, list) else None
        # smp_name
        smp_name_list = args_local.get('smp_name', None)
        smp_name = smp_name_list[i] if isinstance(smp_name_list, list) else None
        # group
        group_list = args_local.get('group', None)
        group = group_list[i] if isinstance(group_list, list) else None
        # outdir
        outdir = os.path.join(args_local['outdir'], smp_name)

        # update args
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'smp_name': smp_name,
            'group': group,
            'outdir': outdir,
            'smp_path': None,
            'dirs_ctl': None,
            'dirs_exp': None,
            'pickle': None
        }
        args_local.update(args_init)

        RNAseqSingle(**args_local).run()


    def run(self):
        """
        check each sample
        save to outdir/fq_name
        """
        ## prepare dirs for all single files
        smp_dirs = []
        for i, fq1 in enumerate(self.fq1):
            smp_dirs.append(os.path.join(self.outdir, self.smp_name[i]))

        ## Pool() run in parallel
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.run_rnaseq_single, self.fq1)

        # for a, b in enumerate(self.fq1):
        #     self.run_rnaseq_single(b)

        ## create report
        # self.report()

        return smp_dirs


class RNAseqSingle(object):
    def __init__(self, **kwargs):
        """
        Alignment one fastq file, quantification
        """
        self.update(kwargs) # init self, fresh new
        self.init_args() # update all variables: *.config, *.args


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
        local_config = RNAseqConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes

        # save arguments
        chk1 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk2 = args_logger(self.__dict__, self.config_txt)

        # check        
        chka1 = isinstance(self.fq1, str)
        chka2 = isinstance(self.fq2, str) or self.fq2 is None
        chka3 = isinstance(self.genome, str)
        chka4 = isinstance(self.outdir, str)
        if not all([chka1, chka2, chka3, chka4]):
            raise Exception('value error: \nfq1: {}\ngenome: {}\noutdir: {}'.format(
                self.fq1, self.genome, self.outdir))


    #######################################
    ## main pipeline
    def prep_raw(self, copy=False):
        """
        Copy raw data to dest dir
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list

        # copy
        if copy is True:
            shutil.copy(self.fq1, raw_fq1)
            shutil.copy(self.fq2, raw_fq2)
        else:
            symlink(self.fq1, raw_fq1, absolute_path=True)
            symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self, trimmed=False):
        """
        Trim reads:
        hiseq.trim.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        args_trim = self.__dict__.copy()
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # update args
        args_trim['fq1'] = args_trim['fq'] = fq1
        args_trim['fq2'] = fq2
        args_trim['outdir'] = self.clean_dir

        if trimmed is True:
            ## fq1
            if is_gz(fq1):
                symlink(fq1, clean_fq1, absolute_path=True)
            else:
                gzip_cmd(fq1, clean_fq1, decompress=False, rm=False)
            ## fq2
            if not fq2 is None:
                if is_gz(fq2):
                    symlink(fq2, clean_fq2, absolute_path=True)
                else:
                    gzip_cmd(fq2, clean_fq2, decompress=False, rm=False)

        else:
            if check_file(self.clean_fq_list):
                log.info('trim() skipped, file exists: {}'.format(
                    self.clean_fq_list))
            else:
                Trimmer(**args_trim).run()


    def align(self):
        """
        Alignment PE reads to reference genome, using STAR
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])
        outdir = args_local.get('align_dir', str(pathlib.Path.cwd()))

        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'fq': fq1,
            'outdir': outdir
        }
        args_local.update(args_init)

        # output
        if check_file(args_local.get('bam_raw', None)):
            log.info('align() skipped, file exists: {}'.format(
                args_local.get('bam_raw', None)))
        else:
            Alignment(**args_local).run()


    def get_raw_bam(self):
        """
        Get the align bam file
        from bam_dir/1., 2., ...
        !!! specific: 2.genome/*.bam
        """
        bamlist = listfile(self.align_dir, '*.bam', recursive=True)
        bamlist = [b for b in bamlist if not b.endswith('.raw.bam')]
        # [spike-in]? [rRNA, genome]
        return(bamlist[-1]) # last one


    def save_bam(self):
        """
        Save genome mapped bam file to bam_file/
        from: self.bam_raw, align/smp_name/2.dm6/smp_name.bam
        to: self.bam_out, bam_files/smp_name.bam
        """
        bam_in = self.get_raw_bam() # to genome
        if not os.path.exists(self.bam_out):
            symlink(bam_in, self.bam_out, absolute_path=True)


    def fc_count(self):
        """
        Run FeatureCounts for the bam file
        """
        args = self.__dict__.copy()

        # determine the strandness
        strand, strand_status = RNAseqLibrary(
            bam=self.get_raw_bam(),
            gtf=self.gtf).run(with_status=True)

        # save to status
        with open(self.strandness_status, 'w') as w:
            w.write(strand_status)

        # check strand
        if strand == 'read1':
            strand_fwd, strand_rev = (1, 2)
        elif strand == 'read2':
            strand_fwd, strand_rev = (2, 1)
        else:
            strand_fwd, strand_rev = (0, 0)

        # run sense
        args_fwd = {
            'gtf': self.gtf,
            'bam_list': self.get_raw_bam(),
            'outdir': self.count_dir,
            'strandness': strand_fwd,
            'outname': 'count.sens.txt'}
        _, _, assign_fwd = FeatureCounts(**args_fwd).run()

        # run anti
        args_rev = {
            'gtf': self.gtf,
            'bam_list': self.get_raw_bam(),
            'outdir': self.count_dir,
            'strandness': strand_rev,
            'outname': 'count.anti.txt'}
        _, _, assign_rev = FeatureCounts(**args_rev).run()


    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam_out,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 12,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }

        Bam2bw(**args_local).run()


    def report(self):
        """
        Create alignment report for RNAseq
        1. trim
        2. align
        3. quant
        ...
        """        
        pkg_dir    = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'rnaseq_report.R')
        cmd_file = os.path.join(self.report_dir, 'cmd.sh')
        report_html = os.path.join(
            self.report_dir, 
            'rnaseq_report.html')

        cmd = 'Rscript {} {} {} {}'.format(
            qc_reportR,
            self.outdir,
            self.report_dir,
            self.feature)

        # save cmd
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')

        if check_file(report_html):
            log.info('report() skipped, file exists: {}'.format(
                report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('report() failed.')


    def run(self):
        """
        Run pipeline for RNAseq,
        process
        """
        # sys.exit('!BBBB-4  exit RNAseqSingle()')
        # init dir
        copy_raw_fq = getattr(self, 'copy_raw_fq', False)
        trimmed = getattr(self, 'trimmed', False)        

        # 1. copy raw data
        self.prep_raw(copy_raw_fq)

        # 2. trim
        self.trim(trimmed)

        # 3. align
        self.align()
        self.save_bam() 

        # 4. quant
        self.fc_count()

        # 5. bam2bw
        self.bam_to_bw()

        # 6.report
        self.report()

        return self.outdir


class RNAseqDeseqMultiple(object):
    """
    Run multiple ctl/exp pairs, from design.json file
    deseq_multiple_args = {pairs}
    """
    def __init__(self, **kwargs):
        """
        read design.json
        """
        self.update(kwargs) # fresh new
        self.init_args() # update all variables: *.config, *.args


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
        --design
        """
        local_config = RNAseqConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        assert self.rnaseq_type == 'deseq_multiple'

        # read data from design.json
        self.design_args = Json(self.design).reader()


    def run_pair_single(self, k):
        """
        run 1 RNAseqDeseqSingle() # !!!! parallel-jobs
        input: pair, [index_ctl, index_exp]
        """
        args_k = self.design_args[k] # fq

        ## run locale
        args_local = self.__dict__.copy()
        args_local.pop('rnaseq_type', None)
        # fx = args_local.pop('design_args', None) # multiple 
        # f1 = fx.pop(k, None) # single

        ## ctl
        args_c = args_local.copy()
        args_c['fq1'] = args_k.get('ctl', None)
        args_c['fq2'] = args_k.get('ctl_read2', None)
        args_c['design'] = None
        dirs_ctl = RNAseqMultiple(**args_c).run()

        ## exp
        args_e = args_local.copy()
        args_e['fq1'] = args_k.get('exp', None)
        args_e['fq2'] = args_k.get('exp_read2', None)
        args_e['design'] = None
        dirs_exp = RNAseqMultiple(**args_e).run()

        ## RUN DESeq single
        args_n = args_local.copy()
        args_x = {
            'dirs_ctl': dirs_ctl,
            'dirs_exp': dirs_exp,
            'smp_name': None,
            'group': None,
            'smp_path': None,
            'ctl': None,
            'exp': None,
            'ctl_read2': None,
            'exp_read2': None,
            'design': None,
            'build_design': None}
        args_n.update(args_x)
        RNAseqDeseqSingle(**args_n).run()


    def run(self):
        """
        for group >= 2, paires >= 1
        """
        de_list = list(self.design_args.keys())
        if len(de_list) > 0:
            ## do not allowed, pool with in pool
            # with Pool(processes=self.parallel_jobs) as pool:
            #     pool.map(self.run_pair_single, list(self.design_args.keys()))
            for k in de_list:
                self.run_pair_single(k)
        else:
            log.warning('groups >2 expected, {} found'.format(self.group))


class RNAseqDeseqSingle(object):
    """
    Run RNAseq for single pair:

    required args:
        -dirs_ctl
        -dirs_exp
        -genme
        -outdir
        -feature

    config
    count/txt, rpkm, up, down, criteria
    de
    pdf/ma, pca, volcano
    report/html
    enrich/GO, GSEA, ...
    """
    def __init__(self, **kwargs):
        """
        dirs_ctl
        dirs_exp
        genome
        outdir
        feature
        group_ctl (None)
        group_exp (None)
        merge replicates
        """
        self.update(kwargs) # fresh new
        self.init_args() # update all variables: *.config, *.args


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
        local_config = RNAseqConfig(**self.__dict__) # update, global
        self.update(local_config.__dict__, force=True) # update local attributes
        assert self.rnaseq_type == 'deseq_single'

        # check arguments
        chk1 = args_checker(self.__dict__, self.config_pickle)
        # Json(self.__dict__).writer(self.config_json)
        chk2 = args_logger(self.__dict__, self.config_txt)

    #######################################
    ## main pipeline
    def copy_count_files(self):
        """
        Copy count files to target dirs
        rename the count files by fqname
        """
        # copy count files
        for n, f in zip(self.smp_name, self.count_ctl + self.count_exp):
            f_new = os.path.join(self.count_dir, n + '.' + os.path.basename(f))
            shutil.copy(f, f_new)


    def deseq2(self):
        """
        ## quality control
        ## PCA/Cor/...
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        deseqR = os.path.join(pkg_dir, 'bin', 'run_deseq2.R')
        cmd_file = os.path.join(self.deseqdir, 'cmd.sh')
        cmd = 'Rscript {} {} {}'.format(
            deseqR,
            os.path.join(self.outdir, self.project_name),
            self.feature) #,
            # os.path.join(self.deseqdir, 'mylog.deseq2.out'))  # pvalue cutoff

        # save cmd
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')

        try:
            run_shell_cmd(cmd)
        except:
            log.warning('DESeq2 failed.')


    def report(self):
        """
        make report for deanalysis
        """
        pass


    def run(self):
        """
        Run all
        """
        # 1. copy count
        self.copy_count_files()

        # # 2. design
        # self.get_design()

        # 3. run DESeq2
        self.deseq2()

        # 4. report to file/qc
        self.report()


class RNAseqLibrary(object):
    """
    Determine the RNAseq library type:
    strandness: 1 ++, 1 --, / 2 +-, 2 -+ :
    dUTP, NSR: 1 +-, 1 -+, / 2 ++, 2 -- :
    ###
    or:
    infer_experiment.py from RSeQC package.
    """
    def __init__(self, bam, gtf, size=200000, **kwargs):
        self.update(kwargs) # fresh new
        self.bam = bam
        self.gtf = gtf
        self.size = size
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
        Prepare Bam file (subset)
        """
        args_default = {
            'outdir': self._tmp(),
            'cleanup': True
        }
        self.update(args_default, force=False)

        # input bam
        if isinstance(self.bam, str):
            pass
        elif isinstance(self.bam, list):
            self.bam = self.bam[0] # the first one
        else:
            raise Exception('Expect str and list, get {}'.format(type(bam)))

        # outdir
        if self.outdir is None:
            self.outdir = self._tmp()
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # subset bam
        self.bam_sub = self.subset_bam(self.outdir, self.size) # subset

        ## input a BED file, convert to GTF
        if self.gtf.lower().endswith('.bed'):
            self.bed = self.gtf
            self.gtf = os.path.join(self.outdir,
                os.path.basename(os.path.splitext(self.gtf)[0]) + '.gtf')
            # convert
            self.bed2gtf(self.bed, self.gtf)


    def _tmp(self):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=True)
        return tmp.name


    def bed2gtf(self, bed_in, gtf_out):
        """
        Convert BED to GTF format
        """
        with open(bed_in) as r, open(gtf_out) as w:
            for line in r:
                fields = line.rstrip().split('\t')
                start = int(fields[1]) + 1 # 0-index to 1-index
                end = fields[2]
                des = 'gene_id "{}"; gene_name "{}";'.format(fields[3], fields[3])
                gtf_line = '\t'.join([
                    fields[0],
                    'BED_file',
                    'exon',
                    str(start),
                    end,
                    '.',
                    fields[5],
                    '.',
                    des])
                w.write(gtf_line + '\n')


    def subset_bam(self, dest, size=20000):
        """
        Extract N reads from bam list
        """
        if not os.path.exists(dest):
            os.mkdir(dest)

        ## create bam index
        Bam(self.bam).index()
        bname = os.path.basename(self.bam)
        dest_bam = os.path.join(dest, bname)
        samfile = pysam.AlignmentFile(self.bam, 'rb')
        subfile = pysam.AlignmentFile(dest_bam, 'wb', template=samfile)

        # count
        n = 0
        for read in samfile.fetch():
            n += 1
            if n > size:
                break
            subfile.write(read)

        return dest_bam


    def stat(self, strandness=1):
        """
        strandness; 1 or 2; for featureCounts
        using -s 1, -s 2
        and return the mapping reads
        """
        ## -s 1
        outname = 'count.s_' + str(strandness) + '.txt'
        args = {
            'gtf': self.gtf,
            'bam_list': self.bam_sub,
            'outdir': self.outdir,
            'strandness': strandness,
            'outname': outname}

        total, assign, assign_df = FeatureCounts(**args).run()

        return assign
        

    def run(self, with_status=False):
        """
        Check -s 1, -s 2:
        """
        a = self.stat(strandness=1).to_list()[0] # s=1
        b = self.stat(strandness=2).to_list()[0] # s=2

        ## loginfo
        strand_status = 'Library strandness by featureCounts: -s\n'
        strand_status += 'sample: {}\n'.format(self.bam)
        strand_status += 'gtf: {}\n'.format(self.gtf)
        strand_status += 'assigned_fwd: {:>6.2f}%  -s=1  stranded\n'.format(a)
        strand_status += 'assigned_rev: {:>6.2f}%  -s=2  reversely_stranded\n'.format(b)
        log.info(strand_status)

        ## correct BED
        if a + b < 20:
            log.warning('check out gtf file: {}'.format(self.gtf))
            tag = 'undetermined'
        elif a == b:
            tag = 'read0'
        elif a > b:
            tag = 'read1'
        else:
            tag = 'read2'

        if self.cleanup is True:
            # check ?
            log.warning('Temporary directory in: {}'.format(self.outdir))
            if self.outdir.startswith('/tmp/'):
                # remove files within /tmp folder # linux
                shutil.rmtree(self.outdir) # dangerous !!!

        if with_status:
            return (tag, strand_status)
        else:
            return tag


class FeatureCounts(object):
    """
    Run featureCounts for GTF + BAM(s)
    mission()
    count.txt
    determin: library type: +-/-+; -+/+-;
    map reads, scale
    FPKM/RPKM;
    map pct:
    """
    def __init__(self, **kwargs):
        """
        arguments:
        gtf:
        bam_list: (index)
        outdir:
        strandness: 0=no, 1=sens, 2=anti, 3=both
        smpname: count.txt default
        threads: 4
        overwrite: False
        """
        # required args
        self.update(kwargs) # fresh new
        # self.args = kwargs
        self.check_status = self.init_fc() # global config


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


    def init_fc(self):
        """
        Initiate the config for RNAseq analysis
        """
        args_default = {
            'gtf': None,
            'bam_list': None,
            'outdir': str(pathlib.Path.cwd()),
            'strandness': 0,
            'threads': 4,
            'overwrite': False,
            'outname': 'count.txt'
        }
        self.update(args_default, force=False)

        # required
        check_path(self.outdir)

        # convert to list
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list]

        # check required args
        chk1 = isinstance(self.gtf, str) and os.path.exists(self.gtf)
        chk2 = all([os.path.exists(i) for i in self.bam_list])
        if not chk1 or not chk2:
            log.error('fc() failed, check gtf={}, bam_list={}'.format(self.gtf, self.bam_list))

        ## absolute path
        self.gtf = file_abspath(self.gtf)
        self.bam_list = file_abspath(self.bam_list)
        self.outdir = file_abspath(self.outdir)

        ## index bam files
        [Bam(i).index() for i in self.bam_list]

        # determine the output files
        self.count_txt = os.path.join(self.outdir, self.outname)
        self.summary = self.count_txt + '.summary'
        self.log_file = os.path.join(self.outdir, 'featureCounts.log')

        ## PE reads
        self.pe_flag = self.is_PE_bam()

        ## tool
        self.fc_exe = shutil.which('featureCounts')

        ## update args
        args = self.__dict__
        args.pop('args_input', None)
        args.pop('cmd_input', None)
        args.pop('dict', None)
        args.pop('args', None)

        ## save config
        args_txt = os.path.join(self.outdir, 'arguments.txt')
        args_json = os.path.join(self.outdir, 'arguments.json')
        args_pickle = os.path.join(self.outdir, 'arguments.pickle')
        # save to file

        # Json(args).writer(args_json)
        args_logger(args, args_txt)
        chk3 = args_checker(args, args_pickle); chk3 = True
        chk4 = self.overwrite

        # status
        return all([chk3, chk4])


    def get_cmd(self):
        """
        prepare args for featureCounts
        """
        cmd = '{} -M -O --fraction -g gene_id -t exon '.format(
            self.fc_exe)

        # for PE
        if self.pe_flag:
            cmd += ' -p -C -B '

        # args
        cmd += '-s {} -T {} -a {} -o {} {} 2> {}'.format(
            self.strandness,
            self.threads,
            self.gtf,
            self.count_txt,
            ' '.join(self.bam_list),
            self.log_file)

        return cmd


    def is_PE_bam(self):
        """
        Check whether the input bam files are Paired or Single file
        Bam().isPaired()
        """
        return all([Bam(i).isPaired() for i in self.bam_list])


    def wrap_log(self):
        """
        save output file to log,
        """
        df = pd.read_csv(self.summary, '\t', index_col=0)
        df.columns = [os.path.basename(i) for i in df.columns.tolist()]
        ## total
        total = df.sum(axis=0, skipna=True)
        ## assign
        assign = df.loc['Assigned', ] / total * 100
        ## pct
        assign_df = assign.to_frame('assigned%')
        ## minimum
        min_pct = assign.min()
        if assign.min() < 50:
            log.warning('caution: -s {}, {:.2f}% reads were assigned, {}'.format(
                self.strandness,
                assign.min(),
                self.summary))
        return (total, assign, assign_df)


    def _tmp(self):
        """
        Create a tmp filename
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', delete=True)
        return os.path.basename(tmp.name)


    def run(self):
        """
        run featureCounts
        """
        cmd = self.get_cmd()
        if os.path.exists(self.count_txt) and self.overwrite is False:
            log.info('FeatureCounts() skipped, file exists: {}'.format(
            self.count_txt))
        else:
            run_shell_cmd(cmd)
            
        return self.wrap_log()


class RNAseqReader(object):
    """
    Return the config, files for RNAseq direcotory

    RNAseqSingle: config/
    RNAseqMultiple: config/
    DeseqSingle: config/

    # to-do
    DEseqMultiple: ? (RNAseqMultiple)
    """
    def __init__(self, path, **kwargs):
        self.path = path
        self.feature = kwargs.get('feature', 'gene')
        self.check()


    def update(self, d, force=True, remove=False):
        """
        d: dict
        all: bool, update exists attributes
        remove: bool, remove exists attributes
        Update attributes from dict
        all exists attr
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


    def check(self, return_args=False):
        """
        Check if the directory is of RNAseq 
        [single|multiple|deseq_single|deseq_multiple]
        single - config:
            outdir/feature/config/arguments.pickle
        multiple - config:
            outdir/config/feature/config/arguments.pickle
        deseq - single:
            outdir/feature/config/arguments.pickle
        """
        x1 = os.path.join(self.path, self.feature, 'config', 'arguments.pickle')
        x2 = os.path.join(self.path, 'config', self.feature, 'config', 'arguments.pickle')

        if os.path.exists(x1):
            args = pickle_to_dict(x1)
        elif os.path.exists(x2):
            args = pickle_to_dict(x2)
        else:
            args = {}
            log.error("""
                unknown directory, expect config file:
                RNAseq single: {}
                RNAseq multiple: {}
                RNAseq deseq: {}""".format(x1, x2, x1))
        
        # update self
        self.update(args, force=True) # fresh new


    def is_rnaseq_single(self):
        """
        Check if the directory is of RNAseq single sample
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'rnaseq_single'


    def is_rnaseq_multi(self):
        """
        Check if the directory is of RNAseq multiple sample
        outdir/config/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'rnaseq_multiple'


    def is_deseq_single(self):
        """
        Check if the directory if of RNAseq DESeq analysis
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'deseq_single'


    def is_deseq_multiple(self):
        """
        Check if the directory if of RNAseq DESeq analysis
        outdir/feature/config/arguments.pickle
        """
        return self.rnaseq_type == 'deseq_multiple'


