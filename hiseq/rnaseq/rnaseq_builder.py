
"""
Build RNAseq experiment, from fastq files
specify ctl, exp (fastq files)

"""


import os
# import re
# import hiseq
# import pysam
# import shutil
import tempfile
from Levenshtein import distance
# import pandas as pd
# from multiprocessing import Pool
# import copy # copy objects
from hiseq.utils.helper import *



class RNAseqFqDesign(object):
    """
    Prepare ctl/exp samples for RNAseq analysis
    append/update

    format:
    - ctl: list, rep1, rep2, ...
    - exp: list, rep1, rep2, ...
    """
    def __init__(self, **kwargs):
        self.update(kwargs)
        self.init_args() # init fq
        self.fq_args() # check fq
        self.save() # to file


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
            'ctl': None,
            'exp': None,
            'ctl_read2': None,
            'exp_read2': None,
            'design': None,
            'append': True
        }

        self.update(args_default, force=False)

        # check args
        ## 1. ctl should be list
        self.ctl = self.fq_input(self.ctl)
        self.exp = self.fq_input(self.exp)
        ## 2. read2
        if not self.ctl_read2 is None:
            self.ctl_read2 = self.fq_input(self.ctl_read2)
        if not self.exp_read2 is None:
            self.exp_read2 = self.fq_input(self.exp_read2)
        ## 3. design
        if self.design is None:
            raise Exception('--design x.json requrid')


    def fq_input(self, fq1):
        """
        require fq input:
        list, 
        paired in PE
        """
        flag = False

        # chk
        if fq1 is None:
            log.error('failed, fastq files, list expected, None found')
        elif isinstance(fq1, str):
            log.warning('list expected, str found')
        elif isinstance(fq1, list):
            chk1 = file_exists(fq1)
            # show log
            log.info('Check file exists:')
            for fx, chk in zip(fq1, chk1):
                print('{} : {}'.format(chk, fx))

            if all(chk1):
                flag = file_abspath(fq1) 
            else:
                log.error('failed fastq files error')
        else:
            log.error('failed, list expected, {} found'.format(type(fq1).__name__))

        if not flag:
            raise Exception('fastq files failed')

        return flag


    def fq_args(self):
        """
        Check fq files
        fq1, required
        fq2, optional

        file exists
        fq1 = fq2
        """
        flag = 0

        # length
        chk1 = len(self.ctl) == len(self.exp)

        # fq2
        if not self.ctl_read2 is None:
            chk2a = len(self.ctl) == len(self.ctl_read2)
            chk2b = [self.fq_pair(a, b) for a, b in zip(self.ctl, self.ctl_read2)]
            # chkeck pair
            if not (chk2a and all(chk2b)): 
                flag += 1
                log.error('failed ctl: {}, ctl_read2: {}, files not matched'.format(
                    len(self.ctl), len(self.ctl_read2)))

        if not self.exp_read2 is None:
            chk3a = len(self.exp) == len(self.exp_read2)
            chk3b = [self.fq_pair(a, b) for a, b in zip(self.exp, self.exp_read2)]
            # check pair
            if not (chk3a and all(chk3b)): 
                flag += 1
                log.error('failed exp: {}, exp_read2: {}, files not matched'.format(
                    len(self.exp), len(self.exp_read2)))

        if flag > 0:
            raise Exception('check fastq files')


    def fq_pair(self, fq1, fq2, n_diff=1):
        """
        read1 and read2 
        should be in the same name, _1, _2
        """
        if isinstance(fq1, str) and isinstance(fq2, str):
            if distance(fq1, fq2) <= n_diff:
                return True


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.json',
            delete=False)
        return tmp.name


    def save(self):
        """
        save config
        """
        # json file
        if not file_exists(self.design):
            args_pre = {}
        else:
            args_pre = Json(self.design).reader()

        # save only specific args
        dd = {'ctl': self.ctl,
              'exp': self.exp,
              'ctl_read2': self.ctl_read2,
              'exp_read2': self.exp_read2,
              'design': self.design,
              'apend': self.append}

        # new data
        # append or update
        if self.append:
            key = 'rnaseq_{:03d}'.format(len(args_pre) + 1)
            # args_input = self.__dict__ # new
            # whether exists
            if dd in list(args_pre.values()):
                log.warning('design exists, skipped ...')
            else:
                args_pre[key] = dd # append
        else:
            key = 'rnaseq_001'
            args_pre = {} # empty
            args_pre[key] = dd #

        # save to file or to stdout
        if self.design is None:
            f = self._tmp()
            Json(args_pre).writer(f)
        else:
            Json(args_pre).writer(self.design)


## Test

# fq = [
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_1.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_2.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep2_1.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep2_2.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_treatment_rep1_1.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_treatment_rep1_2.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_treatment_rep2_1.fq.gz',
# '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_treatment_rep2_2.fq.gz'
# ]
        
# args = {
#     'ctl': [fq[0], fq[2]],
#     'exp': [fq[4], fq[6]],
#     'ctl_read2': [fq[1], fq[3]],
#     'exp_read2': [fq[5], fq[7]],
#     'append': True,
#     'design': 'test.json'
# }

# RNAseqFqDesign(**args)