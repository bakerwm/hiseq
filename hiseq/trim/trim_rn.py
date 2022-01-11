#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: multiple fx

1. Trim adapter,
2. remove dup
3. trim n-bases from read

3.1 TruSeq (NSR)  
    - cut 9, -9; from both ends of read (subseq)    

3.2 TruSeq (iCLIP) 
    - cut 9; from read1    

3.3 TruSeq (eCLIP)  
    - cut 10; -7 from read1  
    - cut 7, -10 from read2

optional
1. --rm-untrim, --save-too-short, ...
"""

import os
import shutil
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.utils.file import (
    check_path, file_abspath, check_fx, check_fx_args, fx_name
)
from hiseq.utils.utils import log, update_obj, Config
from hiseq.trim.trim_r1 import TrimR1, get_args_trim


class TrimRn(object):
    """Processing fastq files: N (multi)
    Trim adapters for multiple fastq files
    should be SE or PE
    """
    def __init__(self, **kwargs):        
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
#             'library_type': None,
            'fq1': None,
            'fq2': None,
            'smp_name': None,
#             'outdir': None,
#             'rmdup': False,
#             'cut_after_trim': None,
#             'cut_to_length': 0,
#             'recursive': False, # for ATACseq, smRNA, adapter-sliding,
#             'keep_tmp': False,
#             'overwrite': False,
#             'adapter3': None, # read1
#             'adapter5': None, # read1
#             'Adapter3': None, # read2
#             'Adapter5': None, # read2
#             'qual_min': 20,
#             'error_rate': 0.1,
#             'overlap': 3,
#             'percent': 80,
#             'rm_polyN': False,
#             'rm_untrim': False,
#             'save_untrim': False,
#             'save_too_short': False,
#             'save_too_long': False,
#             'cut_before_trim': 0,
#             'discard_tooshort': False,
#             'overwrite': False,
#             'parallel_jobs': 1,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_rn' #
        self.is_paired = self.fq2 is not None
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        # init smp_name
        if isinstance(self.smp_name, list) and len(self.fq1) == len(self.smp_name):
            pass
        else:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.rep_list = [os.path.join(
            self.outdir, i) for i in self.smp_name]
        self.init_files()
        self.init_fq()
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_fq(self): # se or pe
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 no valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
                

    def init_files(self):
        self.project_dir = self.outdir
        self.config_dir = self.project_dir + '/config'
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path(self.config_dir)
                
                
    def run_r1(self, f1):
        args_local = self.__dict__.copy()
        i = self.fq1.index(f1) # index
        args_local.update({
            'fq1': f1,
            'fq2': self.fq2[i] if self.is_paired else None,
            'outdir': self.outdir,
            'smp_name': self.smp_name[i],
        })
        TrimR1(**args_local).run()


    def run(self):
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_r1, self.fq1)
        else:
            [self.run_r1(fq) for fq in self.fq1]

            
def get_args():
    return get_args_trim(get_args_io())
    
            
def get_args_io():
    example = '\n'.join([
        'Examples:',
        '1. auto-detect adapters',
        '$ python trim.py -1 fq1 -2 fq2 -o outdir -m 20',
        '2. specific 3-apdapter',
        '$ python trim.py -1 fq1 -2 fq2 -o outdir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT',
    ])
    parser = argparse.ArgumentParser(
        prog='trim',
        description='trim adapters',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='The read2 of pair-end reads')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')
    return parser
            

def main():
    args = vars(get_args().parse_args())
    TrimRn(**args).run()


if __name__ == '__main__':
    main()

#