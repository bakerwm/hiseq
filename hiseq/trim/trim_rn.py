#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: multiple fx

1. Trim adapter,
2. remove dup
3. trim n-bases from read

3.1 TruSeq (NSR)  
    - cut 7, -7; from both ends of read (subseq)    

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
from hiseq.utils.file import check_path, file_abspath, check_fx, \
    check_fx_paired, fx_name
from hiseq.utils.utils import log, update_obj, Config
from hiseq.trim.trim_r1 import TrimR1


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
            'library_type': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'rmdup': False,
            'smp_name': None,
            'cut_after_trim': None,
            'cut_to_length': 0,
            'recursive': False, # for ATACseq, smRNA, adapter-sliding,
            'keep_tmp': False,
            'overwrite': False,
            'adapter3': None, # read1
            'adapter5': None, # read1
            'Adapter3': None, # read2
            'Adapter5': None, # read2
            'qual_min': 20,
            'error_rate': 0.1,
            'overlap': 3,
            'percent': 80,
            'rm_polyN': False,
            'rm_untrim': False,
            'save_untrim': False,
            'save_too_short': False,
            'save_too_long': False,
            'cut_before_trim': 0,
            'discard_tooshort': False,
            'overwrite': False,
            'parallel_jobs': 1,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_rn' #
        self.is_paired = self.fq2 is not None
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        # init smp_name
        # if self.smp_name is None or not len(self.fq1) == len(self.smp_name):
        if isinstance(self.smp_name, list) and len(self.fq1) == len(self.smp_name):
            pass
        else:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        self.init_fq()
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_fq(self): # se or pe
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        c1 = isinstance(self.fq1, list)
        c1q = all(check_fx(self.fq1)) # force
        if self.fq2 is None:
            c2 = c2q = c2pe = False
            c2x = True
        else:
            c2 = isinstance(self.fq2, list)
            c2q = all(check_fx(self.fq2)) # force
            c2pe = all(check_fx_paired(self.fq1, self.fq2))
            c2x = all([c2, c2q, c2pe])
        # out
        c0 = all([c1, c1q, c2x])
        # show message
        msg = '\n'.join([
            '='*80,
            'Run TrimRn() with parameters:',
            '{:>14} : {}'.format('fq1', self.fq1),
            '{:>14} : {}'.format('fq2', self.fq2),
            '-'*40,
            '{:>14} : {}'.format('fq1 [list]', c1),
            '{:>14} : {}'.format('fq1 [fastq]', c1q),
            '{:>14} : {}'.format('fq2 [list]', c2),
            '{:>14} : {}'.format('fq2 [fastq]', c2q),
            '{:>14} : {}'.format('fq [paired]', c2pe),
            '{:>14} : {}'.format('rmdup', self.rmdup),
            '{:>14} : {}'.format('cut after trim', self.cut_after_trim),
            '-'*40,
            'Status: {}'.format(c0),
            '='*80,
        ])
        # print(msg)
        if not c0:
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
    example = '\n'.join([
        'Examples:',
        '1. auto-detect adapters',
        '$ python trim.py -1 fq1 -2 fq2 -o outdir -m 20',
        '2. specific 3-apdapter',
        '$ python trim.py -1 fq1 -2 fq2 -o outdir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT',
    ])
    parser = argparse.ArgumentParser(
        prog='trim',
        description='trim: for multiple reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='The read2 of pair-end reads')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')

    parser.add_argument('--library-type', dest='library_type', default=None,
        type=str, choices=['TruSeq', 'Nextera', 'smallRNA'],
        help='Type of the library structure, \
        TruSeq, TruSeq standard library \
        Nextera, Tn5 standard library, \
        smallRNA, small RNA library')
    parser.add_argument('-m', '--len_min', default=15, metavar='len_min',
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('--cut-to-length', default=0, dest='cut_to_length',
        type=int,
        help='cut reads to from right, default: [0], full length')
    parser.add_argument('--rmdup', action='store_true',
        help='remove duplicates')
    parser.add_argument('--cut-after-trim', dest='cut_after_trim',
        help='cut reads after triming, positive value, \
        cut from left; minus value, cut from right, eg: 3 or -4 or 3,-4, \
        default [0]')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')    
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1, type=int,
        help='Number of jobs run in parallel, default: [1]')

    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    ## global arguments
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', '--overlap', type=int, default=3,
        help='At least overlap between adapter and sequence, default: [3]')
    ## specific
    parser.add_argument('--rm-polyN', action='store_true',
        dest='rm_polyN',
        help='remove polyN')
    parser.add_argument('--rm-untrim', action='store_true',
        dest='rm_untrim',
        help='discard reads without adapter')
    parser.add_argument('--save-untrim', action='store_true',
        dest='save_untrim',
        help='Save untrim reads to file')
    parser.add_argument('--save-too-short', action='store_true',
        dest='save_too_short',
        help='Save too short reads to file')
    parser.add_argument('--save-too-long', action='store_true',
        dest='save_too_long',
        help='Save too short reads to file')
    parser.add_argument('--cut-before-trim', default='0',
        dest='cut_before_trim',
        help='cut n-bases before trimming adapter; positive value, \
        cut from left; minus value, cut from right, eg: 3 or -4 or 3,-4, \
        default [0]')

    parser.add_argument('-a', '--adapter3', default=None,
        help='3-Adapter sequence, default [].')
    parser.add_argument('-g', '--adapter5', default='',
        help='5-Adapter, default: None')

    ## PE arguments
    parser.add_argument('-A', '--Adapter3', default=None,
        help='The 3 adapter of read2, default []')
    parser.add_argument('-G', '--Adapter5', default=None,
        help='The 5 adapter of read1, default: None')
    return parser


def main():
    args = vars(get_args().parse_args())
    TrimRn(**args).run()


if __name__ == '__main__':
    main()

#