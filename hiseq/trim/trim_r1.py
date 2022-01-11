#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: single fx

1. Trim adapter,
2. remove pcr dup
3. trim n-bases from read

3.1 TruSeq (NSR)  
    - cut 7,-7; from both ends of read (subseq)    

3.2 TruSeq (iCLIP) 
    - cut 9; from read1    

3.3 TruSeq (eCLIP)  
    - cut 10,-7 from read1  
    - cut 7,-10 from read2
"""

import os
import re
import shutil
import pathlib
import argparse
import pandas as pd
import pyfastx
from xopen import xopen
from hiseq.utils.seq import Fastx
from hiseq.utils.file import (
    check_path, check_fx_args, file_exists, file_abspath,
    check_fx, check_fx_paired, copy_file, symlink_file, remove_file, fx_name
)
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd, get_date
from hiseq.trim.cutadapt import Cutadapt


class TrimR1(object):
    """
    for single fastq file
    1. trim adapter
    2. remove dup
    3. cut n-bases
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args_local = self.__dict__.copy()
        ob = TrimR1Config(**args_local)
        self = update_obj(self, ob.__dict__.copy(), force=True)
        Config().dump(self.__dict__.copy(), self.config_yaml)
        
    
    def get_msg(self):
        msg = '\n'.join([
            '-'*80,
            '{:>20s} : {}'.format('Date', get_date()),
            '{:>20s} : {}'.format('program', 'Trim'),
            '{:>20s} : {}'.format('config', self.config_yaml),
            '{:>20s} : {}'.format('fq1', self.fq1),
            '{:>20s} : {}'.format('fq2', self.fq2),
            '{:>20s} : {}'.format('outdir', self.outdir),
            '{:>20s} : {}'.format('library_type', self.library_type),
            '{:>20s} : {}'.format('len_min', self.len_min),
            '{:>20s} : {}'.format('len_max', self.len_max),
            '{:>20s} : {}'.format('rmdup', self.rmdup),
            '{:>20s} : {}'.format('cut_before_trim', self.cut_before_trim),
            '{:>20s} : {}'.format('cut_after_trim', self.cut_after_trim),            
            '{:>20s} : {}'.format('3-adapter (read1)', self.adapter3),
            '{:>20s} : {}'.format('3-adapter (read2)', self.Adapter3),
            '{:>20s} : {}'.format('5-adapter (read1)', self.adapter5),
            '{:>20s} : {}'.format('5-adapter (read2)', self.Adapter5),
            '{:>20s} : {}'.format('remove poly(N)', self.rm_polyN),
            '{:>20s} : {}'.format('remove untrim', self.rm_untrim),
            '-'*80,
        ])
        return msg
    
    
    def run(self):
        trim_args = self.__dict__.copy()
        print(self.get_msg())
        Trim = TrimR1pe if self.is_paired else TrimR1se
        Trim(**trim_args).run()
    

class TrimR1Config(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'library_type': None, # auto
            'fq1': None, # str
            'fq2': None, # str
            'outdir': None, # str
            'smp_name': None,
            'len_min': 15,
            'len_max': 0,
            'threads': 4,
            'rmdup': False,
            'cut_after_trim': 0, # skip, str
            'cut_before_trim': 0,
            'cut_to_length': 0,
            'keep_tmp': False,
            'overwrite': False,
            'recursive': False, # for ATACseq, smRNA, adapter-sliding,
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
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_r1' #
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        self.init_fq()
        self.init_files()
        # fix cut_after_trim
        if self.cut_after_trim is None:
            self.cut_after_trim = '0' # default
        self.cut2_l, self.cut2_r = self.parse_cut2_arg(self.cut_after_trim)
        self.cut2_skip = (self.cut2_l + abs(self.cut2_r) == 0)


    def init_fq(self):
        """
        required:
        1. fq1:str, fq2:str/None
        2. paired
        """
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not valide, check above message')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)

        
    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        self.config_dir = self.project_dir + '/config'
        self.cutadapt_dir = self.project_dir + '/tmp/01_cutadapt'
        self.rmdup_dir = self.project_dir + '/tmp/02_rmdup'
        self.cut2_dir = self.project_dir + '/tmp/03_cut_after_trim'
        self.log_dir = self.project_dir + '/log'
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'rmdup_fq': self.rmdup_dir + '/' + self.smp_name + '.fq.gz',
            'rmdup_fq1': self.rmdup_dir + '/' + self.smp_name + '_1.fq.gz',
            'rmdup_fq2': self.rmdup_dir + '/' + self.smp_name + '_2.fq.gz',
            'cut2_fq': self.cut2_dir + '/' + self.smp_name + '.fq.gz',
            'cut2_fq1': self.cut2_dir + '/' + self.smp_name + '_1.fq.gz',
            'cut2_fq2': self.cut2_dir + '/' + self.smp_name + '_2.fq.gz',
            'clean_fq': self.project_dir + '/' + self.smp_name + '.fq.gz',
            'clean_fq1': self.project_dir + '/' + self.smp_name + '_1.fq.gz',
            'clean_fq2': self.project_dir + '/' + self.smp_name + '_2.fq.gz',
            'trim_stat': self.project_dir + '/' + self.smp_name + '.trim.stat',
            'trim_json': self.project_dir + '/' + self.smp_name + '.trim.json',
            'trim_yaml': self.project_dir + '/' + self.smp_name + '.trim.yaml',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path([
            self.config_dir, 
            self.cutadapt_dir, 
            self.rmdup_dir, 
            self.cut2_dir, 
            self.log_dir,
        ])

    
    def is_valid_cut2_arg(self, s):
        """
        cut-after-trim: str
        for 1-indexed
        '10': [10:]
        '-9': [:-9]
        '9,-9': [9:-9]
        """
        # sanitize arg
        try:
            s = str(s)
        except:
            log.error('unknown --cut-after-trim expect str, got {}'.format(
                type(s).__name__))
        if isinstance(s, str):
            s = s.strip().replace(' ', '') # remove white spaces
            p = re.compile('^(\d+)$|^(-\d+)$|^((\d+),(-\d+))$')
            m = p.search(s) # match
            out = isinstance(m, re.Match)
        else:
            out = False
        return out

    
    def parse_cut2_arg(self, s):
        """
        for 1-indexed
        cut2: '10', '9,-9', '5', ...
        '10': [10:]
        '-9': [:-9]
        '9, -9': [9:-9]
        make sure: valid 
        return: start, end
        """
        if self.is_valid_cut2_arg(s):
            s = str(s) # to str
            s = s.strip().replace(' ', '') # remove white spaces
            p = re.compile('^(\d+)$|^(-\d+)$|^((\d+),(-\d+))$')
            m = p.search(s) # match [L, R, L&R, L, R]
            g = m.groups()
            if g[0]: # left end
                s = int(g[0])
                e = 0
            elif g[1]: # right end
                s = 0
                e = int(g[1])
            elif g[2]: # both ends
                s = int(g[3])
                e = int(g[4])
            else:
                pass
        else:
            log.error('invalid --cut-after-trim, expect: 5; 5,-3; -3, \
                got {}'.format(self.cut_after_trim))
            raise ValueError('check --cut-after-trim')
        return (s, e)


class TrimR1se(object):
    """
    for single fastq file: Single-End (SE)
    1. trim adapter
    2. remove dup
    3. cut n-bases
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args_local = self.__dict__.copy()
        local_obj = TrimR1Config(**args_local)
        self = update_obj(self, local_obj, force=True)
        # Config().dump(self.__dict__.copy(), self.config_yaml)
        
        
    def cutadapt(self):
        """
        Cut adapters
        """
        cut1_args = self.__dict__.copy()
        cut1_args['outdir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**cut1_args)
        cut1.run() # 
        cut1_fq = cut1.clean_fq
        n_total, n_cut1, _ = cut1.parse_log()
        copy_file(cut1.log, self.log_dir) # save log
        return (cut1_fq, n_total, n_cut1)
    
    
    def rm_pcr_dup(self, fq_in, fq_out):
        if isinstance(fq_in, str) and isinstance(fq_out, str):
            if self.rmdup:
                if file_exists(fq_out) and not self.overwrite:
                    log.info('rmdup() skipped, file exists')
                elif file_exists(fq_in):
                    Fastx(fq_in).collapse(fq_out, fq_out=True)
                else:
                    log.error('run_rmdup() failed, file not exists: {}'.format(fq_in))
            else:
                symlink_file(fq_in, fq_out, absolute_path=False)

    
    def cut2(self, fq_in, fq_out):
        """
        --cut-after-trim
        """
        # sub-function
        # make sure
        # s: str
        # l: int, >= 0
        # r: int, <= 0
        def substr(s, l, r):
            s1 = s[l:]
            if r < 0:
                s1 = s1[:r]
            return s1
        # run
        n_cut2_rm = 0
        if not self.cut2_skip:
            try:
                with xopen(fq_out, 'wt') as w:
                    for name,seq,qual,comment in pyfastx.Fastx(fq_in):
                        if len(seq) < sum([self.len_min, self.cut2_l, abs(self.cut2_r)]):
                            n_cut2_rm += 1
                            continue # drop short seq
                        # update name
                        if comment:
                            name = '{} {}'.format(name, comment)
                        # sub
                        s = substr(seq, self.cut2_l, self.cut2_r)
                        q = substr(seq, self.cut2_l, self.cut2_r)
                        # output
                        w.write('@{}\n{}\n+\n{}\n'.format(name, s, q))
            except IOError as e:
                log.error(e)
        else:
            log.info('cut_after_trim(), skipped')
            symlink_file(fq_in, fq_out, absolute_path=False)
        return n_cut2_rm
    
    
    def run(self):
        ## re-run
        if file_exists(self.clean_fq):
            log.info('Trim() skipped, file exists: {}'.format(self.clean_fq))
            return None
        # 1. cut adapter
        cut1_fq, n_total, n_cut1 = self.cutadapt()
        # 2. remove PCR dup
        self.rm_pcr_dup(cut1_fq, self.rmdup_fq)
        if os.path.islink(self.rmdup_fq):
            n_rmdup = n_cut1
        else:
            n_rmdup = Fastx(self.rmdup_fq).number_of_seq()
        # 3. cut, further adapters
        n_cut2_rm = self.cut2(self.rmdup_fq, self.cut2_fq)
        n_cut2 = n_rmdup - n_cut2_rm
        # 4. wrap files
        if n_total < 1:
            n_total = 1
        n_cut1_rm = n_total - n_cut1
        n_rmdup_rm = n_cut1 - n_rmdup
        n_clean_pct = '{:.2f}'.format(n_cut2/n_total*100)
        # header
        header = ['#name', 'total', 'too_short', 'dup', 'too_short2', 
                  'clean', 'percent']
        s = [self.smp_name, n_total, n_cut1_rm, n_rmdup_rm, n_cut2_rm, n_cut2,
            n_clean_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_cut1_rm),
            'dup': int(n_rmdup_rm),
            'too_short2': int(n_cut2_rm),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # save fq files, cutadapt log
        if file_exists(self.cut2_fq) and not file_exists(self.clean_fq):
            copy_file(self.cut2_fq, self.clean_fq)
        # 5. remove temp files
        del_list = [cut1_fq, self.rmdup_fq, self.cut2_fq]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)
    


class TrimR1pe(object):
    """
    for single fastq file: Paired-End (PE)
    1. trim adapter
    2. remove dup
    3. cut n-bases
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args_local = self.__dict__.copy()
        local_obj = TrimR1Config(**args_local)
        self = update_obj(self, local_obj, force=True)
        # Config().dump(self.__dict__.copy(), self.config_yaml)
        
        
    def cutadapt(self):
        """
        Cut adapters
        """
        cut1_args = self.__dict__.copy()
        cut1_args['outdir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**cut1_args)
        cut1.run() # 
        cut1_fq1 = cut1.clean_fq1
        cut1_fq2 = cut1.clean_fq2
        n_total, n_cut1, _ = cut1.parse_log()
        copy_file(cut1.log, self.log_dir) # save log
        return (cut1_fq1, cut1_fq2, n_total, n_cut1)

    
    def rm_pcr_dup(self, fq1_in, fq2_in, fq1_out, fq2_out):
        log.info('rmdup() skipped, only works in single-end mode')
        symlink_file(fq1_in, fq1_out, absolute_path=False)
        symlink_file(fq2_in, fq2_out, absolute_path=False)

    
    def cut2(self, fq1_in, fq2_in, fq1_out, fq2_out):
        """
        --cut-after-trim
        """
        # sub-function
        # make sure
        # s: str
        # l: int, >= 0
        # r: int, <= 0
        def substr(s, l, r):
            s1 = s[l:]
            if r < 0:
                s1 = s1[:r]
            return s1
        # run
        n_cut2_rm = 0
        if not self.cut2_skip:
            seq_min = sum([self.len_min, self.cut2_l, abs(self.cut2_r)])
            try:
                with xopen(fq1_out, 'wt') as w1, xopen(fq2_out, 'wt') as w2:
                    for f1,f2 in zip(pyfastx.Fastx(fq1_in), pyfastx.Fastx(fq2_in)):
                        # name,seq,qual,comment
                        if len(f1[1]) < seq_min or len(f2[1]) < seq_min:
                            n_cut2_rm += 1
                            continue # drop short seq
                        # update name
                        n1 = '{} {}'.format(f1[0], f1[3]) if f1[3] else f1[0]
                        n2 = '{} {}'.format(f2[0], f2[3]) if f2[3] else f2[0]
                        # sub
                        s1 = substr(f1[1], self.cut2_l, self.cut2_r)
                        q1 = substr(f1[2], self.cut2_l, self.cut2_r)
                        s2 = substr(f2[1], self.cut2_l, self.cut2_r)
                        q2 = substr(f2[2], self.cut2_l, self.cut2_r)
                        # output
                        w1.write('@{}\n{}\n+\n{}\n'.format(n1, s1, q1))
                        w2.write('@{}\n{}\n+\n{}\n'.format(n2, s2, q2))
            except IOError as e:
                log.error(e)
        else:
            log.info('cut_after_trim(), skipped')
            symlink_file(fq1_in, fq1_out, absolute_path=False)
            symlink_file(fq2_in, fq2_out, absolute_path=False)
        return n_cut2_rm
    
    
    def run(self):
        ## re-run
        if all(file_exists([self.clean_fq1, self.clean_fq2])):
            log.info('Trim() skipped, file exists: {}, {}'.format(self.clean_fq1,
                                                                 self.clean_fq2))
            return None
        # 1. cut adapter
        cut1_fq1, cut1_fq2, n_total, n_cut1 = self.cutadapt()
        # 2. remove PCR dup
        self.rm_pcr_dup(cut1_fq1, cut1_fq2, self.rmdup_fq1, self.rmdup_fq2)
        n_rmdup = n_cut1 # rmdup skipped for PE reads
        # 3. cut, further adapters
        n_cut2_rm = self.cut2(self.rmdup_fq1, self.rmdup_fq2, 
                              self.cut2_fq1, self.cut2_fq2)
        n_cut2 = n_rmdup - n_cut2_rm
        # 4. wrap files
        if n_total < 1:
            n_total = 1
        n_cut1_rm = n_total - n_cut1
        n_rmdup_rm = n_cut1 - n_rmdup
        n_clean_pct = '{:.2f}'.format(n_cut2/n_total*100)
        # header
        header = ['#name', 'total', 'too_short', 'dup', 'too_short2', 
                  'clean', 'percent']
        s = [self.smp_name, n_total, n_cut1_rm, n_rmdup_rm, n_cut2_rm, n_cut2,
            n_clean_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_cut1_rm),
            'dup': int(n_rmdup_rm),
            'too_short2': int(n_cut2_rm),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        copy_file(self.cut2_fq1, self.clean_fq1)
        copy_file(self.cut2_fq2, self.clean_fq2)
        # 5. remove temp files
        del_list = [cut1_fq1, cut1_fq2, self.rmdup_fq1, self.rmdup_fq2,
                    self.cut2_fq1, self.cut2_fq2]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)


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
    parser.add_argument('-1', '--fq1', required=True,
        help='reads in FASTQ files, support (*.gz), 1 file.')
    parser.add_argument('-2', '--fq2', default=None,
        help='The read2 of pair-end reads')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')
    return parser


def get_args_trim(parser):
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
    parser.add_argument('--cut-after-trim', dest='cut_after_trim', default='0',
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
    parser.add_argument('-g', '--adapter5', default=None,
        help='5-Adapter, default: None')

    ## PE arguments
    parser.add_argument('-A', '--Adapter3', default=None,
        help='The 3 adapter of read2, default []')
    parser.add_argument('-G', '--Adapter5', default=None,
        help='The 5 adapter of read1, default: None')
    return parser


def main():
    args = vars(get_args().parse_args())
    TrimR1(**args).run()


if __name__ == '__main__':
    main()

#