#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: single fx

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
import re
import shutil
import pathlib
import argparse
import pandas as pd
import pyfastx
from xopen import xopen
from hiseq.utils.seq import Fastx
from hiseq.utils.file import check_path, file_exists, file_abspath, \
    check_fx, check_fx_paired, copy_file, symlink_file, remove_file, fx_name
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd
from hiseq.trim.cutadapt import Cutadapt



class TrimR1(object):
    """
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
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'library_type': None, # guess
            'fq1': None, # str
            'fq2': None, # str
            'outdir': None, # str
            'smp_name': None,
            'len_min': 15,
            'cut_to_length': 0,
            'threads': 4,
            'rmdup': False,
            'cut_after_trim': None, # skip, str
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
            'cut_before_trim': 0,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_r1' # 
        self.is_paired = file_exists(self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        self.init_fq()
        self.parse_cut2_args()
        Config().dump(self.__dict__.copy(), self.config_yaml)


    def init_fq(self): # se or pe
        c1 = isinstance(self.fq1, str)
        c1q = check_fx(self.fq1) # force
        if self.fq2 is None or self.fq2 == 'None':
            c2 = c2q = c2pe = False
            c2x = True
            self.fq2 = None
            self.rmdup_fq1, self.rmdup_fq2 = (self.rmdup_fq, None)
            self.cut2_fq1, self.cut2_fq2 = (self.cut2_fq, None)
            self.clean_fq1, self.clean_fq2 = (self.clean_fq, None)
        else:
            c2 = isinstance(self.fq2, str)
            c2q = check_fx(self.fq2) # force
            c2pe = check_fx_paired(self.fq1, self.fq2)
            c2x = all([c2, c2q, c2pe])
        # out
        c0 = all([c1, c1q, c2x])
        # show message
        msg = '\n'.join([
            '='*80,
            'Run TrimR1() with parameters:',
            '{:>14} : {}'.format('fq1', self.fq1),
            '{:>14} : {}'.format('fq2', self.fq2),
            '-'*40,
            '{:>14} : {}'.format('fq1 [str]', c1),
            '{:>14} : {}'.format('fq1 [fastq]', c1q),
            '{:>14} : {}'.format('fq2 [str]', c2),
            '{:>14} : {}'.format('fq2 [fastq]', c2q),
            '{:>14} : {}'.format('fq [paired]', c2pe),
            '{:>14} : {}'.format('rmdup', self.rmdup),
            '{:>14} : {}'.format('cut after trim', self.cut_after_trim),
            '-'*40,
            'Status: {}'.format(c0),
            '='*80,
        ])
        print(msg)
        if not c0:
            raise ValueError('fq1, fq2 no valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        
        
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
            self.log_dir])


    def run_rmdup(self, fq_in, fq_out):
        if self.rmdup:
            if file_exists(fq_out) and not self.overwrite:
                log.info('rmdup() skipped, file exists')
            elif file_exists(fq_in):
                Fastx(fq_in).collapse(fq_out, fq_out=True)
            else:
                log.error('run_rmdup() failed, file not exists: {}'.format(fq_in))
        else:
            symlink_file(fq_in, fq_out)


    def parse_cut2_args(self):
        """
        cut2: '10', '9,-9', '5', ...
        convert to 0-index position, python-style
        '10': [10:]
        '-9': [:-9]
        '9, -9': [9:-9]
        
        make sure: valid 
        
        return: start, end
        """
        if not isinstance(self.cut_after_trim, str):
            return None
        cut = self.cut_after_trim.strip().replace(' ', '') # trim ws
        p = re.compile('^(\d+)$|^(-\d+)$|^((\d+),(-\d+))$')
        m = p.search(cut) # match
        g = m.groups() # 1, 2, 3(4,5), 4, 5
        if g[0]:
            s = int(g[0])
            e = 0
        elif g[1]:
            s = 0
            e = int(g[1])
        elif g[2]:
            s = int(g[3])
            e = int(g[4])
        else:
            s = e = 0
            log.error('invalide --cut-after-trim, expect: 5; 5,-3; -3, \
                got {}'.format(cut))
            raise ValueError('check --cut-after-trim')
        self.cut2_s, self.cut2_e = (s, e)

    
    def run_cut2_se(self, fq_in, fq_out):
        """
        run cut_after_trim
          5    , [5:],   cut 5 nt from left
          -3   , [:-3],  cut 3 nt from right
          7,-8 , [7:-8], cut 7 nt from left, 8 nt from right
        """
        try:
            with xopen(fq_out, 'wt') as w:
                for name,seq,qual,comment in pyfastx.Fastx(fq_in):
                    if self.cut2_s > 0 and self.cut2_e < 0:
                        s = seq[self.cut2_s:self.cut2_e]
                        q = qual[self.cut2_s:self.cut2_e]
                    elif self.cut2_s > 0:
                        s = seq[self.cut2_s:]
                        q = qual[self.cut2_s:]
                    elif self.cut2_e < 0:
                        s = seq[:self.cut2_e]
                        q = qual[:self.cut2_e]
                    else:
                        continue
                    # min length
                    if len(s) < self.len_min:
                        continue
                    # update name
                    if comment:
                        name = '{} {}'.format(name, comment)
                    f = '@{}\n{}\n+\n{}'.format(name, s, q)
                    w.write(f+'\n')
        except IOError as e:
            log.error(e)


    def run_cut2_pe(self, fq1_in, fq2_in, fq1_out, fq2_out):
        """
        run cut_after_trim
          5    , [5:],   cut 5 nt from left
          -3   , [:-3],  cut 3 nt from right
          7,-8 , [7:-8], cut 7 nt from left, 8 nt from right
        """
        m = self.len_min + self.cut2_s + abs(self.cut2_e) # minimum length
        try:
            with xopen(fq1_out, 'wt') as w1, xopen(fq2_out, 'wt') as w2:
                for f1,f2 in zip(pyfastx.Fastx(fq2_in), pyfastx.Fastx(fq2_in)):
                    # minimum length
                    if len(f1[1]) < m:
                        continue
                    if self.cut2_s > 0 and self.cut2_e < 0:
                        s1 = f1[1][self.cut2_s:self.cut2_e]
                        q1 = f1[2][self.cut2_s:self.cut2_e]
                        s2 = f2[1][self.cut2_s:self.cut2_e]
                        q2 = f2[2][self.cut2_s:self.cut2_e]
                    elif self.cut2_s > 0:
                        s1 = f1[1][self.cut2_s:]
                        q1 = f1[1][self.cut2_s:]
                        s2 = f2[1][self.cut2_s:]
                        q2 = f2[1][self.cut2_s:]
                    elif self.cut2_e < 0:
                        s1 = f1[1][:self.cut2_e]
                        q1 = f1[1][:self.cut2_e]
                        s2 = f2[1][:self.cut2_e]
                        q2 = f2[1][:self.cut2_e]
                    else:
                        continue
                    # update name
                    if f1[-1]:
                        name1 = '{} {}'.format(f1[0], f1[-1])
                        name2 = '{} {}'.format(f2[0], f2[-1])
                    else:
                        name1 = f1[0]
                        name2 = f2[0]                        
                    # output
                    fq1 = '@{}\n{}\n+\n{}'.format(f1[1], s1, q1)
                    fq2 = '@{}\n{}\n+\n{}'.format(f2[1], s2, q2)
                    w1.write(fq1+'\n')
                    w2.write(fq2+'\n')
        except IOError as e:
            log.error(e)


    def run_se(self):
        """
        1. cut adapter
        2. remove duplicates
        3. cut after rmdup (CLIP, NSR)
        """
        args_local = self.__dict__.copy()
        # 1. cutadapt
        args_cutadapt = args_local.copy()
        args_cutadapt['outdir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**args_cutadapt)
        cut1.run() # 
        cut1_fq = cut1.clean_fq
        n_total, n_cut1, _ = cut1.parse_log()
        # 2. remove dup
        if self.rmdup:
            log.info('rmdup(), running')
            self.rm_dup(cut1_fq, self.rmdup_fq)
            n_rmdup = Fastx(self.rmdup_fq).number_of_seq()
        else:
            log.info('rmdup(), skipped')
            n_rmdup = n_cut1
            symlink_file(cut1_fq, self.rmdup_fq, absolute_path=False)
        # 3. cut after trim
        if self.cut_after_trim:
            log.info('cut_after_trim(), {}'.format(self.cut_after_trim))
            self.run_cut2_se(self.rmdup_fq, self.cut2_fq)
            n_cut2 = Fastx(self.cut2_fq).number_of_seq()
        else:
            log.info('cut_after_trim(), skipped')
            n_cut2 = n_rmdup
            symlink_file(self.rmdup_fq, self.cut2_fq)
        # 4. organize files
        # save: 1.cutadapt log/stat; 2. save nodup log/stat; 3. save cut2 log/stat
        # remove: 1. cutadapt fq; 2. nodup fq; 3. cut2() fq;
        # save stat
        if n_total < 1:
            n_total = 1
        n_rm_too_short = n_total - n_cut1
        n_rm_rmdup = n_cut1 - n_rmdup
        n_rm_cut2 = n_rmdup - n_cut2
        n_clean_pct = '{:.2f}'.format(n_cut2/n_total*100)
        # header
        header = ['#name', 'total', 'too_short', 'dup', 'too_short2', 'clean', 'percent']
        s = [self.smp_name, n_total, n_rm_too_short, n_rm_rmdup, n_rm_cut2, n_cut2, 
            n_clean_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_rm_too_short),
            'dup': int(n_rm_rmdup),
            'too_short2': int(n_rm_cut2),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # save fq files, cutadapt log
        copy_file(cut1.log, self.log_dir)
        if file_exists(self.cut2_fq) and not file_exists(self.clean_fq):
            copy_file(self.cut2_fq, self.clean_fq)
        # 5. remove temp files
        del_list = [cut1_fq, self.rmdup_fq, self.cut2_fq]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)


    def run_pe(self):
        """
        1. cut adapter
        2. remove duplicates
        3. cut after rmdup (CLIP, NSR)
        """
        args_local = self.__dict__.copy()
        # 1. cutadapt
        args_cutadapt = args_local.copy()
        args_cutadapt['outdir'] = self.cutadapt_dir # update
        cut1 = Cutadapt(**args_cutadapt)
        cut1.run()
        cut1_fq1 = cut1.clean_fq1
        cut1_fq2 = cut1.clean_fq2
        cut1_log = cut1.log
        n_total, n_cut1, _ = cut1.parse_log()
        # 2. remove dup
        if self.rmdup:
            log.info('rmdup(), running')
            log.info('remove dup for read1 only')
            self.rm_dup(cut1_fq1, self.rmdup_fq1)
            self.rm_dup(cut1_fq2, self.rmdup_fq2)
            n_rmdup = Fastx(self.rmdup_fq1).number_of_seq()
            # n_rmdup2 = Fastx(rmdup_fq2).number_of_seq()
        else:
            log.info('rmdup(), skipped')
            n_rmdup = n_cut1
            symlink_file(cut1_fq1, self.rmdup_fq1, absolute_path=False)
            symlink_file(cut1_fq2, self.rmdup_fq2, absolute_path=False)
        # 3. cut after trim
        if self.cut_after_trim:
            log.info('cut_after_trim(), {}'.format(self.cut_after_trim))
            self.run_cut2_pe(self.rmdup_fq1, self.rmdup_fq2, 
                self.cut2_fq1, self.cut2_fq2)
            n_cut2 = Fastx(self.cut2_fq1).number_of_seq()
        else:
            log.info('cut_after_trim(), skipped')
            n_cut2 = n_rmdup
            symlink_file(self.rmdup_fq1, self.cut2_fq1)
            symlink_file(self.rmdup_fq2, self.cut2_fq2)
        # 4. organize files
        # save: 1.cutadapt log/stat; 2. save nodup log/stat; 3. save cut2 log/stat
        # remove: 1. cutadapt fq; 2. nodup fq; 3. cut2() fq;
        # save stat
        
        n_rm_too_short = n_total - n_cut1
        n_rm_rmdup = n_cut1 - n_rmdup
        n_rm_cut2 = n_rmdup - n_cut2
        try:
            n_clean_pct = '{:.2f}'.format(n_cut2/n_total*100)
        except ZeroDivisionError as e:
            log.error(e)
            n_clean_pct = 0
        # header
        header = ['#name', 'total', 'too_short', 'dup', 'too_short2', 'clean', 'percent']
        s = [self.smp_name, n_total, n_rm_too_short, n_rm_rmdup, n_rm_cut2, n_cut2, 
            n_clean_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_rm_too_short),
            'dup': int(n_rm_rmdup),
            'too_short2': int(n_rm_cut2),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # save fq files, cutadapt log
        copy_file(cut1.log, self.log_dir)
        copy_file(self.cut2_fq1, self.clean_fq1)
        copy_file(self.cut2_fq2, self.clean_fq2)
        # 5. remove temp files
        del_list = [cut1_fq1, cut1_fq2, 
            self.rmdup_fq1, self.rmdup_fq2,
            self.cut2_fq1, self.cut2_fq2]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)
        # 6. output
        log.info('Output: {} of {} ({}%)'.format(n_cut2, n_total, n_clean_pct))


    def run(self):
        if self.fq2:
            if all(file_exists([self.clean_fq1, self.clean_fq2])) and not self.overwrite:
                log.info('trim() skipped, file exists: {}'.format(self.clean_fq1))
            else:
                self.run_pe()
        else:
            if file_exists(self.clean_fq) and not self.overwrite:
                log.info('trim() skipped, file exists: {}'.format(self.clean_fq))
            else:
                self.run_se()

                
                

def get_args():
    example = '\n'.join([
        'Examples:',
        '1. auto-detect adapters',
        '$ python cutadapt.py -1 fq1 -2 fq2 -o outdir -m 20',
        '2. specific 3-apdapter',
        '$ python cutadapt.py -1 fq1 -2 fq2 -o outdir -m 20 -a AGATCGGAAGAGCACACGTCTGAACT',
    ])
    parser = argparse.ArgumentParser(
        prog='cutadapt',
        description='cutadapt: for single PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-2', '--fq2', default=None,
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
    TrimR1(**args).run()


if __name__ == '__main__':
    main()

#