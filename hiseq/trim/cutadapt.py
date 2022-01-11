#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrap cutadapt commands
"""


import os
import re
import shutil
import pathlib
import argparse
import pandas as pd
from hiseq.utils.seq import Fastx
from hiseq.utils.file import (
    check_path, file_exists, file_abspath, check_fx, check_fx_paired, fx_name
)
from hiseq.utils.utils import log, update_obj, Config, run_shell_cmd


class Cutadapt(object):
    """
    single file (pair)
    Trimmer

    for cutadapt:
    1. adapters; quality; min_len

    1. trim adapters

    1.1 TruSeq
        - read1: AGATCGGAAGAGCACACGTCTGAACT 
        - read2: AGATCGGAAGAGCGTCGTGTAGGGAA
        - common: AGATCGGAAGAGC

    1.2 TruSeq (NSR) 
        - read1: CT..TC AGATCGGAAGAGCACACGTCTGAACT 
        - read2: GA..AG AGATCGGAAGAGCGTCGTGTAGGGAA
        trim 7, -7; from both ends

    1.3 Nextera
        - read1: CTGTCTCTTATACACATCT 
        - read2: CTGTCTCTTATACACATCT

    1.4 small RNA  
        - read1: TGGAATTCTCGGGTGCCAAGG
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
            'smp_name': None,
            'len_min': 15,
            'len_max': 0,
            'cut_to_length': 0,
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
            'threads': 4,
            'overwrite': False,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cutadapt_r1' #
        self.is_paired = file_exists(self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        self.init_fq()
        self.init_adapter()
        Config().dump(self.__dict__.copy(), self.config_yaml)
        
        
    def init_fq(self): # se or pe
        self.fq1 = file_abspath(self.fq1)
        c1 = isinstance(self.fq1, str)
        c1q = check_fx(self.fq1) # force
        if self.fq2 is None or self.fq2 == 'None':
            # update for SE
            self.fq2 = None
            self.clean_fq1, self.clean_fq2 = (self.clean_fq, None)
            self.untrim_fq1, self.untrim_fq2 = (self.untrim_fq, None)
            self.too_short_fq1, self.too_short_fq2 = (self.too_short_fq, None)
            self.too_long_fq1, self.too_long_fq2 = (self.too_long_fq, None)
            self.clean_fq2 = None
            c2 = c2q = c2pe = False
            c2x = True
        else:
            self.fq2 = file_abspath(self.fq2)
            c2 = isinstance(self.fq2, str)
            c2q = check_fx(self.fq2) # force
            c2pe = check_fx_paired(self.fq1, self.fq2)
            c2x = all([c2, c2q, c2pe])
        # out
        c0 = all([c1, c1q, c2x])
        # show message
        msg = '\n'.join([
            '='*80,
            'Run Cutadapt() with parameters:',
            '{:>14} : {}'.format('fq1', self.fq1),
            '{:>14} : {}'.format('fq2', self.fq2),
            '-'*40,
            '{:>14} : {}'.format('fq1 [str]', c1),
            '{:>14} : {}'.format('fq1 [fastq]', c1q),
            '{:>14} : {}'.format('fq2 [str]', c2),
            '{:>14} : {}'.format('fq2 [fastq]', c2q),
            '{:>14} : {}'.format('fq [paired]', c2pe),
            '-'*40,
            'Status: {}'.format(c0),
            '='*80,
        ])
        print(msg)
        if not c0:
            raise ValueError('fq1, fq2 no valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)


    def init_adapter(self):
        """
        Guess the adapters
        1. library_type
        2. from adapter3, Adapter3
        3. guess from the first 1M reads
        """
        lib = {
                'truseq': ['AGATCGGAAGAGC', 'AGATCGGAAGAGC'],
                'nextera': ['CTGTCTCTTATACA', 'CTGTCTCTTATACA'],
                'smallrna': ['TGGAATTCTCGGGTGCCAAGG', 'TGGAATTCTCGGGTGCCAAGG']
            }
        if self.library_type is None and self.adapter3 is None:
            log.info('Auto detect the adapters:')
            try:
                d = Fastx(self.fq1).detect_adapter()
                df = pd.DataFrame.from_dict(d, orient='index').sort_values(0, 
                    ascending=False)
                self.library_type = df.index.to_list()[0] # first one
            except KeyError as e:
                log.error(e)
                self.library_type = 'truseq'
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
        elif isinstance(self.library_type, str):
            self.library_type = self.library_type.lower()
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
        else:
            pass
            # from args


    def init_files(self):
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        self.config_dir = self.project_dir + '/config'
        prefix = self.project_dir + '/' + self.smp_name
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'auto_adapter': prefix + '.auto_adapter.csv',
            'clean_fq': prefix + '.fq.gz', # update for SE
            'clean_fq1': prefix + '_1.fq.gz',
            'clean_fq2': prefix + '_2.fq.gz', 
            'untrim_fq': prefix + '.untrim.fq.gz', # update for SE
            'untrim_fq1': prefix + '.untrim.1.fq.gz',
            'untrim_fq2': prefix + '.untrim.2.fq.gz',
            'too_short_fq': prefix + '.too_short.fq.gz', # update for SE
            'too_short_fq1': prefix + '.too_short.1.fq.gz',
            'too_short_fq2': prefix + '.too_short.2.fq.gz',
            'too_long_fq': prefix + '.too_long.fq.gz', # update for SE
            'too_long_fq1': prefix + '.too_long.1.fq.gz',
            'too_long_fq2': prefix + '.too_long.2.fq.gz',
            'log': prefix + '.cutadapt.log',
            'trim_stat': prefix + '.cutadapt.trim.stat',
            'trim_yaml': prefix + '.cutadapt.trim.yaml',
            'trim_json': prefix + '.cutadapt.trim.json',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path(self.config_dir, create_dirs=True)


    def get_ad3(self):
        """
        If recursive, 0, 1, 2, 3, 4
        """
        # adapter3
        if self.recursive:
            if isinstance(self.adapter3, str):
                ad3_list = [self.adapter3[i:i+14] for i in range(0, 10)]
            else:
                ad3_list = []
            if isinstance(self.Adapter3, str):
                Ad3_list = [self.Adapter3[i:i+14] for i in range(0, 10)]
            else:
                Ad3_list = []
        else:
            if isinstance(self.adapter3, str):
                ad3_list = [self.adapter3]
            else:
                ad3_lsit = []
            if isinstance(self.Adapter3, str):
                Ad3_list = [self.Adapter3]
            else:
                Ad3_list = []
        return (ad3_list, Ad3_list)


    def init_cmd_se(self):
        """
        cutadapt
        For single-end read:                               
        cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
        For paired-end reads:  
        cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
        """
        # for SE, specific args
        arg_untrim = '--untrimmed-output={}'.format(self.untrim_fq) if self.save_untrim else ''
        arg_rm_untrim = '--discard-untrimmed' if self.rm_untrim else ''
        arg_too_short = '--too-short-output={}'.format(self.too_short_fq) if self.save_too_short else ''
        arg_too_long = '--too-long-output={}'.format(self.too_long_fq) if self.save_too_long else ''
        arg_cut_to_length = '--length {}'.format(self.cut_to_length) if self.cut_to_length else ''
        arg_adapter_5 = '-g {}'.format(self.adapter5) if self.adapter5 else ''
        arg_Adapter_5 = '-G {}'.format(self.Adapter5) if self.Adapter5 else ''
        arg_rm_polyN = r'-a "A{50}" -a "C{50}" -a "G{50}" -a "T{50}"' if self.rm_polyN else ''
        arg_cut_before_trim = '--cut {}'.format(self.cut_before_trim) if self.cut_before_trim else ''
        # adapter3
        ad3_list, _ = self.get_ad3()
        ad3_arg = ' '.join([
            '-a {}'.format(i) for i in ad3_list])
        arg_len_max = '-M {}'.format(self.len_max) if self.len_max > 0 else ''
        cmd = ' '.join([
            '{}'.format(shutil.which('cutadapt')),
            ad3_arg,
            arg_adapter_5,
            arg_len_max,
            '-j {}'.format(self.threads),
            '-m {}'.format(self.len_min),
            '-q {}'.format(self.qual_min),
            '-O {}'.format(self.overlap),
            '-e {}'.format(self.error_rate),
            '-n 3 --trim-n --max-n=0.1',
            arg_rm_polyN,
            arg_untrim,
            arg_rm_untrim,
            arg_too_short,
            arg_too_long,
            arg_cut_to_length,
            arg_cut_before_trim,
            '-o {}'.format(self.clean_fq),
            '{}'.format(self.fq1),
            '1>{}'.format(self.log)
            ])
        return cmd


    def init_cmd_pe(self):
        """
        cutadapt
        For single-end read:                               
        cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
        For paired-end reads:  
        cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
        """
        arg_untrim = '--untrimmed-output={}'.format(self.untrim_fq) if self.save_untrim else ''
        arg_rm_untrim = '--discard-untrimmed' if self.rm_untrim else ''
        arg_cut_to_length = '--length {}'.format(self.cut_to_length) if self.cut_to_length else ''
        arg_adapter_5 = '-g {}'.format(self.adapter5) if self.adapter5 else ''
        arg_Adapter_5 = '-G {}'.format(self.Adapter5) if self.Adapter5 else ''
        arg_rm_polyN = r'-a "A{50}" -a "C{50}" -a "G{50}" -a "T{50}"' if self.rm_polyN else ''
        if self.save_too_short:
            arg_too_short = '--too-short-output={} --too-short-paired-output={}'.format(
                self.too_short_fq1, self.too_short_fq2)
        else:
            arg_too_short = ''
        if self.save_too_long:
            arg_too_long = '--too-long-output={} --too-long-paired-output={}'.format(
                self.too_long_fq1, self.too_long_fq2)
        else:
            arg_too_long = ''
        ## adapter3
        ad3_list, Ad3_list = self.get_ad3()
        ad3_arg = ' '.join([
            '-a {}'.format(i) for i in ad3_list])
        Ad3_arg = ' '.join([
            '-A {}'.format(i) for i in Ad3_list])
        cmd = ' '.join([
            '{}'.format(shutil.which('cutadapt')),
            ad3_arg,
            Ad3_arg,
            arg_adapter_5,
            arg_Adapter_5,
            '-j {}'.format(self.threads),
            '-m {}'.format(self.len_min),
            '-q {}'.format(self.qual_min),
            '-O {}'.format(self.overlap),
            '-e {}'.format(self.error_rate),
            '-n 3 --trim-n --max-n=0.1',
            arg_rm_polyN,
            arg_untrim,
            arg_rm_untrim,
            arg_too_short,
            arg_too_long,
            arg_cut_to_length,
            '-o {}'.format(self.clean_fq1),
            '-p {}'.format(self.clean_fq2),
            '{} {}'.format(self.fq1, self.fq2),
            '1>{}'.format(self.log)
            ])
        return cmd


    def parse_log(self):
        """
        Parse the log file: .cutadapt.log

        SE:

        Total reads processed:               1,000,000
        Reads with adapters:                    78,522 (7.9%)
        Reads written (passing filters):       997,922 (99.8%)

        PE:

        Total read pairs processed:          1,000,000
          Read 1 with adapter:                  78,522 (7.9%)
          Read 2 with adapter:                  43,801 (4.4%)
        Pairs written (passing filters):       995,338 (99.5%)
        """
        p = re.compile('processed:\s+([\d,]+).*\n.*passing filters\):\s+([\d,]+)', re.DOTALL)
        n_total = 0
        n_clean = 0
        n_pct = 0
        if file_exists(self.log):
            with open(self.log, 'rt') as r:
                log_txt = r.read()
            m = p.search(log_txt) # match
            if m:
                n_total = int(m.group(1).replace(',', ''))
                n_clean = int(m.group(2).replace(',', ''))
                n_pct = '{:.2f}'.format(n_clean/n_total*100)
        header = ['#name', 'total', 'clean', 'percent']
        s = [self.smp_name, n_total, n_clean, n_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to yaml
        df_stat = {
            'name': self.smp_name,
            'total': n_total,
            'clean': n_clean,
            'percent': n_pct,
        }
        Config().dump(df_stat, self.trim_yaml)
        Config().dump(df_stat, self.trim_json)
        # to global
        self.n_total = n_total
        self.n_clean = n_clean
        self.n_pct = n_pct
        return (n_total, n_clean, n_pct)


    def run(self):
        cmd = self.init_cmd_pe() if self.is_paired else self.init_cmd_se()
        out_fq = self.clean_fq1 if self.is_paired else self.clean_fq
        ## run cmd
        cmd_txt = self.project_dir + '/cmd.sh'
        if file_exists(out_fq) and not self.overwrite:
            log.info('Cutadapt() skipped, file exists: {}'.format(out_fq))
        else:
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')
            try:
                run_shell_cmd(cmd)
                self.parse_log()
            except:
                log.error('Cutadapt() failed, see: {}'.format(self.log))


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
    parser.add_argument('-m', '--len-min', dest='len_min', default=15,
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-M', '--len-max', dest='len_max', default=0,
        type=int, help='Maxmimum length of reads after trimming, defualt [0], ignore')
    parser.add_argument('--cut-to-length', default=0, dest='cut_to_length',
        type=int,
        help='cut reads to from right, default: [0], full length')
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
    Cutadapt(**args).run()


if __name__ == '__main__':
    main()

#