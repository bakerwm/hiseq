
# -*- coding: utf-8 -*-

"""Quality control for fastq files

1. Trim N-bases (optional)
2. Trim adapter (guess adapter, see: trim_galore)
3. Trim N-bases (optional)

optional
1. --rm-untrim, --save-too-short, ...

"""

import os
import sys
import re
import shutil
import pathlib
import logging
import argparse
import pandas as pd
from multiprocessing import Pool
from Levenshtein import distance
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions
from hiseq.utils.file import check_file, check_path, file_exists, file_abspath, \
    check_fx, check_fx_paired, copy_file, symlink_file, remove_file, gzip_file
from hiseq.utils.utils import log, update_obj, Config
from hiseq.align.utils import check_fx_args


class Trim(object):
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
            'smp_name': None,
            'rmdup': False,
            'cut_after_trim': None,
            'cut_to_length': 0,
            'discard_tooshort': False,
            'overwrite': False,
            'parallel_jobs': 1
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_rn' #
        self.is_paired = self.fq2 is not None
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        self.fq1, self.fq2, self.outdir = file_abspath([self.fq1, self.fq2, self.outdir])
        # convert fq1, fq2 into list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2] # or None
        if not check_fx_args(self.fq1, self.fq2, check_empty=True):
            raise ValueError('fq1, fq2 not valid')
        # smp_name
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        elif isinstance(self.smp_name, list):
            if not len(self.fq1) == len(self.smp_name):
                raise ValueError('length of fq1 and smp_name not identical')

        
    def run_r1(self, f1):
        args_local = self.__dict__
        i = self.fq1.index(f1) # index
        args_local.update({
            'fq1': f1,
            'fq2': self.fq2[i] if self.is_paired else None,
            'outdir': self.outdir, # os.path.join(self.outdir, self.smp_name[i]),
            'smp_name': self.smp_name[i],
        })
        TrimR1(**args_local).run()


    def run(self):
        if len(self.fq1) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_r1, self.fq1)
        else:
            for f1 in self.fq1:
                self.run_r1(f1)


class TrimR1(object):
    """
    1. Trim adapter
    2. cut N-bases from ends
    3. remove dup
    2. trim n-bases from read

    2.1 TruSeq (NSR)  
        - cut 7, -7; from both ends of read (subseq)    

    2.2 TruSeq (iCLIP) 
        - cut 9; from read1    

    2.3 TruSeq (eCLIP)  
        - cut 10; -7 from read1  
        - cut 7, -10 from read2
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
            'len_min': 15,
            'threads': 4,
            'rmdup': False,
            'cut_after_trim': 0,
            'keep_tmp': False,
            'overwrite': False,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_r1' # 
        self.is_paired = file_exists(self.fq2)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
        self.init_files()
        Config().dump(self.__dict__, self.config_toml)


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
            'config_toml': self.config_dir + '/config.toml',
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
            'trim_toml': self.project_dir + '/' + self.smp_name + '.trim.toml',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path([
            self.config_dir, 
            self.cutadapt_dir, 
            self.rmdup_dir, 
            self.cut2_dir, 
            self.log_dir], create_dirs=True)


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


    def cut_to_region(self, cut):
        """
        Convert format: cut-after-trim, to region
        Cut n bases from read
        convert cut_after_trim to region

        cut_after_trim:
          5    , 6:-1 , cut 5 nt from left
          -3   , 1:-4 , cut 3 nt from right
          7,-8 , 8:-9 , cut 7 nt from left, 8 nt from right
        """
        cut = cut.strip().replace(' ', '') # remove white spaces
        p = re.compile('^(\d+)$|^(-\d+)$|^((\d+),(-\d+))$')
        m = p.search(cut) # match
        # 1, 2, 3, 4, 5
        if m.group(1):
            region = '{}:-1'.format(int(m.group(1))+1)
        elif m.group(2):
            region = '1:{}'.format(int(m.group(2))-1)
        elif m.group(3):
            region = '{}:{}'.format(int(m.group(4))+1, int(m.group(5))-1)
        else:
            region = None
            log.error('unknown cut format, expect 5|5:-3|-3, got {}'.format(cut))
        return region

    
    def run_cut2(self, fq_in, fq_out, cut_after_trim):
        """
        Cut n bases from read
        convert cut_after_trim to region

        cut_after_trim:
          5    , 6:-1 , cut 5 nt from left
          -3   , 1:-4 , cut 3 nt from right
          7,-8 , 8:-9 , cut 7 nt from left, 8 nt from right
        """
        region = self.cut_to_region(cut_after_trim)
        if region:
            if file_exists(fq_out) and not self.overwrite:
                log.info('run_cut2() skipped, file exists')
            elif file_exists(fq_in):
                Fastx(fq_in, len_min=self.len_min).subseq(fq_out, region=region)
            else:
                pass


    def run_cut2_pe(self, fq1_in, fq2_in, fq1_out, fq2_out, cut_after_trim):
        """
        Cut n bases from read
        convert cut_after_trim to region

        cut_after_trim:
    
          5    , 6:-1 , cut 5 nt from left
          -3   , 1:-4 , cut 3 nt from right
          7,-8 , 8:-9 , cut 7 nt from left, 8 nt from right

        """
        # format of region
        region = self.cut_to_region(cut_after_trim)
        if region:            
            if file_exists([fq1_out, fq2_out]) and not self.overwrite:
                log.info('run_cut2() skipped, file exists')
            elif all(file_exists([fq1_in, fq2_in])):
                Fastx(fq1_in, len_min=self.len_min).subseq_pe(
                    fq2_in, fq1_out, fq2_out, region=region)
            else:
                pass


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
            self.rm_dup(cut1_fq, self.rmdup_fq)
            n_rmdup = Fastx(self.rmdup_fq).number_of_seq()
        else:
            n_rmdup = n_cut1
            symlink_file(cut1_fq, self.rmdup_fq, absolute_path=False)
        # 3. cut after trim
        if self.cut_after_trim:
            self.run_cut2(self.rmdup_fq, self.cut2_fq, self.cut_after_trim)
            n_cut2 = Fastx(self.cut2_fq).number_of_seq()
        else:
            n_cut2 = n_rmdup
            symlink_file(self.rmdup_fq, self.cut2_fq)
        # 4. organize files
        # save: 1.cutadapt log/stat; 2. save nodup log/stat; 3. save cut2 log/stat
        # remove: 1. cutadapt fq; 2. nodup fq; 3. cut2() fq;
        # save stat
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
        # save to toml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_rm_too_short),
            'dup': int(n_rm_rmdup),
            'too_short2': int(n_rm_cut2),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_toml)
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
            log.info('remove dup for read1 only')
            self.rm_dup(cut1_fq1, self.rmdup_fq1)
            self.rm_dup(cut1_fq2, self.rmdup_fq2)
            n_rmdup = Fastx(self.rmdup_fq1).number_of_seq()
            # n_rmdup2 = Fastx(rmdup_fq2).number_of_seq()
        else:
            n_rmdup = n_cut1
            symlink_file(cut1_fq1, self.rmdup_fq1, absolute_path=False)
            symlink_file(cut1_fq2, self.rmdup_fq2, absolute_path=False)
        # 3. cut after trim
        if self.cut_after_trim:
            self.run_cut2_pe(self.rmdup_fq1, self.rmdup_fq2, 
                self.cut2_fq1, self.cut2_fq2, self.cut_after_trim)
            n_cut2 = Fastx(self.cut2_fq1).number_of_seq()
        else:
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
        n_clean_pct = '{:.2f}'.format(n_cut2/n_total*100)
        # header
        header = ['#name', 'total', 'too_short', 'dup', 'too_short2', 'clean', 'percent']
        s = [self.smp_name, n_total, n_rm_too_short, n_rm_rmdup, n_rm_cut2, n_cut2, 
            n_clean_pct]
        s = list(map(str, s))
        msg = '\t'.join(header) + '\n' + '\t'.join(s)
        with open(self.trim_stat, 'wt') as w:
            w.write(msg + '\n')
        # save to toml
        df_stat = {
            'name': self.smp_name,
            'total': int(n_total),
            'too_short': int(n_rm_too_short),
            'dup': int(n_rm_rmdup),
            'too_short2': int(n_rm_cut2),
            'clean': int(n_cut2),
            'percent': float(n_clean_pct),
        }
        Config().dump(df_stat, self.trim_toml)
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
        self.init_files()
        self.init_adapter()
        Config().dump(self.__dict__, self.config_toml)


    def init_args(self):
        args_init = {
            'library_type': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'len_min': 15,
            'cut_to_length': 0,
            'adapter3': None, # read1
            'adapter5': None, # read1
            'Adapter3': None, # read2
            'Adapter5': None, # read2
            'qual_min': 20,
            'error_rate': 0.1,
            'overlap': 3,
            'percent': 80,
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
            if not check_fx_paired(self.fq1, self.fq2):
                raise ValueError('--fq1, --fq2, file not paired')
        # output
        self.outdir = file_abspath(self.outdir)
        # smp_name
        self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)


    def init_adapter(self):
        """
        Determine the adapters
        1. library_type
        2. from adapter3, Adapter3
        3. guess from the first 1M reads
        """
        lib = {
                # 'truseq': 'AGATCGGAAGAGC',
                # 'nextera': 'CTGTCTCTTATACACATCT',
                # 'smallrna': 'TGGAATTCTCGG'
                'truseq': ['AGATCGGAAGAGCACACGT', 'AGATCGGAAGAGCGTCGTG'],
                'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'],
                'smallrna': ['TGGAATTCTCGGGTGCCAAGG', 'TGGAATTCTCGGGTGCCAAGG']
            }
        if isinstance(self.library_type, str):
            self.library_type = self.library_type.lower()
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
        elif isinstance(self.adapter3, str):
            pass
        else:
            log.info('Auto detect the adapters:')
            d = Fastx(self.fq1).detect_adapter()
            df = pd.DataFrame.from_dict(d, orient='index').sort_values(0, 
                ascending=False)
            self.library_type = df.index.to_list()[0] # first one
            self.adapter3, self.Adapter3 = lib.get(
                self.library_type, [None, None])
            ## save Auto-detect
            # df.to_csv(self.auto_adapter)


    def init_files(self):
        """
        Init the files
        """
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        project_prefix = self.project_dir + '/' + self.smp_name
        self.config_dir = self.project_dir + '/config'
        # files
        default_files = {
            # 'config_txt': self.config_dir + '/config.txt',
            # 'config_pickle': self.config_dir + '/config.pickle',
            # 'config_json': self.config_dir + '/config.json',
            'config_toml': self.config_dir + '/config.toml',
            'auto_adapter': project_prefix + '.auto_adapter.csv',
            'clean_fq': project_prefix + '.fq.gz',
            'clean_fq1': project_prefix + '_1.fq.gz',
            'clean_fq2': project_prefix + '_2.fq.gz',
            'untrim_fq': project_prefix + '.untrim.fq.gz',
            'untrim_fq1': project_prefix + '.untrim.1.fq.gz',
            'untrim_fq2': project_prefix + '.untrim.2.fq.gz',
            'too_short_fq': project_prefix + '.too_short.fq.gz',
            'too_short_fq1': project_prefix + '.too_short.1.fq.gz',
            'too_short_fq2': project_prefix + '.too_short.2.fq.gz',
            'too_long_fq': project_prefix + '.too_long.fq.gz',
            'too_long_fq1': project_prefix + '.too_long.1.fq.gz',
            'too_long_fq2': project_prefix + '.too_long.2.fq.gz',
            'log': project_prefix + '.cutadapt.log',
            'trim_stat': project_prefix + '.cutadapt.trim.stat',
            'trim_toml': project_prefix + '.cutadapt.trim.toml',
            'trim_json': project_prefix + '.cutadapt.trim.json',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path(self.config_dir, create_dirs=True)


    def get_ad3(self):
        """
        If recursive, 0, 1, 2, 3, 4
        """
        # adapter3
        if self.recursive:
            ad3_list = [self.adapter3[i:] for i in range(0, 10)]
            Ad3_list = [self.Adapter3[i:] for i in range(0, 10)]
        else:
            ad3_list = [self.adapter3]
            Ad3_list = [self.Adapter3]
        return (ad3_list, Ad3_list)


    def init_cmd_se(self):
        """
        cutadapt
        For single-end read:                               
        cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
        For paired-end reads:  
        cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
        """
        ## for SE
        self.arg_save_untrim = '--untrimmed-output={}'.format(self.untrim_fq) if self.save_untrim else ''
        self.arg_rm_untrim = '--discard-untrimmed' if self.rm_untrim else ''
        self.arg_too_short = '--too-short-output={}'.format(self.too_short_fq) if self.save_too_short else ''
        self.arg_too_long = '--too-long-output={}'.format(self.too_long_fq) if self.save_too_long else ''
        self.arg_cut_to_length = '--length {}'.format(self.cut_to_length) if self.cut_to_length else ''
        self.arg_adapter_5 = '-g {}'.format(self.adapter5) if self.adapter5 else ''
        self.arg_Adapter_5 = '-G {}'.format(self.Adapter5) if self.Adapter5 else ''
        ## adapter3
        ad3_list, _ = self.get_ad3()
        ad3_arg = ' '.join([
            '-a {}'.format(i) for i in ad3_list])
        cmd = ' '.join([
            '{}'.format(shutil.which('cutadapt')),
            '-j {}'.format(self.threads),
            ad3_arg,
            # '-a {}'.format(self.adapter3),
            '-j {}'.format(self.threads),
            '-m {}'.format(self.len_min),
            '-q {}'.format(self.qual_min),
            '-O {}'.format(self.overlap),
            r'-a "A{150}"',
            r'-a "C{150}"',
            r'-a "G{150}"',
            r'-a "T{150}"',
            '-e {}'.format(self.error_rate),
            '-n 3 --trim-n --max-n=0.1',
            self.arg_save_untrim,
            self.arg_rm_untrim,
            self.arg_too_short,
            self.arg_too_long,
            self.arg_cut_to_length,
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
        self.arg_save_untrim = '--untrimmed-output={}'.format(self.untrim_fq) if self.save_untrim else ''
        self.arg_rm_untrim = '--discard-untrimmed' if self.rm_untrim else ''
        self.arg_too_short = '--too-short-output={} --too-short-paired-output={}'.format(
            self.too_short_fq1, self.too_short_fq2) if self.save_too_short else ''
        self.arg_too_long = '--too-long-output={} --too-long-paired-output={}'.format(
            self.too_long_fq1, self.too_long_fq2) if self.save_too_long else ''
        self.arg_cut_to_length = '--length {}'.format(self.cut_to_length) if self.cut_to_length else ''
        self.arg_adapter_5 = '-g {}'.format(self.adapter5) if self.adapter5 else ''
        self.arg_Adapter_5 = '-G {}'.format(self.Adapter5) if self.Adapter5 else ''
        ## adapter3
        ad3_list, Ad3_list = self.get_ad3()
        ad3_arg = ' '.join([
            '-a {}'.format(i) for i in ad3_list])
        Ad3_arg = ' '.join([
            '-A {}'.format(i) for i in Ad3_list])
        cmd = ' '.join([
            '{}'.format(shutil.which('cutadapt')),
            '-j {}'.format(self.threads),
            # '-a {}'.format(self.adapter3),
            # '-A {}'.format(self.Adapter3),
            ad3_arg,
            Ad3_arg,
            self.arg_adapter_5,
            self.arg_Adapter_5,
            '-j {}'.format(self.threads),
            '-m {}'.format(self.len_min),
            '-q {}'.format(self.qual_min),
            '-O {}'.format(self.overlap),
            r'-a "A{150}"',
            r'-a "C{150}"',
            r'-a "G{150}"',
            r'-a "T{150}"',
            '-e {}'.format(self.error_rate),
            '-n 3 --trim-n --max-n=0.1',
            self.arg_save_untrim,
            self.arg_rm_untrim,
            self.arg_too_short,
            self.arg_too_long,
            self.arg_cut_to_length,
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
        # save to toml
        df_stat = {
            'name': self.smp_name,
            'total': n_total,
            'clean': n_clean,
            'percent': n_pct,
        }
        Config().dump(df_stat, self.trim_toml)
        Config().dump(df_stat, self.trim_json)
        # to global
        self.n_total = n_total
        self.n_clean = n_clean
        self.n_pct = n_pct
        return (n_total, n_clean, n_pct)


    def run(self):
        """
        SE or PE
        """
        cmd = self.init_cmd_pe() if self.fq2 else self.init_cmd_se()
        out_fq = self.clean_fq1 if self.fq2 else self.clean_fq
        ## run cmd
        cmd_txt = self.project_dir + '/cmd.sh'
        if file_exists(out_fq) and not self.overwrite:
            log.info('cutadapt() skipped, file exists: {}'.format(out_fq))
        else:
            with open(cmd_txt, 'wt') as w:
                w.write(cmd + '\n')
            try:
                run_shell_cmd(cmd)
                self.parse_log()
            except:
                log.error('cutadapt() failed, see: {}'.format(self.log))



def get_args():
    """
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser = argparse.ArgumentParser(
        description='hiseq qc, trim adapters and qc')
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
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')

    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1,
        type=int, help='Number of jobs to run in parallel, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    ## global arguments
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')

    ## specific
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
    parser.add_argument('--cut-after-trim', default='0',
        dest='cut_after_trim',
        help='cut n-bases after trimming adapter; positive value, \
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
    Trim(**args).run()


if __name__ == '__main__':
    main()
    
# 


## deprecated, re-write
# class TrimR1(object):
#     """Processing fastq files: 1 (single)

#     General
#       a. low quality bases at 3'
#       b. trim-n

#     1. CLIP reads
#       a. trim adapter (sliding)
#       b. cut inline-barcode
#       c. cut N-bases (5', 3' ends)
#       e. collapse (remove duplicates)
#       e. cut random barcode

#     2. RNAseq reads
#       a. trim adapter
#       b. cut N-bases

#     3. ChIPseq reads
#       a. trim adapter
#       b. cut N-bases (5', 3' ends)

#     4. smRNAseq reads
#       a. trim adapter
#       b. discard untrimmed reads
#       c. cut N-bases

#     5. ATACseq reads
#       a. trim adapter (sliding)
#     """
#     def __init__(self, **kwargs):
#         self = update_obj(self, kwargs, force=True)
#         self.init_args()

    
#     def init_args(self):
#         args_init = {
#             'fq1': None,
#             'fq2': None,
#             'outdir': None,
#             'smp_name': None,
#             'rmdup': False,
#             'cut_after_trim': None,
#             'cut_to_length': 0,
#             'discard_tooshort': False,
#             'overwrite': False,
#         }
#         self = update_obj(self, args_init, force=False)
#         self.hiseq_type = 'trim_r1'
#         self.is_paired = isinstance(self.fq2, str) # paired
#         if not check_fx_args(self.fq1, self.fq2, check_empty=True):
#             raise ValueError('fq1, fq2 not valid')
#         if self.outdir is None:
#             self.outdir = str(pathlib.Path.cwd())
#         self.outdir = file_abspath(os.path.expanduser(self.outdir))
#         if self.smp_name is None:
#             self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired)
#         # output
#         self.init_fx()
#         self.init_files()
    
    
#     def init_fx(self):
#         if self.is_paired:
#             self.fq_out1 = self.prefix + '_1.fq.gz'
#             self.fq_out2 = self.prefix + '_2.fq.gz'
#         else:
#             self.fq_out1 = self.prefix + '.fq.gz'
#             self.fq_out2 = None


#     def init_files(self):
#         self.config_toml = os.path.join(self.outdir, 'config.toml')
#         self.project_dir = os.path.join(self.outdir, self.smp_name)
#         self.prefix = os.path.join(self.project_dir, self.smp_name)
#         # subdirs
#         self.log_dir = os.path.join(self.project_dir, 'log')
#         self.temp_dir = os.path.join(self.project_dir, 'temp')
#         self.cutadapt_dir = os.path.join(self.temp_dir, '01.cutadapt_output')
#         self.rmdup_dir = os.path.join(self.temp_dir, '02.rmdup_output')
#         self.cut2_dir = os.path.join(self.temp_dir, '03.cut_after_trim')        
#         self.cutadapt_log = os.path.join(self.cutadapt_dir, self.smp_name+'.cutadapt.log')
#         self.rmdup_fq1 = os.path.join(self.rmdup_dir, os.path.basename(self.fq_out1))
#         self.cut2_fq1 = os.path.join(self.cut2_dir, os.path.basename(self.fq_out1))
#         self.log_cutadapt = os.path.join(self.smp_name + '.cutadapt.log')
#         self.log_rmdup = os.path.join(self.smp_name + '.rmdup.log')
#         self.log_cut2 = os.path.join(self.smp_name + '.rmdup_cut.log')
#         self.report_stat = self.prefix + '.qc.stat'
#         if self.is_paired:
#             self.rmdup_fq2 = os.path.join(self.rmdup_dir, os.path.basename(self.fq_out2))
#             self.cut2_fq2 = os.path.join(self.cut2_dir, os.path.basename(self.fq_out2))
#         else:
#             self.rmdup_fq2 = None
#             self.cut2_fq2 = None
        
    
#     def report(self):
#         n_input = Fastx(self.fq1).number_of_seq()
#         n_output = Fastx(self.fq_out1).number_of_seq()
#         head = '\t'.join(['#sample', 'input', 'output', 'percent'])
#         line = '\t'.join([
#             self.smp_name,
#             '{}'.format(n_input),
#             '{}'.format(n_output),
#             '{:.2f}%'.format(int(n_output)/int(n_input)*100)
#         ])
#         with open(self.report_stat, 'wt') as w:
#             w.write('\n'.join([head, line])+'\n')
#         return (n_input, n_output)


#     def run_cutadapt(self):
#         args_cutadapt = self.__dict__
#         args_cutadapt.update({
#             'outdir': self.cutadapt_dir
#         })
#         return Cutadapt(**args_cutadapt).run()
    
    
#     def run_rmdup(self, f1, f2):
#         if self.rmdup:
#             if file_exists(self.rmdup_fq1) and not self.overwrite:
#                 pass
#             else:
#                 Fastx(f1).collapse(self.rmdup_fq1, fq_out=True)
#                 if self.is_paired and file_exists(f2):
#                     Fastx(f2).collapse(self.rmdup_fq2, fq_out=True)
#         else:
#             symlink_file(f1, self.rmdup_fq1)
#             if self.is_paired and file_exists(f2):
#                 symlink_file(f2, self.rmdup_fq2)
        
        
#     def run_cut2(self):
#         if isinstance(self.cut_after_trim, str):
#             if file_exists(self.cut2_fq1) and not self.overwrite:
#                 pass
#             else:
#                 args_cut2 = {
#                     'cut': self.cut_after_trim,
#                     'cut_to_length': self.cut_to_length,
#                     'discard_tooshort': self.discard_tooshort
#                 }
#                 Fastx(self.rmdup_fq1).cut(self.cut2_fq1, **args_cut2)
#                 if self.is_paired and file_exists(self.rmdup_fq2):
#                     Fastx(self.rmdup_fq2).cut(self.cut2_fq2, **args_cut2)
#         else:
#             symlink_file(self.rmdup_fq1, self.cut2_fq1)
#             if self.is_paired and file_exists(self.rmdup_fq2):
#                 symlink_file(self.rmdup_fq2, self.cut2_fq2)
                
    
#     def run(self):
#         # to-be-create-dirs
#         dir_list = [
#             self.log_dir, self.temp_dir, self.cutadapt_dir, 
#             self.rmdup_dir, self.cut2_dir
#         ]
#         check_path(dir_list, create_dirs=True)
#         Config().dump(self.__dict__, self.config_toml) # save config
#         f1, f2 = self.run_cutadapt()
#         self.run_rmdup(f1, f2)
#         self.run_cut2()
#         self.wrap()
#         self.n_input, self.n_output = self.report()
#         return self.out_files
          
    
## deprecated, replaced by: hiseq.utils.file.check_fx_paired()
# def fq_paired(fq1, fq2):
#     """
#     Make sure fq1 and fq2, proper paired
#     """
#     flag = False
#     if isinstance(fq1, str) and isinstance(fq2, str):
#         flag = distance(fq1, fq2) == 1
#     elif isinstance(fq1, list) and isinstance(fq2, list):
#         if len(fq1) == len(fq2):
#             flag = all([fq_paired(i, j) for i, j in zip(fq1, fq2)])
#     return flag


## deprecated, replaced by: hiseq.utils.utils.update_obj()
# def update_obj(obj, d, force=True, remove=False):
#     """
#     d: dict
#     force: bool, update exists attributes
#     remove: bool, remove exists attributes
#     Update attributes from dict
#     force exists attr
#     """
#     # fresh start
#     if remove is True:
#         for k in obj.__dict__:
#             delattr(obj, k)
#     # add attributes
#     if isinstance(d, dict):
#         for k, v in d.items():
#             if not hasattr(obj, k) or force:
#                 setattr(obj, k, v)

#     return obj

    
    
    

# class Cutadapt(object):

#     # def __init__(self, fq, outdir, fq2=None, **kwargs):
#     def __init__(self, **kwargs):
#         """
#         Trim adapter and low quality bases for fastq file(s) using ``cutadapt``
#         program version 2.7. the program will guess the adapters: TruSeq,
#         Nextera, small RNAseq, ...


#         *fq* is typically the name of fastq file

#         *outdir* is the path to the directory saving the results

#         *fq2* is the second read file of PE sequencing, optinoal

#         Usage is to specify the fastq file, adapter and other arguments:

#             Cutadapt('read1.fq', 'output', 'read2.fq')

#         """
#         # prepare args
#         self.args = ArgumentsInit(kwargs, trim=True).dict.__dict__
#         self.args.pop('args_input', None)
#         self.args.pop('cmd_input', None)
#         self.args.pop('dict', None)
#         self.init_cut()


#     def init_cut(self):
#         """
#         Check arguments: fq, fq2, outdir
#         """
#         args = self.args.copy()

#         # global vars
#         self.fq1 = self.fq = args.get('fq1', None) # checked in ArgumentsInit
#         self.fq2 = args.get('fq2', None) # as above
#         self.outdir = args.get('outdir', None) # as above

#         assert isinstance(self.fq1, str)
#         assert isinstance(self.outdir, str)
#         assert is_path(self.outdir)

#         # optional
#         self.fqname = file_prefix(self.fq1)[0]
#         if not self.fq2 is None:
#             self.fqname = re.sub('[._][rR]?1$', '', self.fqname)
#             self.fqname = re.sub('_1$', '', self.fqname)
#         if args.get('smp_name', None):
#             self.fqname = args['smp_name']

#         self.fq_out_prefix = os.path.join(self.outdir, self.fqname) #self.args['fq_out_prefix']
#         self.fq_out = self.fq_out_prefix + '.fq'    # SE
#         self.fq_out1 = self.fq_out_prefix + '_1.fq' # PE1
#         self.fq_out2 = self.fq_out_prefix + '_2.fq' # PE2
#         self.log  = self.fq_out_prefix + '.cutadapt.log'
#         self.cutadapt = shutil.which('cutadapt') # in PATH

#         # files
#         se_out = [self.fq_out, None]
#         pe_out = [self.fq_out1, self.fq_out2]
#         self.out_files = se_out if self.fq2 is None else pe_out


#     def cut_cut(self, cutadapt_arg=True, read2=False):
#         """
#         Number of bases to cut from left or right of sequence
#         5: from left
#         -3: from right
#         for cutadapt command:
#         --cut={cut}
#         """
#         args = self.args.copy()

#         cut_before_trim = args['cut_before_trim']
#         cut_opts = cut_before_trim.split(',') if ',' in cut_before_trim else [cut_before_trim]

#         if cutadapt_arg:
#             if read2:
#                 p = '-U' # read2
#             else:
#                 p = '-u' # read1
#             out = ' '.join(['{} {}'.format(p, i) for i in cut_opts])
#         else:
#             out = cut_opts

#         return out


#     def adapter_sliding(self, adapter, step=2, window=15):
#         """
#         For some situation, the adapter (3') in library was differ, especially
#         for ligation strategy libraries, (eg: iCLIP, eCLIP).
#         We need to make sliding windows of the adatpers for adapter removing
#         """
#         args = self.args.copy()

#         # ## always, the 3' adapter
#         # adapter = args['adapter3']

#         ## sliding
#         print(adapter, window, step)
#         adapter_sliders = [adapter[i:i+window] for i in range(0, len(adapter)-window, step) if i >= 0]
#         if not adapter_sliders:
#             adapter_sliders = [adapter] # full length

#         return adapter_sliders


#     def get_cmd(self):
#         """
#         Create basic command line for cutadapt program: SE read

#         -a <ad> -m 15 -q 20 --trim-n --max-n=0.1 --error-rate=0.1

#         """
#         args = self.args.copy()

#         # 3' adapter
#         ad3_list = self.adapter_sliding(args['adapter3']) if args['adapter_sliding'] else [args['adapter3']]
#         arg_ad3 = ' '.join(['-a {}'.format(i) for i in ad3_list])

#         # SE mode
#         if self.fq2 is None:
#             arg_AD3 = '' # SE
#             arg_out = '-o {}'.format(self.fq_out)
#             # cut
#             if args['cut_before_trim']:
#                 arg_cut = self.cut_cut()

#             # untrimmed
#             if args['save_untrim']:
#                 arg_untrim = '--untrimmed-output={}'.format(self.fq_out_prefix + '.untrimmed.fq')
#                 args['threads'] = 1
#             elif args['rm_untrim']:
#                 arg_untrim = '--discard-untrimmed'
#             else:
#                 arg_untrim = '' # save untrimmed

#             # too-short
#             if args['save_too_short']:
#                 arg_short = '--too-short-output={}'.format(self.fq_out_prefix + '.too_short.fq')
#             else:
#                 arg_short = ''

#             # too-long
#             if args['save_too_long']:
#                 arg_long = '--too-long-output={}'.format(self.fq_out_prefix + '.too_long.fq')
#             else:
#                 arg_long = ''

#             # cut to length
#             if args['cut_to_length'] >= args['len_min']:
#                 arg_cut_to_length = '--length={}'.format(args['cut_to_length'])
#             else:
#                 arg_cut_to_length = ''

#         else:
#             AD3_list = self.adapter_sliding(args['AD3']) if args['adapter_sliding'] else [args['AD3']]
#             arg_AD3 = ' '.join(['-A {}'.format(i) for i in AD3_list])
#             arg_out = '-o {} -p {}'.format(self.fq_out1, self.fq_out2)

#             # cut
#             if args['cut_before_trim']:
#                 arg_cut1 = self.cut_cut()
#                 arg_cut2 = self.cut_cut(read2=True)
#                 arg_cut = arg_cut1 + ' ' + arg_cut2

#             # untrimmed
#             if args['save_untrim']:
#                 arg_untrim = '--untrimmed-output={} \
#                     --untrimmed-paired-output={}'.format(
#                         self.fq_out_prefix + '.untrimmed_1.fq',
#                         self.fq_out_prefix + '.untrimmed_2.fq')
#                 args['threads'] = 1
#             elif args['rm_untrim']:
#                 arg_untrim = '--discard-untrimmed'
#             else:
#                 arg_untrim = '' # save untrimmed

#             # too-short
#             if args['save_too_short']:
#                 arg_short = '--too-short-output={} \
#                     --too-short-paired-output={}'.format(
#                         self.fq_out_prefix + '.too_short_1.fq',
#                         self.fq_out_prefix + '.too_short_2.fq')
#             else:
#                 arg_short = ''

#             # too-long
#             if args['save_too_long']:
#                 arg_long = '--too-long-output={} \
#                     --too-long-paired-output={}'.format(
#                         self.fq_out_prefix + '.too_long_1.fq',
#                         self.fq_out_prefix + '.too_long_2.fq')
#             else:
#                 arg_long = ''

#             # cut to length
#             if args['cut_to_length'] >= args['len_min']:
#                 arg_cut_to_length = '--length={}'.format(args['cut_to_length'])
#             else:
#                 arg_cut_to_length = ''

#         ## save log
#         arg_log = ''
#         ## command line

#         # command line
#         arg_main = '{} -m {} -q {} --trim-n --max-n=0.1 --error-rate={} \
#             --times={} --cores {}'.format(
#                 self.cutadapt,
#                 args['len_min'],
#                 args['qual_min'],
#                 args['error_rate'],
#                 args['trim_times'],
#                 args['threads'])

#         cmd_line = ' '.join([
#             arg_main,
#             arg_ad3,
#             arg_AD3,
#             arg_cut,
#             arg_untrim,
#             arg_short,
#             arg_long,
#             arg_cut_to_length])

#         return cmd_line


#     def run_se(self):
#         args = self.args.copy()

#         cmd_line = self.get_cmd()

#         if os.path.exists(self.fq_out) and args['overwrite'] is False:
#             log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
#         else:
#             cmd_line += ' -o {} {} 1>{}'.format(
#                 self.fq_out,
#                 self.fq,
#                 self.log)
#             run_shell_cmd(cmd_line)

#         return [self.fq_out, None]


#     def run_pe(self):
#         args = self.args.copy()

#         cmd_line = self.get_cmd()

#         if os.path.exists(self.fq_out1) and os.path.exists(self.fq_out2) and args['overwrite'] is False:
#             log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
#         else:
#             cmd_line += ' -o {} -p {} {} {} 1>{}'.format(
#                 self.fq_out1,
#                 self.fq_out2,
#                 self.fq1,
#                 self.fq2,
#                 self.log)
#             run_shell_cmd(cmd_line)

#         return [self.fq_out1, self.fq_out2]


#     def run(self):
#         args = self.args.copy()

#         # SE
#         if self.fq2 is None:
#             fq_out = self.run_se()
#         # PE
#         else:
#             fq_out = self.run_pe()

#         return fq_out





