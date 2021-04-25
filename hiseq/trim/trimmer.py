
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
import logging
import pandas as pd
from multiprocessing import Pool
from hiseq.utils.args import args_init, ArgumentsInit, Adapter
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions
from Levenshtein import distance


def fq_paired(fq1, fq2):
    """
    Make sure fq1 and fq2, proper paired
    """
    flag = False
    if isinstance(fq1, str) and isinstance(fq2, str):
        flag = distance(fq1, fq2) == 1
    elif isinstance(fq1, list) and isinstance(fq2, list):
        if len(fq1) == len(fq2):
            flag = all([fq_paired(i, j) for i, j in zip(fq1, fq2)])

    return flag


def update_obj(obj, d, force=True, remove=False):
    """
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    # fresh start
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)

    return obj


class TrimRn(object):
    """
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
            'parallel_jobs': 1
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'hiseq_rn' # 

        # convert fq1, fq2 into list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]

        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]

        # check fq1
        if not isinstance(self.fq1, list):
            raise ValueError('--fq1, expect list, got {}'.format(type(self.fq1).__name__))

        if not all(file_exists(self.fq1)):
            raise ValueError('--fq1, file not exists')

        # check fq2
        if self.fq2 is None:
            pass
        elif isinstance(self.fq2, list):
            if not all(file_exists(self.fq2)):
                raise ValueError('--fq2, file not exists')

            # paired
            chk = [fq_paired(i, j) for i, j in zip(self.fq1, self.fq2)]
            if not all(chk):
                raise ValueError('--fq1, --fq2, not paired')
        else:
            raise ValueError('--fq2, expect list, got {}'.format(type(self.fq2).__name__))


    def pick_fq_samples(self, i):
        """
        Pick samples for each group
        """
        if isinstance(i, int):
            # fq1
            if isinstance(self.fq1, list):
                fq1 = self.fq1[i]
            else:
                fq1 = None

            # fq2
            if isinstance(self.fq2, list):
                fq2 = self.fq2[i]
            else:
                fq2 = None

            return (fq1, fq2)
        else:
            raise ValueError('int expected, {} found'.format(type(g).__name__))


    def run_fq_single(self, i):
        """
        Run for single fq

        i is the index of fq list
        """
        args_local = self.__dict__.copy()
        # update, fq1, fq2
        args_local['fq1'], args_local['fq2'] = self.pick_fq_samples(i)
        
        Trim(**args_local).run()


    def run(self):
        # run each fq
        if isinstance(self.fq1, list):
            i_list = list(range(len(self.fq1)))
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_fq_single, i_list)


class Trim(object):
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
        self.init_files()
        Toml(self.__dict__).to_toml(self.config_toml)
        # ## save arguments
        # args_checker(self.__dict__, self.config_pickle)
        # args_logger(self.__dict__, self.config_txt)


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
        self.hiseq_type = 'hiseq_r1' # 

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

        # output
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = fq_name(self.fq1, pe_fix=False)


    def init_files(self):
        """
        Init the files
        """
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        self.config_dir = self.project_dir + '/config'
        self.cutadapt_dir = self.project_dir + '/tmp/01_cutadapt'
        self.rmdup_dir = self.project_dir + '/tmp/02_rmdup'
        self.cut2_dir = self.project_dir + '/tmp/03_cut_after_trim'
        self.log_dir = self.project_dir + '/log'

        # files
        default_files = {
            'config_txt': self.config_dir + '/config.txt',
            'config_pickle': self.config_dir + '/config.pickle',
            'config_json': self.config_dir + '/config.json',
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
            'stat': self.project_dir + '/' + self.smp_name + '.trim.stat'
        }
        self = update_obj(self, default_files, force=True) # key

        check_path([
            self.config_dir, 
            self.cutadapt_dir, 
            self.rmdup_dir, 
            self.cut2_dir, 
            self.log_dir])


    def run_rmdup(self, fq_in, fq_out):
        """
        Remove duplicates, by sequence
        """
        if file_exists(fq_in):
            Fastx(fq_in).collapse(fq_out, fq_out=True)


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
        # format of region
        region = self.cut_to_region(cut_after_trim)

        if region:
            if file_exists(fq_in):
                Fastx(fq_in, len_min=self.len_min).subseq(fq_out, region=region)


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
            if all(file_exists([fq1_in, fq2_in])):
                Fastx(fq1_in, len_min=self.len_min).subseq_pe(fq2_in, fq1_out, fq2_out, region=region)


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
            file_symlink(cut1_fq, self.rmdup_fq, absolute_path=False)

        # 3. cut after trim
        if self.cut_after_trim:
            self.run_cut2(self.rmdup_fq, self.cut2_fq, self.cut_after_trim)
            n_cut2 = Fastx(self.cut2_fq).number_of_seq()
        else:
            n_cut2 = n_rmdup
            file_symlink(self.rmdup_fq, self.cut2_fq)

        # 4. organize files
        # save: 1.cutadapt log/stat; 2. save nodup log/stat; 3. save cut2 log/stat
        # remove: 1. cutadapt fq; 2. nodup fq; 3. cut2() fq;
        # save stat
        n_rm_too_short = n_total - n_cut1
        n_rm_rmdup = n_cut1 - n_rmdup
        n_rm_cut2 = n_rmdup - n_cut2
        with open(self.stat, 'wt') as w:
            w.write('\t'.join([
                self.smp_name, 
                str(n_total),
                str(n_rm_too_short),
                str(n_rm_rmdup),
                str(n_rm_cut2),
                str(n_cut2),
                '{:2f}'.format(n_cut2/n_total*100)
                ]) + '\n')

        # save fq files, cutadapt log
        file_copy(cut1.log, self.log_dir)
        if file_exists(self.cut2_fq) and not file_exists(self.clean_fq):
            file_copy(self.cut2_fq, self.clean_fq)

        # 5. remove temp files
        del_list = [cut1_fq, self.rmdup_fq, self.cut2_fq]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)


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
            file_symlink(cut1_fq1, self.rmdup_fq1, absolute_path=False)
            file_symlink(cut1_fq2, self.rmdup_fq2, absolute_path=False)

        # 3. cut after trim
        if self.cut_after_trim:
            self.run_cut2_pe(self.rmdup_fq1, self.rmdup_fq2, 
                self.cut2_fq1, self.cut2_fq2, self.cut_after_trim)
            n_cut2 = Fastx(self.cut2_fq1).number_of_seq()
        else:
            n_cut2 = n_rmdup
            file_symlink(self.rmdup_fq1, self.cut2_fq1)
            file_symlink(self.rmdup_fq2, self.cut2_fq2)

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

        with open(self.stat, 'wt') as w:
            w.write('\t'.join(header) + '\n')
            w.write(msg + '\n')

        # save fq files, cutadapt log
        file_copy(cut1.log, self.log_dir)
        file_copy(self.cut2_fq1, self.clean_fq1)
        file_copy(self.cut2_fq2, self.clean_fq2)

        # 5. remove temp files
        del_list = [cut1_fq1, cut1_fq2, 
            self.rmdup_fq1, self.rmdup_fq2,
            self.cut2_fq1, self.cut2_fq2]
        if not self.keep_tmp:
            file_remove(del_list, ask=False)

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
        Toml(self.__dict__).to_toml(self.config_toml)

        # # save config
        # chk0 = args_checker(self.__dict__, self.config_pickle)
        # chk1 = args_logger(self.__dict__, self.config_txt)


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
        self.hiseq_type = 'hiseq_r1' # 

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

        # output
        self.outdir = file_abspath(self.outdir)

        # smp_name
        self.smp_name = fq_name(self.fq1, pe_fix=False)


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
            'stat': project_prefix + '.cutadapt.stat'
        }
        self = update_obj(self, default_files, force=True) # key

        check_path(self.config_dir)


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

        # save to stat
        with open(self.stat, 'wt') as w:
            w.write('\t'.join([self.smp_name, str(n_total), str(n_clean), str(n_pct)]) + '\n')

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



@Logger('INFO')
class Trimmer(object):
    """Processing fastq files

    # for only 1 input fastq: (for smp_name)

    General
      a. low quality bases at 3'
      b. trim-n

    1. CLIP reads

      a. trim adapter (sliding)
      b. cut inline-barcode
      c. cut N-bases (5', 3' ends)
      e. collapse (remove duplicates)
      e. cut random barcode

    2. RNAseq reads

      a. trim adapter
      b. cut N-bases

    3. ChIPseq reads
      a. trim adapter
      b. cut N-bases (5', 3' ends)

    4. smRNAseq reads
      a. trim adapter
      b. discard untrimmed reads
      c. cut N-bases

    5. ATACseq reads
      a. trim adapter (sliding)

    """

    def __init__(self, **kwargs):
        """
        Main function
        """
        self.args = args_init(kwargs, trim=True)
        self.init_trimmer()


    def init_trimmer(self):
        # global vars
        self.fq1 = self.args.get('fq1', None) # checked in ArgumentsInit
        self.fq2 = self.args.get('fq2', None) # as above
        self.outdir = self.args.get('outdir', None) # as above
        assert isinstance(self.fq1, str)
        assert isinstance(self.outdir, str)
        assert is_path(self.outdir)

        # optional
        self.fqname = file_prefix(self.fq1)[0]
        if not self.fq2 is None:
            self.fqname = re.sub('[._][rR]?1$', '', self.fqname)
            self.fqname = re.sub('_1$', '', self.fqname)
        if self.args.get('smp_name', None):
            self.fqname = self.args['smp_name']

        self.fq_out_prefix = os.path.join(self.outdir, self.fqname) #self.args['fq_out_prefix']
        self.fq_out = self.fq_out_prefix + '.fq'    # SE
        self.fq_out1 = self.fq_out_prefix + '_1.fq' # PE1
        self.fq_out2 = self.fq_out_prefix + '_2.fq' # PE2

        # gzip
        if self.args.get('gzip', False):
            self.fq_out += '.gz'
            self.fq_out1 += '.gz'
            self.fq_out2 += '.gz'

        # files
        se_out = [self.fq_out, None]
        pe_out = [self.fq_out1, self.fq_out2]
        self.out_files = se_out if self.fq2 is None else pe_out


    def report(self):
        """
        Report the number of reads in each step
        """
        args = self.args.copy()
        rpt_file = self.fq_out_prefix + '.qc.stat'
        n_input = Fastx(self.fq1).count
        n_output = Fastx(self.out_files[0]).count # first one of output

        rpt_line = '#sample\tinput\toutput\tpercent\n'
        rpt_line += '{}\t{:d}\t{:d}\t{:.2f}%'.format(
            self.fqname,
            int(n_input),
            int(n_output),
            n_output / n_input * 100)

        with open(rpt_file, 'wt') as w:
            w.write(rpt_line + '\n')

        return [n_input, n_output]


    def wrap(self):
        """
        1. save log files: log/
        2. remove temp files: outdir/temp
        3. gzip output files:
        """
        args = self.args.copy()

        ## output
        log_dir = os.path.join(self.outdir, 'log')
        is_path(log_dir)

        log_cutadapt, log_rmdup, log_cut2 = [
            os.path.join(log_dir, self.fqname + i) for i in ['.cutadapt.log', '.rmdup.log', '.rmdup_cut.log']]

        ## 1. log file, cutadapt log        
        log_cutadapt_from = os.path.join(self.outdir, 'temp', '01.cutadapt_output', self.fqname + '.cutadapt.log')
        shutil.copy(log_cutadapt_from, log_cutadapt)

        ## 2. rmdup, cut
        with open(log_rmdup, 'wt') as w:
            w.write('remove_duplicate: {}'.format(args['rmdup']))

        ## 3. cut2
        with open(log_cut2, 'wt') as w:
            w.write('cut_after_trim: {}'.format(args['cut_after_trim']))

        ## 4. temp files
        temp_dir = os.path.join(self.outdir, 'temp')
        if not args['save_temp'] is True:
            shutil.rmtree(temp_dir, ignore_errors=True)


    def check(self):
        """
        Check arguments, target file exists, 
        """
        args = self.args.copy()

        ## save parameters
        args_pickle = self.fq_out_prefix + '.arguments.pickle'

        ## chk
        chk1 = args_checker(args, args_pickle)
        chk2 = args['overwrite'] is False

        if self.fq2 is None:
            chk3 = os.path.exists(self.out_files[0]) # SE
        else:
            chk3 = os.path.exists(self.out_files[0]) and os.path.exists(self.out_files[1]) # PE

        return all([chk1, chk2, chk3])


    def run(self):
        args = self.args.copy()

        ## update 'library-type'
        args = Adapter(libtype=args['library_type'], **args).trim_args

        if self.check():
            log.warning('{:>20} : file exists, skipped ...'.format(self.fqname))
            return self.out_files
        else:
            args_file = self.fq_out_prefix + '.arguments.txt'
            args_logger(args, args_file, True) # update arguments.txt

        ## 1.cutadapt directory
        args['outdir'] = os.path.join(self.outdir, 'temp', '01.cutadapt_output')
        f1, f2 = Cutadapt(**args).run()

        ## 2.rmdup
        rmdup_dir = os.path.join(self.outdir, 'temp', '02.rmdup_output')
        rmdup_f1 = os.path.join(rmdup_dir, os.path.basename(f1))
        rmdup_f2 = None if f2 is None else os.path.join(rmdup_dir, os.path.basename(f2))

        ## 3.cut-after-trim
        cut2_dir = os.path.join(self.outdir, 'temp', '03.cut_after_trim')
        cut2_f1 = os.path.join(rmdup_dir, 'cut.' + os.path.basename(f1))
        cut2_f2 = None if f2 is None else os.path.join(rmdup_dir, 'cut.' + os.path.basename(f2))

        ## create dir
        assert is_path(rmdup_dir)
        assert is_path(cut2_dir)

        # SE
        if self.fq2 is None:
            ## rmdup
            if args['rmdup']:
                Fastx(f1).collapse(rmdup_f1, fq_out=True)
            else:
                rmdup_f1 = f1

            ## cut2
            if eval(args['cut_after_trim']):
                args['cut'] = args.get('cut_after_trim', 0) # eg: 5,-3
                args['cut_to_length'] = args.get('cut_to_length', 0) # eg: 5,-3
                args['discard_tooshort'] = args.get('discard_tooshort', True)
                Fastx(rmdup_f1).cut(cut2_f1, **args)
            else:
                cut2_f1 = rmdup_f1

            ## clean
            if args['gzip']:
                gzip_cmd(cut2_f1, self.fq_out, decompress=False)
            else:
                shutil.move(cut2_f1, self.fq_out)
        # PE
        else:
            ## rmdup
            if args['rmdup']:
                Fastx(f1).collapse(rmdup_f1, fq_out=True)
                Fastx(f2).collapse(rmdup_f2, fq_out=True)
            else:
                rmdup_f1 = f1
                rmdup_f2 = f2

            # cut2
            if eval(args['cut_after_trim']): # default: '0'
                args['cut'] = args.get('cut_after_trim', 0) # eg: 5,-3
                args['cut_to_length'] = args.get('cut_to_length', 0) # eg: 5,-3
                args['discard_tooshort'] = args.get('discard_tooshort', True)
                if args['rmdup']:
                    ## if rmdup + cut2
                    Fastx(rmdup_f1).cut(cut2_f1, **args)
                    Fastx(rmdup_f2).cut(cut2_f2, **args)
                else:
                    ## if only cut2, require: PE pairs
                    Fastx(rmdup_f1).cut_pe(input2=rmdup_f2, out1=cut2_f1, out2=cut2_f2, **args)
            else:
                cut2_f1 = rmdup_f1
                cut2_f2 = rmdup_f2

            ## clean
            if args['gzip']:
                gzip_cmd(cut2_f1, self.fq_out1, decompress=False)
                gzip_cmd(cut2_f2, self.fq_out2, decompress=False)
            else:
                shutil.move(cut2_f1, self.fq_out1)
                shutil.move(cut2_f2, self.fq_out2)

        # wrap
        self.wrap()

        # report
        self.n_input, self.n_output = self.report()

        return self.out_files


# @Logger('INFO')
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





