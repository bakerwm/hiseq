#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
=========
Attention
=========
This script was designed for small RNA seqeuencing reads, by paired-end mode with 150bp read length.
expect read through the insert, and including both 5'UMI, barcode, 3'adapter in reads


Trimming reads for piRNA seq data

1. smRNA library structure (read1)
{5'-UMI}-{insert}-{3'inline-barcode}-{3'adapter}

=========
Sequences
5'-UMI:
  - 15nt
  - at the left most of read1
  - in read1: NNNCGANNNTACNNN, NNNATCNNNAGTNNN 

3'-inline-barcode
  - 7nt
  - at the left most of read2 (for Paired-end sequencing)
  - between insert and p7-adapter
  - in read2: CGTGAT{T}
  - in read1: {A}ATCACG
  
3'-adapter
  - TruSeq library structure
  - AGATCGGAAGAGC (for both read1,read2)

===================
Analysis in details

1. remove 3' adapter, discard reads not including adapters 
2. extract UMI+barcode, remove PCR dup
3. wrapping stat, raw, no-adapter, no-UMI, PCR-dup, clean

=====
Tools

1. hiseq.qc.hiseq_lib.HiseqLib(): guess UMI and barcode (top 100000 reads)
1. hiseq.trim.trim_r1.trimR1(): remove adapters
2. umitools: from Weng lab, https://github.com/weng-lab/umitools 
3. python scripts: wrap stat

========
Examples

1. hiseq trim -1 data/raw_data/*1.fq.gz -o data/cut_ad -m 44 --threads 8 --rm-untrim 

2. umitools reformat_sra_fastq -5 NNNCGANNNTACNNN,NNNATCNNNAGTNNN -3 NNNNAATCACG,NNNNACGATGT -i $fq -o $fout -d $fdup 2>$log
"""


import sys
import os
import re
import argparse
import pathlib
import shutil
import tempfile
import hiseq # for report file
from hiseq.qc.hiseq_lib import HiseqLibSmRNA
from hiseq.trim.trim_r1 import TrimR1
from hiseq.utils.utils import log, Config, update_obj, run_shell_cmd, read_hiseq
from hiseq.utils.file import check_path, file_exists, \
    file_abspath, fx_name, symlink_file, check_fx_args
from multiprocessing import Pool


class TrimSmRNA(object):
    """
    TrimSmRNA for multiple fastq files
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
    
    def init_args(self):
        args_init = {
            'fq1': None, # required
            'fq2': None, # optional, for checking 3' UMI/barcode only
            'outdir': None,
            'top_n': 100000,
            'cutoff': 0.8,
            'verbose': False,
            'threads': 4,
            'parallel_jobs': 1,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_rn'
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.args_valid = self.check_fq() # valid
        self.init_files()
        self.rep_list = [os.path.join(self.outdir, i) for i in fx_name(self.fq1, fix_pe=True)]
        Config().dump(self.__dict__, self.config_yaml)
        

    def check_fq(self):
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq1 = [self.fq2]
        return check_fx_args(self.fq1, self.fq2)
    
        
    def init_files(self):
        self.project_dir = self.outdir
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.config_yaml = os.path.join(self.config_dir, 'config.yaml')
        self.report_dir = os.path.join(self.project_dir, 'report')
        check_path([self.config_dir, self.report_dir])
        
        
    def run_single(self, fq1):
        if isinstance(self.fq2, list):
            i = self.fq1.index(fq1)
            fq2 = self.fq2[i]
        else:
            fq2 = None
        args_local = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.outdir,
            'top_n': self.top_n,
            'verbose': self.verbose,
            'threads': self.threads,
            'overwrite': self.overwrite,
        }
        TrimSmRNAR1(**args_local).run()
        
    
    def run_multi(self):
        if self.parallel_jobs > 1 and len(self.fq1) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, self.fq1)
        else:
            for fq1 in self.fq1:
                self.run_single(fq1)
    

    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        hiseq_report_R = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        report_html = os.path.join(
            self.report_dir,
            'HiSeq_report.html')
        report_stdout = os.path.join(self.report_dir, 'report.stdout')
        report_stderr = os.path.join(self.report_dir, 'report.stderr')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            hiseq_report_R,
            self.project_dir,
            self.report_dir,
            '1>{}'.format(report_stdout),
            '2>{}'.format(report_stderr),
            ])
        # save command
        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # report_html
        if file_exists(report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                report_html))
        else:
            run_shell_cmd(cmd)
    
    def run(self):
        if self.args_valid:
            self.run_multi()
            self.report()
        else:
            log.error('TrimSmRNA() exit, illegal fastq files')


class TrimSmRNAR1(object):
    """
    1. trim adapter
    2. extract UMI
    
    directory:
    - trim_ad/
    - extract_umi/
    - trim.json
    - *.fq.gz
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args_init = {
            'fq1': None, # required
            'fq2': None, # optional, for checking 3' UMI/barcode only
            'outdir': None,
            'top_n': 100000,
            'cutoff': 0.8,
            'verbose': False,
            'threads': 4,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'trim_r1'
        self.prefix = fx_name(self.fq1, fix_pe=True)
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.args_valid = self.check_required()
        self.init_files()
        self.guess_smRNA() # update umi5, umi3, barcode
        Config().dump(self.__dict__, self.config_yaml)
        
        
    def check_required(self):
        tools = {
            'fq1': self.fq1,
            'hiseq': shutil.which('hiseq'),
            'cutadapt': shutil.which('cutadapt'),
            'umitools': shutil.which('umitools'),
        }
        flag = all([isinstance(v, str) for k,v in tools.items()])
        msg_list = ['{:14s}: {}'.format(k, v) for k,v in tools.items()]
        msg = '\n'.join([
            '-'*80,
            'Check tools:',
            '\n'.join(msg_list),
            '-'*80,
        ])
        print(msg)
        return flag and self.init_fq(self.fq1)
    
    
    def init_fq(self, fq):
        flag = False
        if isinstance(fq, str):
            if file_exists(self.fq1):
                if re.search('f(ast)?q(.gz)?$', fq, flags=re.IGNORECASE):
                    flag = True
                else:
                    log.error('not like a fastq file: {}'.format(fq))
            else:
                log.error('file not exists: {}'.format(fq))
        else:
            log.error('fq, expect str, got {}'.format(type(fq).__name__))
        if not flag:
            log.error('illegal fastq file, see message above')
        return flag
    
        
    def init_files(self):
        self.project_dir = os.path.join(self.outdir, self.prefix)
        self.config_dir = os.path.join(self.project_dir, 'config')
        default_files = {
            'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
            'guess_umi_dir': os.path.join(self.project_dir, '0_guess_umi'),
            'trim_dir': os.path.join(self.project_dir, '1_trim_ad'),
            'umi_dir': os.path.join(self.project_dir, '2_extract_umi'),
            'exception_dir': os.path.join(self.project_dir, '3_exceptions'),
            'clean_fq': os.path.join(self.project_dir, self.prefix+'.fq.gz'),
            'trim_json': os.path.join(self.project_dir, self.prefix+'.trim.json'),
            'report_dir': os.path.join(self.project_dir, 'report'),
        }
        self = update_obj(self, default_files, force=True)
        check_path([
            self.config_dir, self.trim_dir, self.umi_dir, self.guess_umi_dir,
            self.report_dir
        ])
        
    
    def guess_smRNA(self):
        args_local = {
            'fq1': self.fq1,
            'fq2': self.fq2,
            'outdir': self.guess_umi_dir,
            'top_n': self.top_n,
            'verbose': self.verbose,
        }
        s = HiseqLibSmRNA(**args_local)
        s.run()
        # update: umi5,umi3,barcode
        self.umi5 = s.umi5
        self.umi3 = s.umi3
        self.barcode = s.barcode
    
    
    # step-1.
    def trim_ad(self, fq):
        args_local = {
            'fq1': fq,
            'fq2': None, # skip fq2
            'outdir': self.trim_dir,
            'smp_name': self.prefix,
            'len_min': 44, # 15+18+11
            'len_max': 66, # IP: 15+30+11, total: 15+40+11
            'qual_min': 35,
            'overlap': 3,
            'rm_untrim': True, # remove no-adapter reads
            'threads': self.threads,
        }
        trim = TrimR1(**args_local)
        trim.run()
        return trim.clean_fq1 # no ad
    
    
    # step-2.
    def extract_umi(self, fq):
        print('!AAAA-1', self.umi5, self.umi3, self.barcode)
        # guess umi5, umi3
        if isinstance(self.umi5, list):
            u5 = ','.join(self.umi5)
        elif isinstance(self.umi5, str):
            u5 = self.umi5
        else:
            u5 = None
            log.warning('UMI5 not found')
        # guess umi3
        if isinstance(self.umi3, list):
            u3 = ','.join(self.umi3)
        elif isinstance(self.umi3, str):
            u3 = self.umi3
        else:
            u3 = None
            log.warning('UMI3 not found')
        # guess barcode
        if isinstance(self.barcode, list):
            bc = ','.join(self.barcode)
        elif isinstance(self.barcode, str):
            bc = self.barcode
        else:
            bc = None
            log.warning('barcode not found')
        # choose process
        if u5:
            u3 = u3 if u3 else bc if bc else 'N' # 3end
        else:
            log.error('UMI5 not detected, skipped')
            return fq
        # the main process        
        out_fq = os.path.join(self.umi_dir, self.prefix+'.fq.gz')
        dup_fq = os.path.join(self.umi_dir, self.prefix+'.dup.fq.gz')
        stdout = os.path.join(self.umi_dir, 'run_umitools.stdout')
        stderr = os.path.join(self.umi_dir, 'run_umitools.stderr')
        cmd_txt = os.path.join(self.umi_dir, 'cmd.sh')
        cmd = ' '.join([
            '{}'.format(shutil.which('umitools')),
            'reformat_sra_fastq',
            '-5 {}'.format(u5),
            '-3 {}'.format(u3),
            '-i {}'.format(fq),
            '-o {}'.format(out_fq),
            '-d {}'.format(dup_fq),
            '1> {}'.format(stdout),
            '2> {}'.format(stderr),
        ])
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # check args
        if file_exists(out_fq) and not self.overwrite:
            log.info('extract_umi() skipped, file exists: {}'.format(out_fq))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('extract_umi() failed')
#         # for clean data
#         if file_exists(out_fq):
#             symlink_file(out_fq, self.clean_fq)
        return out_fq
    
    
    # step-3. remove p5 
    def guess_exceptions(self, fq):
        """
        check any exceptions in fq:
        
        - p5 adapters
        - umi
        - ...
        
        Possible: 5-UMI-adapter not in correct position
        
        sequence: 5'-TTCCCTACACGACGCTCTTCCGATCTNNNCGANNNTACNNN-3'
        correct: 5'-5-UMI-adapter-{insert}-3'
        wrong:   5'-5-UMI-adapter-{insert}-[5-UMI-adapter]-3'
        
        How to?
        Trim sequences at 3' end, at least 15+18 = 33 nt
        """
        # check if barcode/umi is None
        if self.umi3 is None and self.barcode is None:
            # use 'N' instead
            # trim 10nt at 3' end,
            cut_after_trim = '0,-10'
        else:
            cut_after_trim = 0
        args_local = {
            'fq1': fq,
            'fq2': None, # skip fq2
            'outdir': self.exception_dir,
            'smp_name': self.prefix,
            'len_min': 18, # 18-40
            'cut_after_trim': cut_after_trim,
            'qual_min': 35,
            'save_untrim': False,
            'overlap': 2,
            'adapter3': 'TTCCCTACACGACGCTCTTCCGATCT', # ACTCT
            'recursive': True,
            'rm_untrim': False, # remove no-adapter reads
            'threads': self.threads,
        }
        trim = TrimR1(**args_local)
        trim.run()
        # return trim.clean_fq1 # no ad
        # for clean data
        if file_exists(trim.clean_fq1):
            symlink_file(trim.clean_fq1, self.clean_fq)
        return trim.clean_fq1
    
    
    def parse_value(self, x):
        """
        Extract value after ":"
        return, [key, value] pair
        
        Example output of: umitools reformat_sra_fastq
        
        ------------------------------------------------------------------------
        Summary of UMI patterns at 5' and 3' ends:
        ------------------------------------------------------------------------
        5' UMI pattern: NNNCGANNNTACNNN
        3' UMI pattern: NNNNACGATGT
        ------------------------------------------------------------------------
        Number of UMI errors allowed: 0
        Input is gzipped.

        Stats: 
        Total input reads:      1607291
        Reads dropped due to improper UMI:      55843
        Final proper read:      1551448
                Reads that are duplicates:      253870
                Reads that are non-duplicates:  1297578
                
        output:
        parse values after ":"
        output: 
         - str, characters
         - int, only digits
         - list, contains ","
        """
        p1 = re.compile('^[A-Za-z]+$') # characters
        p2 = re.compile('^[\\d+,]+$') # digits, 1,234  2158
        out = None
        if isinstance(x, str):
            if ':' in x:
                l = x.strip().split(':')
                v = l.pop() # last item, value
                k = l.pop() # possible key
                v = v.lstrip() # remove spacess on left
                if p2.match(v):
                    v = v.replace(',', '') # remove possible comma, 
                    v = int(v)
                # check if "," exists, multiple items
                elif ',' in v:
                    v = [i.lstrip().rstrip() for i in v.split(',')]
                else:
                    v = str(v)
                out = [k, v]
        return out
    
    
    def wrap_umitools_log(self, x):
        """
        parse values from log, convert to dict
        to-do: allow pattern recognition
        
        {"Summary of UMI patterns at 5' and 3' ends": '',
         "5' UMI pattern": 'NNNCGANNNTACNNN',
         "3' UMI pattern": 'NNNNACGATGT',
         'Number of UMI errors allowed': 0,
         'Stats': '',
         'Total input reads': 1607291,
         'Reads dropped due to improper UMI': 55843,
         'Final proper read': 1551448,
         'Reads that are duplicates': 253870,
         'Reads that are non-duplicates': 1297578}
        """
        if not file_exists(x):
            log.error('umitools log not exists: {}'.format(x))
            return {} # empty
        try:
            with open(x) as r:
                lines = r.readlines()
            s = [self.parse_value(i) for i in lines]
            d = {i[0]:i[1] for i in s if isinstance(i, list)} # to dict
        except IOError as e:
            log.error('failed reading umi log, {}'.format(umi_log))
            d = {} # empty
        out = {
            'umi5': d.get("5' UMI pattern", ""),
            'umi3': d.get("3' UMI pattern", ""),
            'errors': d.get('Number of UMI errors allowed', 0),
            'input': d.get('Total input reads', 0),
            'no_umi': d.get('Reads dropped due to improper UMI', 0),
            'dup': d.get('Reads that are duplicates', 0),
            'clean': d.get('Reads that are non-duplicates', 0)
        }
        return out
    
    
    def wrap_stat(self):
        # 1. trim_ad
        fq_trim_dir = os.path.join(self.trim_dir, self.prefix)
        a = read_hiseq(fq_trim_dir)
        d1 = Config().load(a.trim_json)
        # 2. umitools log
        umi_log = os.path.join(self.umi_dir, 'run_umitools.stderr')
        d2 = self.wrap_umitools_log(umi_log)
        # 3. trim exceptions
        fq_excepion_dir = os.path.join(self.exception_dir, self.prefix)
        b = read_hiseq(fq_excepion_dir)
        d3 = Config().load(b.trim_json)
        # 4. combine
        # name, total, clean, percentage, too_short, too_short2, no_umi, dup
        out = {
            'name': d1.get('name', ''),
            'total': d1.get('total', 0),
            'clean': d3.get('clean', 0),
            'too_short': d1.get('total', 0) - d1.get('clean', 0),
            'too_short2': d3.get('total', 0) - d3.get('clean', 0), # rm_p5
            'no_umi': d2.get('no_umi', 0),
            'dup': d2.get('dup', 0),
        }
        # for divided by zero
        try:
            out['percentage'] = out.get('clean')/out.get('total')
        except ZeroDivisionError as e:
            log.error(e)
            out['percentage'] = 0
        # save to stat
        Config().dump(out, self.trim_json)
    

    def report(self):
        pkg_dir = os.path.dirname(hiseq.__file__)
        hiseq_report_R = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        report_html = os.path.join(
            self.report_dir,
            'HiSeq_report.html')
        report_stdout = os.path.join(self.report_dir, 'report.stdout')
        report_stderr = os.path.join(self.report_dir, 'report.stderr')
        cmd = ' '.join([
            '{}'.format(shutil.which('Rscript')),
            hiseq_report_R,
            self.project_dir,
            self.report_dir,
            '1>{}'.format(report_stdout),
            '2>{}'.format(report_stderr),
            ])
        # save command
        cmd_txt = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # report_html
        if file_exists(report_html) and not self.overwrite:
            log.info('report() skipped, file exists: {}'.format(
                report_html))
        else:
            run_shell_cmd(cmd)


    def run(self):
        if self.args_valid:
            fq_noad = self.trim_ad(self.fq1)
            fq_clean = self.extract_umi(fq_noad)
            fq_x = self.guess_exceptions(fq_clean) # remove P5 adapters
            self.wrap_stat()
            self.report()
        else:
            log.error('TrimSmRNA() faiiled, args illegal')
        

def get_args():
    parser = argparse.ArgumentParser(
        description='trim_smRNA')
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='read1 of FASTQ files')
    parser.add_argument('-2', '--fq2', nargs='+', required=False,
        default=None, help='read2 of FASTQ files, optional')
    parser.add_argument('-o', '--outdir', default=None, required=False,
        help='The directory to save results.')
    parser.add_argument('-n', '--top-n', dest='top_n', type=int, default=100000,
        help='Guess smRNA structure from top_n reads, default: [100000]')
    parser.add_argument('-c', '--cutoff', type=float, default=0.8,
        help='Cutoff for matching UMI/barcode/p7, default: [0.8]')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads for cutadapt --cores, default: [4]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs to run in parallel, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists files')
    parser.add_argument('--debug', dest='verbose', action='store_true',
        help='Print messages in detail, to help debugging')
    return parser
        

def main():
    args = vars(get_args().parse_args())
    TrimSmRNA(**args).run()

    
if __name__ == '__main__':
    main()
        
# 