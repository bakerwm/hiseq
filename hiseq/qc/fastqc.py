#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality control for fastq files

1. fastqc
2. library structure
3. bacteria

"""

import os
import sys
import re
import pathlib
import shutil
import logging
import argparse
import hiseq
from multiprocessing import Pool
from hiseq.utils.seq import Fastx
from hiseq.utils.utils import update_obj, log, run_shell_cmd
from hiseq.utils.file import list_fx, fx_name, file_abspath, \
    check_path, file_exists
# from hiseq.utils.helper import * # all help functions
import hiseq


class Fastqc(object):
    """"
    Run fastqc for fastq files
    process data using fastqcr/hiseqr ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        default values:
        output
        """
        args_init = {
            'fq': None,
            'outdir': str(pathlib.Path.cwd()),
            'parallel_jobs': 1,
            'threads': 1,
            'overwrite': False,
            'fastqc_exe': shutil.which('fastqc')
        }
        self = update_obj(self, args_init, force=False)
        self.outdir = file_abspath(self.outdir)
        check_path(self.outdir)
        # fastq files
        self.fq = self.init_fq(self.fq)
        self.fq = file_abspath(self.fq) # absolute path


    def init_fq(self, fq):
        """
        input is fastq file,
        input is dir of fastq file
        """
        out = []
        if isinstance(fq, str):
            if os.path.exists(fq):
                if os.path.isdir(fq):
                    out = list_fx(fq)
                elif os.path.isfile(fq):
                    out = [fq]
                else:
                    log.warning('file or path expected, got {}'.format(fq))
            else:
                log.warning('path not exists. {}'.format(fq))
        elif isinstance(fq, list):
            for i in fq:
                out += self.init_fq(i)
        else:
            log.warning('str, list expected, got {}'.format(type(fq).__name__))
        return out

    
    def fastqc_output(self, fq):
        """
        the output of fastq files: html, zip, log
        """
        f_html = os.path.join(self.outdir, fx_name(fq) + '_fastqc.html')
        f_zip = os.path.join(self.outdir, fx_name(fq) + '_fastqc.zip')
        f_log = os.path.join(self.outdir, fx_name(fq) + '_fastqc.log')
        return (f_html, f_zip, f_log)


    def run_single_fq(self, fq):
        """
        Run for single fastq file
        """
        # the name of output files
        f_html, f_zip, f_log = self.fastqc_output(fq)
        cmd = ' '.join([
            self.fastqc_exe,
            '-t {} -o {}'.format(self.threads, self.outdir),
            '{} 2>{}'.format(fq, f_log)
        ])
        # check
        if file_exists(f_zip) and self.overwrite is False:
            log.info('fastqc() skipped, file exists: {}'.format(fx_name(fq)))
        else:
            run_shell_cmd(cmd)
    
    
    def run_multiple_fq(self):
        """"
        Run fastqc for multiple fastq files, parallel-jobs
        """
        if self.parallel_jobs > 1 and len(self.fq) > 1:
            ## Pool() run in parallel
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fq, self.fq)
        else:
            [self.run_single_fq(fq) for fq in self.fq]


    def report(self):
        """
        Create report for fastq files
        including *zip files
        """
        # the template: R package
        log.info('Organize all fastqc files')
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'fastqc_report.R')
        report_html = os.path.join(self.outdir, 'fastqc_report.html')
        log_stdout = os.path.join(self.outdir, 'report.stdout')
        log_stderr = os.path.join(self.outdir, 'report.stderr')
        cmd_file = os.path.join(self.outdir, 'report.sh')
        cmd = ' '.join([
            shutil.which('Rscript'),
            qc_reportR,
            self.outdir,
            self.outdir,
            '1> {}'.format(log_stdout),
            '2> {}'.format(log_stderr),
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(report_html) and self.overwrite is False:
            log.info('fastqc_report() skipped, file exists')
        else:
            run_shell_cmd(cmd)
        if not os.path.exists(report_html):
            log.error('fastqc_report() failed, check {}'.format(log_stderr))
        # finish
        log.info('output - {}'.format(report_html))


    def run(self):
        # 1. run fastqc
        self.run_multiple_fq()
        # 2. generate report
        self.report()


def get_args():
    parser = argparse.ArgumentParser(
        description='hiseq qc, fastqc')
    parser.add_argument('-i', '--fq', nargs='+', required=True,
        help='reads in FASTQ files, or directory contains fastq files')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')
    parser.add_argument('--fastqc', default='fastqc',
        help='The path to the fastqc command, default: [fastqc]')
    parser.add_argument('-f', '--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    return parser


def main():
    args = vars(get_args().parse_args())
    Fastqc(**args).run()
    
    
if __name__ == '__main__':
    main()
    
#