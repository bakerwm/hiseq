

# -*- coding: utf-8 -*-

"""Quality control for fastq files
fastqc

"""

import os
import sys
import re
import shutil
import logging
import hiseq
from multiprocessing import Pool
from hiseq.utils.seq import Fastx
from hiseq.utils.helper import * # all help functions
import hiseq


class Fastqc(object):
    """"
    Run quality control for fastq files
    run fastqc
    process data using fastqcr/hiseqr ...
    """
    def __init__(self, **kwargs):
        self.update(kwargs) # fresh new
        self.init_args() # update


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
        """
        default values:
        output
        """
        args_default = {
            'fq': None,
            'outdir': str(pathlib.Path.cwd()),
            'parallel_jobs': 1,
            'threads': 1,
            'overwrite': False,
            'fastqc_exe': shutil.which('fastqc'),
            }
        self.update(args_default, force=False) # update missing attrs
        self.outdir = file_abspath(self.outdir)
        check_path(self.outdir)
        # fastq files
        self.fq = self.init_fq(self.fq)
        self.fq = file_abspath(self.fq) # absolute path


    def init_fq(self, fq):
        """
        input is fastq file,
        input is dir of fsatq file
        """
        fq_files = []
        if isinstance(fq, str):
            if os.path.exists(fq):
                fq = file_abspath(fq)
                if os.path.isdir(fq):
                    fx = listfile(fq, "*.fastq")
                    fx += listfile(fq, "*.fastq.gz")
                    fx += listfile(fq, "*.fq")
                    fx += listfile(fq, "*.fq.gz")
                    fq_files += fx
                elif os.path.isfile(fq):
                    fq_files.append(fq)
                else:
                    log.warning('file or path expected, {} found'.format(fq))
            else:
                log.warning('path not exists. {}'.format(fq))
                pass
        elif isinstance(fq, list):
            for i in fq:
                fx = self.init_fq(i)
                fq_files += fx
        else:
            log.warning('str, list expected, {} found'.format(type(fq).__name__))

        return fq_files


    def mission(self):
        """
        Determine the analysis mission: single/multiple
        """
        tag = None
        if isinstance(self.fq, str):
            # file exists
            if check_file(self.fq, emptycheck=True):
                if Fastx(self.fq).is_fastq():
                    tag = 'single'
                    self.fastqc_single(self.fq)
                else:
                    log.warning('not fastq format, {}'.format(x))
            else:
                log.warning('file not exists, or is empty, {}'.format(x))
        elif isinstance(self.fq, list):
            f_out = []
            for fq in self.fq:
                if check_file(fq, emptycheck=True):
                    if Fastx(fq).is_fastq():
                        f_out.append(fq)
                    else:
                        log.warning('not fastq format, {}'.format(fq))
                else:
                    log.warning('file not exists, or is empty, {}'.format(fq))
            # check
            if len(f_out) >= 1:
                tag = 'multiple'
                self.fastqc_multiple(f_out)
            else:
                log.warning('no fq files found')
        else:
            log.warning('str, list expected, {} found'.format(type(self.fq).__name__))

        return tag


    def fastqc_output(self, fq):
        """
        the output of fastq files: html, zip, log
        """
        f_html = os.path.join(self.outdir, fq_name(fq) + '_fastqc.html')
        f_zip = os.path.join(self.outdir, fq_name(fq) + '_fastqc.zip')
        f_log = os.path.join(self.outdir, fq_name(fq) + '_fastqc.log')
        return (f_html, f_zip, f_log)


    def fastqc_single(self, fq):
        """
        Run fastqc for fastq file, single mode
        with parameters
        """
        assert isinstance(fq, str)
        f_html, f_zip, f_log = self.fastqc_output(fq)

        cmd = ' '.join([
            self.fastqc_exe,
            '-t {} -o {}'.format(self.threads, self.outdir),
            '{} 2>{}'.format(fq, f_log)])

        # check
        if file_exists(f_zip) and self.overwrite is False:
            log.info('file exists, fastqc skipped: {}'.format(fq_name(fq)))
        else:
            run_shell_cmd(cmd)


    def fastqc_multiple(self, fq_list):
        """"
        Run fastqc for multiple fastq files, parallel-jobs
        """
        ## Pool() run in parallel
        with Pool(processes=self.parallel_jobs) as pool:
            pool.map(self.fastqc_single, fq_list)


    def report(self):
        """
        Create report for fastq files
        including *zip files
        """
        log.info('Organize all fastqc files')
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'fastqc_report.R')
        report_html = os.path.join(self.outdir, 'fastqc_report.html')
        cmd_file = os.path.join(self.outdir, 'cmd.sh')

        cmd = ' '.join([
            shutil.which('Rscript'),
            qc_reportR,
            self.outdir,
            self.outdir])

        with open(cmd_file, 'wt') as w:
            w.write(cmd + '\n')

        if os.path.exists(report_html) and self.overwrite is False:
            log.info('file exists, skip generating html.')
        else:
            run_shell_cmd(cmd)

        if not os.path.exists(report_html):
            log.error('failed, generating html file')

        # finish
        log.info('output - {}'.format( report_html))


    def run(self):
        """
        Run all
        """
        tag = self.mission() # mode
        if tag:
            self.report()

