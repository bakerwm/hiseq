#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function, calculate the IDR, between peak files

Using the command: idr
"""

import os


class PeakIDR(object):
    """Calculate the IDR for peak files
    Irreproducibility Discovery Rate
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'idr_cmd': shutil.which('idr'),
            'peak': None,
            'outdir': None,
            'prefix': None,
            'input_type': 'narrowPeak', # broadPeak, bed, gff
            'cor_method': 'pearson', # spearman
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.check_flag = False
        if file_exists(self.idr_cmd) or is_cmd(self.idr_cmd):
            pass
        else:
            log.error('idr command not exists: {}'.format(self.idr_cmd))
            self.check_flag = False
        if isinstance(self.peak, list):
            if file_exists(self.peak):
                if len(self.peak) > 1:
                    self.check_flag = True
                else:
                    log.info('IDR skipped, require at least 2 files')
            else:
                log.warning('peak not exists, {}'.format(self.peak))
        else:
            log.warning('peak not valid, expect list, got {}'.format(
                type(self.peak).__name__))
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        if not isinstance(self.prefix, str):
            self.prefix = 'idr'
        input_fmt = ['narrowPeak', 'broadPeak', 'bed', 'gff']
        if self.input_type not in input_fmt:
            log.warning('input_type not valid, expect {}, got {}'.format(
                input_fmt, self.input_type))
            self.check_flag = False
        if self.cor_method not in ['pearson', 'spearman']:
            log.warning('cor_method not valid, expect {}, got {}'.format(
                ['pearson', 'spearman'], self.cor_method))
        self.config_toml = os.path.join(self.outdir, 'config.toml')
        

    def run_idr(self, peakA=None, peakB=None):
        """Prepare the command line for two peaks
        example:
        idr --samples peakA peakB --input-file-type narrowPeak --plot --out-file idr.txt --log-output-file log.stdout
        """
        nameA, nameB = file_prefix([peakA, peakB])
        prefix = '{}.{}.vs.{}.idr'.format(self.prefix, nameA, nameB)
        # files
        cmd_shell = os.path.join(self.outdir, prefix + '.cmd.sh')
        idr_txt = os.path.join(self.outdir, prefix + '.txt')
        idr_png = idr_txt + '.png' # os.path.join(self.outdir, prefix + '.png')
        log_stdout = os.path.join(self.outdir, prefix + '.log.stdout')
        lot_stderr = os.path.join(self.outdir, prefix + '.log.stderr')
        cmd = ' '.join([
            '{}'.self.idr_cmd,
            '--samples {} {}'.format(peakA, peakB),
            '--input-file-type {}'.format(self.input_type),
            '--plot', #.format(idr_png),
            '--output-file {}'.format(idr_txt),
            '--log-output-file {}'.format(log_stdout),
            '2> {}'.format(log_stderr)
            ])
        with open(cmd_shell, 'wt') as w:
            w.write(cmd + '\n')
        if os.path.exists(idr_png) and not self.overwrite:
            log.warning('idr skipped, file exists: {}'.format(idr_png))
        else:
            # require, peaks > 100
            rowA = file_nrows(peakA)
            rowB = file_nrows(peakB)
            if rowA >= 100 and rowB >= 100:
                run_shell_cmd(cmd)
            else:
                log.warning('at least 100 peaks required, got {}, {}'.format(
                    rowA, rowB))
        # check
        if not file_exists(idr_png):
            log.error('output not exists: {}'.format(idr_png))


    def run(self):
        """
        Run A, B
        """
        if self.check_flag:
            for peakA, peakB in combinations(self.peak, 2):
                self.run_idr(peakA, peakB)
        else:
            log.error('args not valid')
                
                