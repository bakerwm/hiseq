#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Call peaks using MACS2 or SEACR
"""

import os
import re
import shutil
from hiseq.utils.file import check_path, file_exists, file_abspath, \
    file_prefix, remove_file
from hiseq.utils.utils import log, update_obj, Config, get_date, run_shell_cmd
from hiseq.utils.bam import Bam
from hiseq.utils.genome import Genome


class CallPeak(object):
    """
    Call peaks using MACS2 or SEACR

    MACS2: 
    macs2 callpeak -t ip.bam -c input.bam \
      -g hs -f BAMPE -n macs2_peak_q0.1 \
      --outdir outdir -q 0.1 
      --keep-dup all 
      2>log.txt
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        default_args = {
            'method': 'macs2',
            'ip': None,
            'input': None,
            'genome': None,
            'genome_size': None,
            'genome_size_file': None,
            'prefix': None,
            'outdir': None,
            'overwrite': False,
            'hiseq_type': 'callpeak_r1',
            'keep_tmp': False
            }
        self = update_obj(self, default_args, force=False)
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
            log.warning('set pwd as outdir: {}'.format(self.outdir))
        self.outdir = file_abspath(self.outdir)
        if not file_exists(self.ip):
            raise ValueError('ip, required')
        self.ip = file_abspath(self.ip)
        self.input = file_abspath(self.input)
        # file name
        self.ip_name = file_prefix(self.ip)
        if self.input is None:
            self.input_name = ''
        else:
            self.input_name = file_prefix(self.input) 
        if not isinstance(self.prefix, str):
            self.prefix = self.ip_name
        self.prefix_top001 = self.prefix + '.top0.01'
        # genome and size
        self.init_genome()
        self.init_files()
        check_path(self.outdir, create_dirs=True)


    def init_genome(self):
        """
        update:
        genome_size
        genome_size_file
        """
        g = {
          'dm6': 'dm',
          'dm3': 'dm',
          'mm9': 'mm',
          'mm10': 'mm',
          'hg19': 'hs',
          'hg38': 'hs'
        }
        if isinstance(self.genome, str):
            self.genome_size = g.get(self.genome, 0)
            self.genome_size_file = Genome(genome=self.genome).fasize()
            if self.genome_size is None:
                raise ValueError('unknown genome: {}'.format(self.genome))
        elif isinstance(self.genome_size, int):
            if self.genome_size < 1:
                raise ValueError('genome_size, failed, {}'.format(
                    self.genome_size))
            elif not isinstance(self.genome_size_file, str):
                raise ValueError('genome_size_file, failed, {}'.format(
                    self.genome_size_file))
            else:
                pass
        else:
            raise ValueError('--genome, --genome-size, required')


    def init_files(self):
        default_files = {
            'config_toml': self.outdir + '/config.toml',
            'macs2_peak': self.outdir + '/' + self.prefix + '_peaks.narrowPeak',
            'macs2_stdout': self.outdir + '/' + self.prefix + '.macs2.stdout',
            'macs2_stderr': self.outdir + '/' + self.prefix + '.macs2.stderr',
            'macs2_peak_xls': self.outdir + '/' + self.prefix + '_peaks.xls',
            'ip_bed': self.outdir + '/' + self.ip_name + '.frag.bed',
            'ip_bg': self.outdir + '/' + self.ip_name + '.frag.bg',
            'input_bed': self.outdir + '/' + self.input_name + '.frag.bed',
            'input_bg': self.outdir + '/' + self.input_name + '.frag.bg',
            'seacr_peak': self.outdir + '/' + self.prefix + '.stringent.bed',
            'seacr_peak_top': self.outdir + '/' + self.prefix_top001 + '.stringent.bed',
            'seacr_stdout': self.outdir + '/' + self.prefix + '.seacr.stdout',
            'seacr_stderr': self.outdir + '/' + self.prefix + '.seacr.stderr',
        }
        self = update_obj(self, default_files, force=True) # key
        

    def run_macs2(self):
        input_arg = '' if self.input is None else '-c {}'.format(self.input)
        cmd = ' '.join([
            '{} callpeak'.format(shutil.which('macs2')),
            '-t {} {}'.format(self.ip, input_arg),
            '-g {} -f BAMPE'.format(self.genome_size),
            '-n {}'.format(self.outdir + '/' + self.prefix),
            '-q 0.1 --keep-dup all',
            '--nomodel --shift -100 --extsize 200',
            '1> {}'.format(self.macs2_stdout),
            '2> {}'.format(self.macs2_stderr),
            ])
        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.macs2.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # run
        if file_exists(self.macs2_peak) and not self.overwrite:
            log.info('run_macs2() skipped, file exists:{}'.format(
                self.macs2_peak))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('run_macs2() failed, check {}'.format(
                    self.macs2_log))


    def bampe_to_bg(self, bam, bg):
        """
        Convert bam to bed
        bedpe, sort -n (by name)
        """
        bam_sorted = os.path.splitext(bg)[0] + '.sorted_by_name.bam'
        bed = os.path.splitext(bg)[0] + '.bed'
        cmd = ' '.join([
            'samtools sort -@ 4 -n -o {} {}'.format(bam_sorted, bam),
            '&& bedtools bamtobed -bedpe -i {}'.format(bam_sorted),
            "| awk '$1==$4 && $6-$2 < 1000 {print $0}'",
            '| cut -f 1,2,6',
            '| sort -k1,1 -k2,2n > {}'.format(bed),
            '&& bedtools genomecov -bg -i {}'.format(bed),
            '-g {}'.format(self.genome_size_file),
            '> {}'.format(bg)
        ])
        # save cmd
        cmd_txt = self.outdir + '/' + self.prefix + '.bam2bg.cmd.sh'
        with open(cmd_txt, 'wt') as w:
            w.write(cmd + '\n')
        # run
        if file_exists(bg) and not self.overwrite:
            log.info('bampe_to_bg() skipped, file exists:{}'.format(bg))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('bampe_to_bg() failed, check {}'.format(bg))
        # temp files
        del_list = [bam_sorted, bed]
        if not self.keep_tmp:
            remove_file(del_list, ask=False)


    def run_seacr(self):
        """
        bash $seacr ip.bedgraph input.bedgraph \
        non stringent output
    
        see: https://github.com/FredHutch/SEACR
        1. bam-to-bed
        bedtools bamtobed -bedpe -i in.bam > out.bed
        awk '$1==$4 && $6-$2 < 1000 {print $0}' out.bed > out.clean.bed
        cut -f1,2,6 out.clean.bed | sort -k1,1 -k2,2n -k3,3n > out.fragments.bed
        bedtools genomecov -bg -i out.fragments.bed -g genome > out.fragments.bg
        
        2. call peak
        # ip-only
        bash seacr out.fragments.bg 0.01 non stringent top0.01.peaks
        # ip-vs-input
        bash seacr out.fragments.bg IgG.bg non stringent output
        """
#         # ip-only
#         if file_exists(self.ip):
#             cmd1 = ' '.join([
#                 'bash {}'.format(shutil.which('SEACR_1.3.sh')),
#                 '{} 0.01'.format(self.ip_bg),
#                 'non stringent {}'.format(self.outdir + '/' + self.prefix_top001),
#                 '1> {}'.format(self.seacr_stdout),
#                 '2> {}'.format(self.seacr_stderr),
#             ])
#             # save cmd
#             cmd1_txt = self.outdir + '/' + self.prefix + '.SEACR.top0.01.sh'
#             with open(cmd1_txt, 'wt') as w:
#                 w.write(cmd1 + '\n')
#             if file_exists(self.seacr_peak_top) and not self.overwrite:
#                 log.info('run_seacr() skipped, file exists:{}'.format(
#                     self.seacr_peak_top))
#             else:
#                 try:
#                     self.bampe_to_bg(self.ip, self.ip_bg)
#                     run_shell_cmd(cmd1)
#                     remove_file(self.ip_bg, ask=False)
#                 except:
#                     log.error('run_seacr() failed, check {}'.format(
#                         self.seacr_peak_top))
        # ip-vs-input
        if file_exists(self.input):
            cmd2 = ' '.join([
                'bash {}'.format(shutil.which('SEACR_1.3.sh')),
                '{} {}'.format(self.ip_bg, self.input_bg),
                'non stringent {}'.format(self.outdir + '/' + self.prefix),
                '1> {}'.format(self.seacr_stdout),
                '2> {}'.format(self.seacr_stderr),
            ])
            # save cmd
            cmd2_txt = self.outdir + '/' + self.prefix + '.SEACR.ip_vs_input.sh'
            with open(cmd2_txt, 'wt') as w:
                w.write(cmd2 + '\n')
            # run
            if file_exists(self.seacr_peak) and not self.overwrite:
                log.info('run_seacr() skipped, file exists:{}'.format(
                    self.seacr_peak))
            else:
                try:
                    self.bampe_to_bg(self.input, self.input_bg)
                    run_shell_cmd(cmd2)
                    remove_file(self.input_bg, ask=False)
                except:
                    log.error('run_seacr() failed, check {}'.format(
                        self.seacr_peak))


    def run(self):
        Config().dump(self.__dict__, self.config_toml)
        msg = '\n'.join([
            '-'*80,
            '{:>14s} : {}'.format('Date', get_date()),
            '{:>14s} : {}'.format('program', 'hiseq.atac.CallPeak'),
            '{:>14s} : {}'.format('config', self.config_toml),
            '{:>14s} : {}'.format('peakcaller', self.method),
            '{:>14s} : {}'.format('genome', self.genome),
            '{:>14s} : {}'.format('outdir', self.outdir),
            '{:>14s} : {}'.format('ip', self.ip),
            '{:>14s} : {}'.format('input', self.input),
            '{:>14s} : {}'.format('genome_size', self.genome_size),
            '{:>14s} : {}'.format('genome_size_file', self.genome_size_file),
            '{:>14s} : {}'.format('keep_tmp', self.keep_tmp),
            '{:>14s} : {}'.format('overwrite', self.overwrite),
            '{:>14s} : {}'.format('hiseq_type', self.hiseq_type),
            '-'*80,
        ])
        print(msg)
        # generate cmd
        if self.method.lower() == 'macs2':
            self.run_macs2()
        elif self.method.lower() == 'seacr':
            self.run_seacr()
        else:
            raise ValueError('unknown method: {macs2|seacr}, got {}'.format(
                self.method))
        # remove files
        del_list = []
