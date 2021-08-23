#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Chipseq-seq pipeline: level-2 (build design.yaml)
create design.yaml

ip:
ip_fq2:
input:
input_fq2:

search files by key words, from command line
"""

import os
import argparse
from hiseq.utils.file import check_fx_args, check_fx_paired, file_abspath, \
    file_exists, fx_name, list_fx
from hiseq.utils.utils import log, update_obj, Config, get_date


class ChipseqRd(object):
    """
    Generate design.yaml
    format:
    ip:
    ip_fq2:
    input:
    input_fq2:
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'design': None,
            'fq_dir': None,
            'ip': None, # str
            'input': None, # str
            'ip_fq1': None,
            'ip_fq2': None,
            'input_fq1': None,
            'input_fq2': None,
            'append': True,
            'as_se': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'chipseq_rd'
        self.build_design = False # force
        self.fq_dir = file_abspath(self.fq_dir)


    def parse_fq_v1(self):
        """
        Version-1:
        parse fastq files by keyword
        fq-dir
        ip (required)
        input (optional)
        
        return: dict (groups)
        """
        # check fq_dir
        if isinstance(self.fq_dir, str):
            f_list = list_fx(self.fq_dir)
        else:
            log.warning('--fq-dir required')
            f_list = []
        # check files
        if len(f_list) < 2:
            log.error('not enough fq files')
            return None
        # check ip
        if not isinstance(self.ip, list):
            log.error('unknown ip, expect list, got {}'.format(
                type(self.ip).__name__))
            return None
        # check input
        if isinstance(self.input, str):
            self.input = [self.input] * len(self.ip)
        elif isinstance(self.input, list):
            if len(self.input) == 1:
                self.input = [self.input] * len(self.ip)
        else:
            log.error('--input expect str, list, got {}'.format(
                type(self.input).__name__))
            return None
        # check wt mut equal
        if not len(self.ip) == len(self.input):
            log.error('ip, input not equal in length')
            return None
        out = {}
        for ka, kb in zip(self.ip, self.input):
            # ip files
            ka_fq = [i for i in f_list if ka in fx_name(i)]
            ka_fq = [i for i in ka_fq if file_exists(i)]
            ka_fq1 = [i for i in ka_fq if fx_name(i).endswith('_1')]
            ka_fq2 = [i for i in ka_fq if fx_name(i).endswith('_2')]
            ka_paired = check_fx_paired(ka_fq1, ka_fq2)
            k_name = fx_name(ka_fq[0], fix_pe=True, fix_rep=True)
            # input files
            kb_fq = [i for i in f_list if kb in fx_name(i)]
            kb_fq = [i for i in kb_fq if file_exists(i)]
            kb_fq1 = [i for i in kb_fq if fx_name(i).endswith('_1')]
            kb_fq2 = [i for i in kb_fq if fx_name(i).endswith('_2')]
            kb_paired = check_fx_paired(kb_fq1, kb_fq2)
            # as SE
            if self.as_se or not all([ka_paired, kb_paired]):
                if len(ka_fq1) == 0:
                    ka_fq1 = ka_fq
                if len(kb_fq1) == 0:
                    kb_fq1 = kb_fq
                ka_fq2 = kb_fq2 = None
            out.update({
                k_name: {
                    'mut_fq1': ka_fq1,
                    'mut_fq2': ka_fq2,
                    'wt_fq1': kb_fq1,
                    'wt_fq2': kb_fq2
                }
            })
        # output
        return out
            
    
    def parse_fq_v2(self):
        """
        Version-2:
        parse fastq files from --ip-fq1, ...
        fq-dir
        ip (required)
        input (optional)
        
        return: dict (groups)
        """
        # ip
        ka = check_fx_args(fq1=self.ip_fq1, fq2=self.ip_fq2)
        # input
        kb = check_fx_args(fq1=self.input_fq1, fq2=self.input_fq2)
        if all([ka, kb]):
            k_name = fx_name(self.ip_fq1[0], fix_pe=True, fix_rep=True)
            out = {
                k_name: {
                'ip_fq1': self.ip_fq1,
                'ip_fq2': self.ip_fq2,
                'input_fq1': self.input_fq1,
                'input_fq2': self.input_fq2
                }
            }
        else:
            log.error('--ip-fq1, --input-fq1 files not valid')
            out = None
        return None

                
    def parse_fq(self):
        """
        locate fastq files, into groups: self.fq_groups
        required:
        1. fq_dir
        2. fq1, fq2
        """
        if isinstance(self.fq_dir, str) and isinstance(self.ip, list) \
            and isinstance(self.input, list):
            fq = self.parse_fq_v1()
        elif isinstance(self.ip_fq1, list) and isinstance(self.ip_fq2, list):
            fq = self.parse_fq_v2()
        else:
            fq = {} # empty
        if len(fq) == 0: # empty
            # raise ValueError('no fastq files found, check: fq_dir, ip, ip_fq1, ')
            log.error('no fastq files found, check: fq_dir, ip, ip_fq1, ')
            return {} # empty
        # loading config
        if file_exists(self.design):
            out = Config().load(self.design)
        else:
            out = {}
        if not isinstance(fq, dict):
            return None
        # update fq
        msg = ['='*80, 'build design']
        for k, v in fq.items():
            if k in out:
                status = 'skipped, file exists'
            else:
                out.update({k:v}) # update
                status = 'add'
            # for log
            ip_fq1 = v.get('ip_fq1', None)
            ip_fq2 = v.get('ip_fq2', None)
            input_fq1 = v.get('input_fq1', None)
            input_fq2 = v.get('input_fq2', None)
            ip_name = fx_name(ip_fq1[0], fix_pe=True, fix_rep=True)
            # fix for SE
            if ip_fq2 is None:
                ip_fq2 = [None]*len(ip_fq1)
            if input_fq2 is None:
                input_fq2 = [None]*len(input_fq1)
            # for ip
            msg.append('-'*80)
            msg.append(['[ip]', ip_name])
            for a, b in zip(ip_fq1, ip_fq2):
                msg.append('{:>60}'.format('-'*56))
                msg.extend([['fq1', a], ['fq2', b]])
            # for input
            if input_fq1 is not None:
                input_name = fx_name(input_fq1[0], fix_pe=True, fix_rep=True)
                msg.append('{:>80}'.format('-'*76))
                msg.append(['[input]', input_name])
                for a, b in zip(input_fq1, input_fq2):
                    msg.append('{:>40}'.format('-'*36))
                    msg.extend([['fq1', a], ['fq2', b]])
            # status            
            msg.append('{:>20}'.format('-'*16))
            msg.append(['[Status]', status])
        # last
        msg.append('='*80)
        # format msg
        msg_fmt = []
        for i in msg:
            if isinstance(i, list):
                msg_fmt.append('{:>14} : {}'.format(i[0], i[1]))
            elif isinstance(i, str):
                msg_fmt.append(i)
            else:
                continue
        # out
        print('\n'.join(msg_fmt))
        return out 


    def run(self):
        self.fq_groups = self.parse_fq()
        Config().dump(self.fq_groups, self.design)


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. Generate design.yaml, with --fq-dir',
        '$ python chipseq_rd.py -d design.yaml --fq-dir data/raw_data --ip K9 --input IgG',
        '2. Generate design.yaml, with -1 and -2',
        '$ python chipseq_rd.py -d design.yaml -1 *1.fq.gz -2 *2.fq.gz',
    ])
    parser = argparse.ArgumentParser(
        prog='chipseq_rd',
        description='Generate design',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', default=None, required=True,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('--ip', nargs='+', dest='ip', default=None,
        help='keyword of IP fastq file, auto-find read1/2')
    parser.add_argument('--input', nargs='+', dest='input', default=None,
        help='keyword of Input fastq file, auto-find read1/2')
    parser.add_argument('--se', dest='as_se', action='store_true',
        help='choose only fastq1 of PE reads')
        
    parser.add_argument('--ip-fq1', nargs='+', dest='ip_fq1', default=None,
        help='filepath or keyword of IP fastq file, read1 of PE')
    parser.add_argument('--ip-fq2', nargs='+', dest='ip_fq2', default=None,
        help='filepath or keyword of IP fastq file, read2 of PE')
    parser.add_argument('--input-fq1', nargs='+', dest='input_fq1', 
        default=None,
        help='filepath or keyword of Input fastq file, read1 of PE')
    parser.add_argument('--input-fq2', nargs='+', dest='input_fq2',
        default=None,
        help='filepath or keyword of Input fastq file, read2 of PE')
    return parser


def main():
    args = vars(get_args().parse_args())
    ChipseqRd(**args).run()


if __name__ == '__main__':
    main()

#