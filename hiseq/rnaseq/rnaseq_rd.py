#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
RNA-seq pipeline: level-2 (build design.yaml)
create design.yaml

mut
wt

mut_fq1
mut_fq2
wt_fq1
wt_fq2

search files by key words, from command line
"""

import os
import argparse
from hiseq.utils.file import file_abspath, file_exists, fx_name, list_fx, \
    check_fx_paired
from hiseq.utils.utils import log, update_obj, Config, get_date


class RnaseqRd(object):
    """
    Generate design.yaml
    format:
    mut
    wt
    mut_fq1
    mut_fq2
    wt_fq1
    wt_fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'design': None,
            'fq_dir': None,
            'mut': None, # str
            'wt': None, # str
            'mut_fq1': None,
            'mut_fq2': None,
            'wt_fq1': None,
            'wt_fq2': None,
            'append': True,
            'as_se': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rd'
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
        # check ip
        if not isinstance(self.ip, list):
            return None
        # check input
        if isinstance(self.input, str):
            self.input = [self.input] * len(self.ip)
        elif isinstance(self.input, list):
            if len(self.input) == 1:
                self.input = [self.input[0]] * len(self.ip)
        else:
            log.error('--input expect str, list, got {}'.format(
                type(self.input).__name__))
            return None
        # split files input groups
        out = {}
        for ka, kb in zip(self.ip, self.input):
            # ip files
            ka_fq = [i for i in f_list if ka in fx_name(i)]
            ka_fq = [i for i in ka_fq if file_exists(i)]
            ka_fq1 = [i for i in ka_fq if fx_name(i).endswith('_1')]
            ka_fq2 = [i for i in ka_fq if fx_name(i).endswith('_2')]
            ka_paired = check_fx_paired(ka_fq1, ka_fq2)
            # input files
            kb_fq = [i for i in f_list if kb in fx_name(i)]
            kb_fq = [i for i in kb_fq if file_exists(i)]
            kb_fq1 = [i for i in kb_fq if fx_name(i).endswith('_1')]
            kb_fq2 = [i for i in kb_fq if fx_name(i).endswith('_2')]
            kb_paired = check_fx_paired(kb_fq1, kb_fq2)
            if len(ka_fq) > 0 and len(kb_fq) > 0:
                if all([ka_paired, kb_paired]) and not self.as_se:
                    k_name = fx_name(ka_fq[0], fix_pe=True, fix_rep=True)
                    out.update({
                        k_name: {
                            'mut_fq1': ka_fq1,
                            'mut_fq2': ka_fq2,
                            'wt_fq1': kb_fq1,
                            'wt_fq2': kb_fq2,
                        }
                    })
                else:
                    k_name = fx_name(ka_fq[0], fix_pe=False, fix_rep=True)
                    out.update({
                        k_name: {
                            'mut_fq1': ka_fq1,
                            'mut_fq2': None,
                            'wt_fq1': kb_fq1,
                            'wt_fq2': None,
                        }
                    })
        # output
        return out
                    
        
    def parse_fq_v1(self):
        """
        Version-1:
        parse fastq files by keyword
        fq-dir
        mut (required)
        wt (optional)
        
        return: dict (groups)
        """
        # check fq_dir
        if isinstance(self.fq_dir, str):
            f_list = list_fx(self.fq_dir)
        else:
            log.warning('--fq-dir required')
            f_list = []
        # check mut
        if not isinstance(self.mut, list):
            return None
        # check wt
        if self.wt is None:
            self.wt = [None] * len(self.mut)
        elif isinstance(self.wt, str):
            self.wt = [self.wt] * len(self.mut)
        elif isinstance(self.wt, list):
            if len(self.wt) == 1:
                self.wt = [self.wt[0]] * len(self.mut)
        else:
            self.wt = [None] * len(self.mut)
        # split files wt groups
        out = {}
        for ka, kb in zip(self.mut, self.wt):
            # mut files
            ka_fq = [i for i in f_list if ka in fx_name(i)]
            ka_fq1 = [i for i in ka_fq if fx_name(i).endswith('_1')]
            ka_fq2 = [i for i in ka_fq if fx_name(i).endswith('_2')]
            ka_paired = check_fx_paired(ka_fq1, ka_fq2) 
            # wt files
            if isinstance(kb, str):
                kb_fq = [i for i in f_list if kb in fx_name(i)]
                kb_fq1 = [i for i in kb_fq if fx_name(i).endswith('_1')]
                kb_fq2 = [i for i in kb_fq if fx_name(i).endswith('_2')]
                kb_paired = check_fx_paired(kb_fq1, kb_fq2)
            else:
                kb_fq1 = kb_fq2 = None
                kb_paired = True #
            k_name = fx_name(ka_fq1[0], fix_pe=True, fix_rep=True)
            # update or not
            # force to SE
            if self.as_se:
                ka_fq2 = kb_fq2 = None
            if all([ka_paired, kb_paired]):
                out.update({
                    k_name: {
                        'mut_fq1': ka_fq1,
                        'mut_fq2': ka_fq2,
                        'wt_fq1': kb_fq1,
                        'wt_fq2': kb_fq2
                    }
                })
            else:
                log.warning('not paired: {}'.format(k_name))
        # output
        return out
            
    
    def parse_fq_v2(self):
        """
        Version-2:
        parse fastq files from --mut-fq1, ...
        fq-dir
        mut (required)
        wt (optional)
        
        return: dict (groups)
        """
        # mut
        if isinstance(self.mut_fq1, list):
            mut_paired = check_fx_paired(self.mut_fq1, self.mut_fq2)
        else:
            mut_paired = False
            log.warning('--mut-fq1 expect list, got {}'.format(
                type(self.mut_fq1).__name__))
        # wt (optional)
        if isinstance(self.wt_fq1, list):
            wt_paired = check_fx_paired(self.wt_fq1, self.wt_fq2)
        else:
            wt_paired = True
            self.wt_fq1 = self.wt_fq2 = None # force
        # group
        if self.as_se:
            self.mut_fq2 = self.wt_fq2 = None
        if all([mut_paired, wt_paired]):
            k_name = fx_name(self.mut_fq1[0], fix_pe=True, fix_rep=True)
            out = {
                k_name: {
                'mut_fq1': self.mut_fq1,
                'mut_fq2': self.mut_fq2,
                'wt_fq1': self.wt_fq1,
                'wt_fq2': self.wt_fq2
                }
            }
        else:
            out = None
        return out

                
    def parse_fq(self):
        """
        locate fastq files, into groups: self.fq_groups

        required:
        1. fq_dir
        2. fq1, fq2
        """
        # convert mut/wt to list
        if isinstance(self.mut, str):
            self.mut = [self.mut]
        if isinstance(self.wt, str):
            self.wt = [self.wt]
        if isinstance(self.mut_fq1, str):
            self.mut_fq1 = [self.mut_fq1]
        if isinstance(self.mut_fq2, str):
            self.mut_fq2 = [self.mut_fq2]
        if isinstance(self.wt_fq1, str):
            self.wt_fq1 = [self.wt_fq1]
        if isinstance(self.wt_fq2, str):
            self.wt_fq2 = [self.wt_fq2]
        # choose sub-command
        if isinstance(self.fq_dir, str) and isinstance(self.mut, list):
            fq = self.parse_fq_v1()
        elif isinstance(self.mut_fq1, list) and isinstance(self.mut_fq2, list):
            fq = self.parse_fq_v2()
        else:
            fq = {} # empty
        if len(fq) == 0: # empty
            raise ValueError('no fastq files found, check: fq_dir, mut, mut_fq1, ')
        # loading config
        if file_exists(self.design):
            out = Config().load(self.design)
        else:
            out = {}
        # update fq
        msg = ['='*80, 'build design']
        for k, v in fq.items():
            if k in out:
                status = 'skipped, file exists'
            else:
                out.update({k:v}) # update
                status = 'add'
            # for log
            mut_fq1 = v.get('mut_fq1', None)
            mut_fq2 = v.get('mut_fq2', None)
            wt_fq1 = v.get('wt_fq1', None)
            wt_fq2 = v.get('wt_fq2', None)
            mut_name = fx_name(mut_fq1[0], fix_pe=True, fix_rep=True)
            # for mut
            msg.append('-'*80)
            msg.append(['[mut]', mut_name])
            if isinstance(mut_fq2, list):
                for a, b in zip(mut_fq1, mut_fq2):
                    msg.append('{:>60}'.format('-'*56))
                    msg.extend([['fq1', a], ['fq2', b]])
            else:
                for a in mut_fq1:
                    msg.append('{:>60}'.format('-'*56))
                    msg.extend([['fq1', a], ['fq2', None]])
            # for wt
            if wt_fq1 is not None:
                wt_name = fx_name(wt_fq1[0], fix_pe=True, fix_rep=True)
                msg.append('{:>80}'.format('-'*76))
                msg.append(['[wt]', wt_name])
                if isinstance(wt_fq2, list):
                    for a, b in zip(wt_fq1, wt_fq2):
                        msg.append('{:>60}'.format('-'*56))
                        msg.extend([['fq1', a], ['fq2', b]])
                else:
                    for a in wt_fq1:
                        msg.append('{:>60}'.format('-'*56))
                        msg.extend([['fq1', a], ['fq2', None]])
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
        '$ python rnaseq_rd.py -d design.yaml --fq-dir data/raw_data --mut K9 --wt IgG',
        '2. Generate design.yaml, with -1 and -2',
        '$ python rnaseq_rd.py -d design.yaml -1 *1.fq.gz -2 *2.fq.gz',
    ])
    parser = argparse.ArgumentParser(
        prog='rnaseq_rd',
        description='Generate design',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', default=None, required=True,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('--mut', nargs='+', dest='mut', default=None,
        help='keyword of mut fastq file, auto-find read1/2')
    parser.add_argument('--wt', nargs='+', dest='wt', default=None,
        help='keyword of wt fastq file, auto-find read1/2')
    parser.add_argument('--se', dest='as_se', action='store_true',
        help='choose only fastq1 of PE reads')
    # details
    parser.add_argument('--mut-fq1', nargs='+', dest='mut_fq1', default=None,
        help='filepath or keyword of mut fastq file, read1 of PE')
    parser.add_argument('--mut-fq2', nargs='+', dest='mut_fq2', default=None,
        help='filepath or keyword of mut fastq file, read2 of PE')
    parser.add_argument('--wt-fq1', nargs='+', dest='wt_fq1', 
        default=None,
        help='filepath or keyword of wt fastq file, read1 of PE')
    parser.add_argument('--wt-fq2', nargs='+', dest='wt_fq2',
        default=None,
        help='filepath or keyword of wt fastq file, read2 of PE')
    return parser


def main():
    args = vars(get_args().parse_args())
    RnaseqRd(**args).run()


if __name__ == '__main__':
    main()

#