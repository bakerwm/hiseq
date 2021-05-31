#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
RNA-seq pipeline: level-2 (build design.toml)
create design.toml

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
    Generate design.toml
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
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rd'
        self.build_design = False # force
        self.fq_dir = file_abspath(self.fq_dir)


    def update_design(self):
        # load design/fq_groups
        if file_exists(self.design) and self.append:
            design_dict = Config().load(self.design)
        else:
            design_dict = {}
        # new fq_groups
        fq_flag = []
        msg = ['='*80, 'Check fastq files']
        for g, d in self.fq_groups.items():
            d['fq1'] = file_abspath(d.get('fq1', []))
            d['fq2'] = file_abspath(d.get('fq2', []))
            # status
            fg_flag = []
            msg.append('-'*80)
            msg.append('>{:<10}'.format(g))
            for f1, f2 in zip(d['fq1'], d['fq2']):
                fs = file_exists([f1, f2])
                fg_flag.append(all(fs))
                fs = ['ok' if i else 'failed' for i in fs]
                msg.append('{:>10} : {}'.format('fq_dir', os.path.dirname(f1)))
                msg.append('{:>10} : {}'.format(fs[0], os.path.basename(f1)))
                msg.append('{:>10} : {}'.format(fs[1], os.path.basename(f2)))
            # status
            msg.append('{:<6} : {}'.format('status', all(fg_flag)))
            if d in design_dict.values():
                msg.append('{:<6} : {}'.format('action', 'skipped, file exists'))
                continue
            if all(fg_flag):
                design_dict.update({g:d})
            fq_flag.append(all(fg_flag))
        msg.append('='*80)
        print('\n'.join(msg))
        if not all(fq_flag):
            raise ValueError('fastq files failed, check above log')
        self.fq_groups = design_dict


    def check_fx_args(self):
        """
        Check the args for fastq
        mut_fq1 : file or keyword (+)
        mut_fq2 : file or keyword (+)
        wt_fq1 : file or keyword (+, optional)
        wt_fq2 : file or keyword (+, optional)
        """
        # file exists or keywords
        # for mut files
        if isinstance(self.mut, list):
            c_mut = all([isinstance(i, str) for i in self.mut])
        elif isinstance(self.mut_fq1, list):
            c2 = isinstance(self.mut_fq2, list)
            c1p = all(check_fx_paired(self.mut_fq1, self.mut_fq2))
            c1s = all([isinstance(i, str) for i in self.mut_fq1 + self.mut_fq2])
            c_mut = all([c2, c1p, c1s])
        else:
            c_mut = False
        ## for wt files
        if isinstance(self.wt, list):
            c_wt = all([isinstance(i, str) for i in self.wt])
        elif self.wt_fq1 is None:
            c_wt = True
        elif isinstance(self.wt_fq1, list):
            c4 = isinstance(self.wt_fq2, list)
            c3p = all(check_fx_paired(self.wt_fq1, self.wt_fq2))
            c3s = all([isinstance(i, str) for \
                       i in self.wt_fq1 + self.wt_fq2])
            c_wt = all([c4, c3p, c3s])
        else:
            c_wt = False
        # message
        chk = all([c_mut, c_wt])
        msg = '\n'.join([
            '='*80,
            'Arugments :',
            '{:>14} : {}'.format('mut', self.mut),
            '{:>14} : {}'.format('wt', self.wt),
            '{:>14} : {}'.format('mut_fq1', self.mut_fq1),
            '{:>14} : {}'.format('mut_fq2', self.mut_fq2),
            '{:>14} : {}'.format('wt_fq1', self.wt_fq1),
            '{:>14} : {}'.format('wt_fq2', self.wt_fq2),
            '-'*40,
            'Status :',
            '{:>20} : {}'.format('mut files or str', c_mut),
            '{:>20} : {}'.format('wt files or str', c_wt),
            '-'*40,
            'Results : {}'.format(chk),
        ])
        print(msg)
        # output
        return chk


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
            ka_paired = all(check_fx_paired(ka_fq1, ka_fq2)) 
            # wt files
            if isinstance(kb, str):
                kb_fq = [i for i in f_list if kb in fx_name(i)]
                kb_fq1 = [i for i in kb_fq if fx_name(i).endswith('_1')]
                kb_fq2 = [i for i in kb_fq if fx_name(i).endswith('_2')]
                kb_paired = all(check_fx_paired(kb_fq1, kb_fq2))
            else:
                kb_fq1 = kb_fq2 = None
                kb_paired = True #
            k_name = fx_name(ka_fq1[0], fix_pe=True, fix_rep=True)
            # update or not
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
            mut_paired = all(check_fx_paired(self.mut_fq1, self.mut_fq2))
        else:
            mut_paired = False
            log.warning('--mut-fq1 expect list, got {}'.format(
                type(self.mut_fq1).__name__))
        # wt (optional)
        if isinstance(self.wt_fq1, list):
            wt_paired = all(check_fx_paired(self.wt_fq1, self.wt_fq2))
        else:
            wt_paired = True
            self.wt_fq1 = self.wt_fq2 = None # force
        # group
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
            for a, b in zip(mut_fq1, mut_fq2):
                msg.append('{:>60}'.format('-'*56))
                msg.extend([['fq1', a], ['fq2', b]])
            # for wt
            if wt_fq1 is not None:
                wt_name = fx_name(wt_fq1[0], fix_pe=True, fix_rep=True)
                msg.append('{:>80}'.format('-'*76))
                msg.append(['[wt]', wt_name])
                for a, b in zip(wt_fq1, wt_fq2):
                    msg.append('{:>60}'.format('-'*56))
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
        '1. Generate design.toml, with --fq-dir',
        '$ python rnaseq_rd.py -d design.toml --fq-dir data/raw_data --mut K9 --wt IgG',
        '2. Generate design.toml, with -1 and -2',
        '$ python rnaseq_rd.py -d design.toml -1 *1.fq.gz -2 *2.fq.gz',
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