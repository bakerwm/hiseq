#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
CnR-seq pipeline: level-2 (build design.yaml)
create design.yaml

ip:
ip_fq2:
input:
input_fq2:

search files by key words, from command line
"""

import os
import argparse
from hiseq.utils.file import file_abspath, file_exists, fx_name, list_fx, \
    check_fx_paired
from hiseq.utils.utils import log, update_obj, Config, get_date


class CnrRd(object):
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
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_rd'
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
        ip_fq1 : file or keyword (+)
        ip_fq2 : file or keyword (+)
        input_fq1 : file or keyword (+, optional)
        input_fq2 : file or keyword (+, optional)
        """
        # file exists or keywords
        ## for ip files
        if isinstance(self.ip, list):
            c_ip = all([isinstance(i, str) for i in self.ip])
        elif isinstance(self.ip_fq1, list):
            c2 = isinstance(self.ip_fq2, list)
            c1p = all(check_fx_paired(self.ip_fq1, self.ip_fq2))
            c1s = all([isinstance(i, str) for i in self.ip_fq1 + self.ip_fq2])
            c_ip = all([c2, c1p, c1s])
        else:
            c_ip = False
        ## for input files
        if isinstance(self.input, list):
            c_input = all([isinstance(i, str) for i in self.input])
        elif self.input_fq1 is None:
            c_input = True
        elif isinstance(self.input_fq1, list):
            c4 = isinstance(self.input_fq2, list)
            c3p = all(check_fx_paired(self.input_fq1, self.input_fq2))
            c3s = all([isinstance(i, str) for \
                       i in self.input_fq1 + self.input_fq2])
            c_input = all([c4, c3p, c3s])
        else:
            c_input = False
        # message
        chk = all([c_ip, c_input])
        msg = '\n'.join([
            '='*80,
            'Arugments :',
            '{:>14} : {}'.format('ip', self.ip),
            '{:>14} : {}'.format('input', self.input),
            '{:>14} : {}'.format('ip_fq1', self.ip_fq1),
            '{:>14} : {}'.format('ip_fq2', self.ip_fq2),
            '{:>14} : {}'.format('input_fq1', self.input_fq1),
            '{:>14} : {}'.format('input_fq2', self.input_fq2),
            '-'*40,
            'Status :',
            '{:>20} : {}'.format('ip files or str', c_ip),
            '{:>20} : {}'.format('input files or str', c_input),
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
        if self.input is None:
            self.input = [None] * len(self.ip)
        elif isinstance(self.input, str):
            self.input = [self.input] * len(self.ip)
        elif isinstance(self.input, list):
            if len(self.input) == 1:
                self.input = [self.input[0]] * len(self.ip)
        else:
            self.input = [None] * len(self.ip)
        # split files input groups
        out = {}
        for ka, kb in zip(self.ip, self.input):
            # ip files
            ka_fq = [i for i in f_list if ka in fx_name(i)]
            ka_fq1 = [i for i in ka_fq if fx_name(i).endswith('_1')]
            ka_fq2 = [i for i in ka_fq if fx_name(i).endswith('_2')]
            ka_paired = all(check_fx_paired(ka_fq1, ka_fq2)) 
            # input files
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
                        'ip_fq1': ka_fq1,
                        'ip_fq2': ka_fq2,
                        'input_fq1': kb_fq1,
                        'input_fq2': kb_fq2
                    }
                })
            else:
                log.warning('not paired: {}'.format(k_name))
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
        if isinstance(self.ip_fq1, list):
            ip_paired = all(check_fx_paired(self.ip_fq1, self.ip_fq2))
        else:
            ip_paired = False
            log.warning('--ip-fq1 expect list, got {}'.format(
                type(self.ip_fq1).__name__))
        # input (optional)
        if isinstance(self.input_fq1, list):
            input_paired = all(check_fx_paired(self.input_fq1, self.input_fq2))
        else:
            input_paired = True
            self.input_fq1 = self.input_fq2 = None # force
        # group
        if all([ip_paired, input_paired]):
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
            out = None
        return out

                
    def parse_fq(self):
        """
        locate fastq files, into groups: self.fq_groups

        required:
        1. fq_dir
        2. fq1, fq2
        """
        if isinstance(self.fq_dir, str) and isinstance(self.ip, list):
            fq = self.parse_fq_v1()
        elif isinstance(self.ip_fq1, list) and isinstance(self.ip_fq2, list):
            fq = self.parse_fq_v2()
        else:
            fq = {} # empty
        if len(fq) == 0: # empty
            raise ValueError('no fastq files found, check: fq_dir, ip, ip_fq1, ')
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
            ip_fq1 = v.get('ip_fq1', None)
            ip_fq2 = v.get('ip_fq2', None)
            input_fq1 = v.get('input_fq1', None)
            input_fq2 = v.get('input_fq2', None)
            ip_name = fx_name(ip_fq1[0], fix_pe=True, fix_rep=True)
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
        '1. Generate design.yaml, with --fq-dir',
        '$ python cnr_rd.py -d design.yaml --fq-dir data/raw_data --ip K9 --input IgG',
        '2. Generate design.yaml, with -1 and -2',
        '$ python cnr_rd.py -d design.yaml -1 *1.fq.gz -2 *2.fq.gz',
    ])
    parser = argparse.ArgumentParser(
        prog='cnr_rd',
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
    CnrRd(**args).run()


if __name__ == '__main__':
    main()

#