#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-2 (build design.toml)
create design.toml
"""

import os
import argparse
from hiseq.utils.file import file_abspath, file_exists, fx_name, list_fx, \
    check_fx_paired
from hiseq.utils.utils import log, update_obj, Config, get_date


class AtacRd(object):
    """
    Generate design.toml
    format:
    input: fq1, fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
          'design': None,
          'fq_dir': None,
          'fq1': None,
          'fq2': None,
          'append': True,
        }
        self = update_obj(self, args_init, force=False)
        self.atacseq_type = 'atacseq_rd'
        self.build_design = False # force


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


    def parse_fq(self):
        """
        locate fastq files, into groups: self.fq_groups

        required:
        1. fq_dir
        2. fq1, fq2
        """
        if isinstance(self.fq_dir, str):
            f_list = list_fx(self.fq_dir)
        elif isinstance(self.fq1, list) and isinstance(self.fq2, list):
            f_list = self.fq1 + self.fq2
        elif isinstance(self.fq1, str) and isinstance(self.fq2, str):
            f_list = [self.fq1, self.fq2]
        else:
            f_list = []
        if len(f_list) < 2: # paired
            raise ValueError('fq_dir, fq1,fq2, failed')
        # recognize samples by name: filename_repx.fq.gz
        out = {}
        f_names = fx_name(f_list, fix_pe=True, fix_rep=True)
        f_names = sorted(list(set(f_names)))
        for f in f_names:
            fq1 = []
            fq2 = []
            for i in f_list:
                iname = fx_name(i, fix_pe=True, fix_rep=True)
                iname2 = fx_name(i, fix_pe=False, fix_rep=False)
                if iname == f:
                    if iname2.endswith('1'):
                        fq1.append(i)
                    elif iname2.endswith('2'):
                        fq2.append(i)
                    else:
                        pass #
            # check fq1, fq2 paired
            if not check_fx_paired(fq1, fq2):
                log.error('fq not paired, {}'.format(f))
                continue
            # save
            out[f] = {
                'fq1': fq1,
                'fq2': fq2,
                'group': f,
            }
        # output
        return out


    def run(self):
        self.fq_groups = self.parse_fq()
        self.update_design()
        Config().dump(self.fq_groups, self.design)


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. Generate design.toml, with --fq-dir',
        '$ python atac_rd.py -d design.toml --fq-dir data/raw_data',
        '2. Generate design.toml, with -1 and -2',
        '$ python atac_rd.py -d design.toml -1 *1.fq.gz -2 *2.fq.gz',
    ])
    parser = argparse.ArgumentParser(
        prog='atac_rd',
        description='Generate design',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', default=None, required=True,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('-1', '--fq1', nargs='+', default=None,
        help='read1 files, (or read1 of PE reads)')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='read2 of PE reads')
    return parser


def main():
    args = vars(get_args().parse_args())
    AtacRd(**args).run()


if __name__ == '__main__':
    main()

#