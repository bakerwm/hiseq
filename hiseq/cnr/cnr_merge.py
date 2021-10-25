#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Combine/merge multiple cnr report

cnr_r1, cnr_rn, ...
"""

import os
import pathlib
import argparse
from hiseq.cnr.cnr_rp import CnrRp
from hiseq.utils.file import check_path, file_abspath, list_dir
from hiseq.utils.utils import log, update_obj, Config, is_hiseq_dir
from hiseq.cnr.utils import copy_hiseq_qc


class CnrMerge(object):
    def __init__(self, **kwargs):
        c = CnrMergeConfig(**kwargs)
        self = update_obj(self, c.__dict__, force=True)
        
        
    def run(self):
        Config().dump(self.__dict__, self.config_yaml) # save config
        Config().dump({'indir': self.indir}, self.indir_json) # save dir list
        copy_hiseq_qc(self.project_dir) # copy files, config
        CnrRp(self.outdir).run() # generate report


class CnrMergeConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'indir': None, # required
            'outdir': None,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_merge'
        if self.outdir is None:
            self.outdir = os.path.join(str(pathlib.Path.cwd()), 'hiseq_merge')
        self.indir = file_abspath(self.indir)
        self.outdir = file_abspath(self.outdir)
        self.hiseq_list = self.parse_dir(self.indir)
        self.init_files()
    
    
    def parse_dir(self, x, hiseq_type='cnr_'):
        """
        Update indir, for cnr_r1, cnr_rn, cnr_rx, ...
        level-1:
        level-2:        
        indir could be dir or a list of dirs saving in a file
        """
        out = []
        if isinstance(x, str):
            if os.path.isdir(x):
                # level-1, skip config, report, ...
                if is_hiseq_dir(x, hiseq_type):
                    out.append(x)
                # level-2
                d = [i for i in list_dir(x, include_dir=True) if is_hiseq_dir(i, hiseq_type)]
                if isinstance(d, list):
                    out += d
            elif os.path.isfile(x):
                d = []
                with open(x) as r:
                    for line in r:
                        line = line.strip()
                        if line.startswith('#') or len(line) == 0:
                            continue
                        d.append(line.strip().split(' ')[0])
                out = [j for i in d for j in self.parse_dir(i)]
            else:
                pass
        elif isinstance(x, list):
            out = [j for i in x for j in self.parse_dir(i)]
        else:
            pass
        return out
    
    
    def init_files(self):
        # dirs
        self.project_dir = self.outdir
        default_dirs = {
            'config_dir': 'config',
            'data_dir': 'data',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        # files
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml', # updated
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
            'indir_json': self.data_dir + '/indir.json',
            # qc files
            'trim_summary_json':self.qc_dir +  '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'dup_summary_json': os.path.join(self.qc_dir, '01.pcr_dup_summary.json'),
            'lendist_csv': self.qc_dir + '/02.length_distribution.csv',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_json': self.qc_dir + '/03.FRiP.json',            
        }
        self = update_obj(self, default_files, force=True) # key
        # dirs
        dir_list = [
            self.config_dir, self.data_dir, self.qc_dir, self.report_dir,
        ]
        check_path(dir_list, create_dirs=True)

    
def get_args():
    """Parsing arguments for cnr_rx
    """
    example = '\n'.join([
        'Examples:',
        '1. Organize multiple hiseq dirs (cnr)',
        '$ python hiseq_merge.py -o ',
    ])
    parser = argparse.ArgumentParser(
        prog='hiseq_merge',
        description='hiseq_merge: merge multiple hiseq dirs',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--indir', nargs='+', required=True,
        help='hiseq project directories')
    parser.add_argument('-o', '--outdir', default=None,
        help='Directory saving results, default: [./hiseq_merge]')
    parser.add_argument('-f', '--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    CnrMerge(**args).run()


if __name__ == '__main__':
    main()

#
