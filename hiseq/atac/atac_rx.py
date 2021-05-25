#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-2 (run pipe)

loading fastq config from `design.toml`, run pipeline, with specific parameters

analysis-module:
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.atac.atac_rn import AtacRn
from hiseq.atac.atac_rp import AtacRp
from hiseq.utils.file import check_path, symlink_file, file_abspath, file_prefix, file_exists
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu


class AtacRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = AtacRxConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def run_single_group(self, i):
        # update arguments
        i = {k:v for k,v in i.items() if k in ['fq1', 'fq2', 'group']}
        i.update({
            'build_design': False,
            'design': None,
        })
        if len(self.fq_groups) > 1:
            i['parallel_jobs'] = 1 # force
        args = self.__dict__.copy()
        args.update(i)
        AtacRn(**args).run()


    def run_multi_group(self):
        """
        self.fq_groups
          - key: string
          - value: fq1/fq2/rep_list
        """
        if self.parallel_jobs > 1 and len(self.fq_groups) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_group, self.fq_groups.values())
        else:
            for d in list(self.fq_groups.values()):
                self.run_single_group(d)


    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_toml)
        # 2. run AtacRn->AtacR1
        self.run_multi_group()
        # 3. generate report
        AtacRp(self.project_dir).run()


class AtacRxConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'design': None, # required
            'outdir': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'atacseq_rx'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        if not file_exists(self.design):
            raise ValueError('file not exists, --design {}'.format(self.design))
        self.init_files()
        self.init_fx()


    def init_files(self):
        # dirs
        self.project_dir = os.path.join(self.outdir)
        self.config_dir = os.path.join(self.project_dir, 'config')
        self.report_dir = os.path.join(self.project_dir, 'report')
        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml',
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
        }
        self = update_obj(self, default_files, force=True) # key
        check_path([self.config_dir, self.report_dir], create_dirs=True)


    def init_fx(self):
        """
        Loading fastq config from design
        """
        self.fq_groups = Config().load(self.design)
        if len(self.fq_groups) == 0:
            raise ValueError('no data in design: {}'.format(self.design))


def get_args():
    """Parsing arguments for atac_rx
    """
    example = '\n'.join([
        'Examples:',
        '1. Run pipeline for design.toml, with different parameters',
        '$ python atac_rx.py -d design.toml -g dm6 -o results',
        '2. Run pipeline with different parameters',
        '$ python atac_rx.py -d design.toml -g dm6 -o results --extra-index te',
    ])
    parser = argparse.ArgumentParser(
        prog='atac_rx',
        description='atac_rx: for multiple groups of PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config; generated by atac_rd.py')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: dm6')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index (bowtie2)')
    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes, for TSS enrichment analysis')

    # optional arguments - 1
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1, type=int,
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    parser.add_argument('--cut-to-length', dest='cut_to_length',
        default=0, type=int,
        help='cut reads to specific length from tail, default: [0]')
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')
    return parser


def main():
    args = vars(get_args().parse_args())
    AtacRx(**args).run()


if __name__ == '__main__':
    main()


#