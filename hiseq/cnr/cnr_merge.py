#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Combine/merge multiple cnr report

cnr_r1, cnr_rn, ...
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.cnr.cnr_rn import CnrRn
from hiseq.cnr.cnr_rp import CnrRp
from hiseq.utils.file import (
    check_fx_args, check_path, symlink_file, file_abspath, file_prefix,
    file_exists, check_fx_paired, fx_name, list_file, list_dir, copy_file, copy_dir
)
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, init_cpu, print_dict, read_hiseq,
    find_longest_common_str, is_hiseq_dir, list_hiseq_file
)
from hiseq.cnr.utils import (
    hiseq_bam2bw,
    cnr_call_peak, cnr_merge_bam, qc_lendist, qc_frip, 
    qc_bam_cor, qc_peak_idr, qc_peak_overlap, qc_bam_fingerprint, 
    qc_tss_enrich, qc_genebody_enrich, cnr_bw_compare, 
    copy_hiseq_qc
)
from hiseq.utils.genome import Genome
from hiseq.align.align_index import AlignIndex, check_index_args



class CnrMerge(object):
    def __init__(self, **kwargs):
        c = CnrMergeConfig(**kwargs)
        self = update_obj(self, c.__dict__, force=True)
        
        
    def run(self):
        Config().dump(self.__dict__, self.config_yaml) # save config
        Config().dump({'indir': self.indir}, self.indir_json) # save dir list
        copy_hiseq_qc(self.project_dir) # copy files, config
        
        

#     def prepare_files(self):
#         # ip - bam,bw
#         ra = read_hiseq(self.ip_dir, 'rn')
#         symlink_file(ra.bam, self.ip_bam)
#         symlink_file(ra.bw, self.ip_bw)
#         # input - bam,bw
#         rb = read_hiseq(self.input_dir, 'rn')
#         if rb.is_hiseq:
#             symlink_file(rb.bam, self.input_bam)
#             symlink_file(rb.bw, self.input_bw)
    
    
#     def run_rx(self):
#         """
#         Run ip over input (IgG), for quality control
#         """
#         # ip-over-input
#         cnr_bw_compare(self.project_dir, 'rx') # generate bw, ip.over.input
#         cnr_call_peak(self.project_dir, 'rx') # call peak
#         # cnr_call_motifs(self.project_dir, 'rx') # to-do
#         # qc - enrich
#         qc_tss_enrich(self.project_dir, 'rx')
#         qc_genebody_enrich(self.project_dir, 'rx')
#         # qc - bam cor
#         qc_bam_cor(self.project_dir, 'rx')
#         # qc - fingerprint
#         qc_bam_fingerprint(self.project_dir, hiseq_type='rx', bam_type='rn')
#         # qc - motifs-fingerprint
        
    
#     def run_ip_only(self):
#         """
#         copy ip files to rx directory: bam,bw,qc
#         """
#         # ip - bam,bw
#         ra = read_hiseq(self.ip_dir, 'rn')
#         symlink_file(a.peak, self.peak)
#         symlink_file(a.peak_seacr, self.peak_seacr)
#         symlink_file(a.peak_seacr_top001, self.peak_seacr_top001)
#         # ip - qc
#         symlink_file(a.tss_enrich_png, self.tss_enrich_png)
#         symlink_file(a.genebody_enrich_png, self.genebody_enrich_png)
#         symlink_file(a.bam_fingerprint_png, self.bam_fingerprint_png)

            
#     def run(self):
#         # 1. save config
#         Config().dump(self.__dict__, self.config_yaml)
#         # 2. run CnrRn
#         self.run_rn()
#         # 3. run CnrRx
#         self.prepare_files()
#         if self.input_fq1 is None:
#             print('!ip - only, ...')
#             self.run_ip_only()
#         else:
#             self.run_rx()
#         # 4. generate report
#         CnrRp(self.project_dir).run()


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
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.indir = self.parse_dir(self.indir)
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
                with open(x) as r:
                    for line in r:
                        line = line.strip()
                        if line.startswith('#') or len(line) == 0:
                            continue
                        d = line.strip().split(' ')
                out = [parse_dir(i) for i in d]
            else:
                pass
        elif isinstance(x, list):
            out = [parse_dir(i) for i in x]
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
        '1. Run pipeline for design.yaml, with different parameters',
        '$ python cnr_rx.py -d design.yaml -g dm6 -o results',
        '2. Run pipeline with different parameters',
        '$ python cnr_rx.py -d design.yaml -g dm6 -o results --extra-index te',
    ])
    parser = argparse.ArgumentParser(
        prog='cnr_rx',
        description='cnr_rx: for multiple groups of PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config; generated by cnr_rd.py')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: dm6')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index (bowtie2)')
    parser.add_argument('--bed', '--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes, for TSS enrichment analysis')
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--cut', action='store_true', 
        help='Cut reads to 50nt, equal to: --cut-to-length 50 --recursive')

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
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')
    return parser


def main():
    args = vars(get_args().parse_args())
    CnrRx(**args).run()


if __name__ == '__main__':
    main()

#
