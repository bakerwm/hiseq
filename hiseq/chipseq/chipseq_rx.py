#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Chipseq pipeline: level-2 (: rx; ip vs IgG)

# or without IgG => ATACseq pipe
loading fastq config from `design.yaml`, run pipeline, with specific parameters

analysis-module:
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.chipseq.chipseq_rn import ChipseqRn
from hiseq.chipseq.chipseq_rp import ChipseqRp
from hiseq.utils.file import (
    check_fx_args, check_fx_paired, check_path, symlink_file, file_abspath,
    file_prefix, file_exists, fx_name
)
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, init_cpu, print_dict, read_hiseq,
    find_longest_common_str
)
from hiseq.chipseq.utils import (
    hiseq_bam2bw, 
    chipseq_call_peak, chipseq_merge_bam, qc_lendist, qc_frip, qc_bam_cor,
    qc_peak_idr, qc_peak_overlap, qc_bam_fingerprint, qc_tss_enrich, 
    qc_genebody_enrich, chipseq_bw_compare
)
from hiseq.align.align_index import AlignIndex, check_index_args
from hiseq.utils.genome import Genome


class ChipseqRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = ChipseqRxConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)

        
    def run_rn(self):
        rq_list = [
            'outdir', 'index', 'genome', 'genome_index', 'extra_index', 
            'spikein', 'spikein_index', 'threads', 'parallel_jobs',
            'overwrite', 'binsize', 'genome_size', 'genome_size_file',
            'keep_tmp', 'trimmed', 'cut_to_length', 'recursive',
            'gene_bed'
        ]
        args_global = self.__dict__
        args_local = {k:v for k,v in args_global.items() if k in rq_list}
        # for ip
        args_local.update({
            'is_ip': True,
            'fq1': self.ip_fq1,
            'fq2': self.ip_fq2,
            'smp_name': self.ip_name,
            'parallel_jobs': 1,
        })
        ChipseqRn(**args_local).run()
        # for input
        if self.input_fq1 is not None:
            args_local.update({
                'is_ip': False,
                'fq1': self.input_fq1,
                'fq2': self.input_fq2,
                'smp_name': self.input_name,
            })
            ChipseqRn(**args_local).run()


    def prepare_files(self):
        # ip - bam,bw
        ra = read_hiseq(self.ip_dir, 'rn')
        symlink_file(ra.bam, self.ip_bam)
        symlink_file(ra.bw, self.ip_bw)
        # input - bam,bw
        rb = read_hiseq(self.input_dir, 'rn')
        if rb.is_hiseq:
            symlink_file(rb.bam, self.input_bam)
            symlink_file(rb.bw, self.input_bw)
    
    
    def run_rx(self):
        """
        Run ip over input (IgG), for quality control
        """
        # ip-over-input
        chipseq_bw_compare(self.project_dir, 'rx') # generate bw, ip.over.input
        chipseq_call_peak(self.project_dir, 'rx') # call peak
        # chipseq_call_motifs(self.project_dir, 'rx') # to-do
        # qc - enrich
        qc_tss_enrich(self.project_dir, 'rx')
        qc_genebody_enrich(self.project_dir, 'rx')
        # qc - bam cor
        qc_bam_cor(self.project_dir, 'rx')
        # qc - fingerprint
        qc_bam_fingerprint(self.project_dir, hiseq_type='rx', bam_type='rn')
        # qc - motifs-fingerprint
        
    
    def run_ip_only(self):
        """
        copy ip files to rx directory: bam,bw,qc
        """
        # ip - bam,bw
        ra = read_hiseq(self.ip_dir, 'rn')
        symlink_file(a.peak, self.peak)
        symlink_file(a.peak_seacr, self.peak_seacr)
        symlink_file(a.peak_seacr_top001, self.peak_seacr_top001)
        # ip - qc
        symlink_file(a.tss_enrich_png, self.tss_enrich_png)
        symlink_file(a.genebody_enrich_png, self.genebody_enrich_png)
        symlink_file(a.bam_fingerprint_png, self.bam_fingerprint_png)

            
    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_yaml)
        # 2. run ChipseqRn
        self.run_rn()
        # 3. run ChipseqRx
        self.prepare_files()
        if self.input_fq1 is None:
            print('!ip - only, ...')
            self.run_ip_only()
        else:
            self.run_rx()
        # 4. generate report
        ChipseqRp(self.project_dir).run()


class ChipseqRxConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
            'ip_fq1': None,
            'ip_fq2': None,
            'input_fq1': None,
            'input_fq2': None,
            'ip_name': None,
            'input_name': None,
            'smp_name': None, #ip.vs.input
            'outdir': None,
            'index': None,
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 100,
            'genome_size': 0,
            'genome_size_file': 0,
            'rmdup': True, # key
            'keep_tmp': None,
            'trimmed': False,
            'cut': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'chipseq_rx'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        self.init_cut()
        self.init_fq()
        self.init_name()
        self.init_index()
        self.init_files()


    def init_cut(self):
        """
        Cut the reads to specific length
        """
        if self.cut:
            self.cut_to_length = 50
            self.recursive = True
        

    def init_fq(self):
        """
        check fastq files: ip_fq1, ip_fq2, input_fq1, input_fq2
        """
        # ip
        if not check_fx_args(self.ip_fq1, self.ip_fq2):
            raise ValueError('ip fastq files failed')
        # input
        if not check_fx_args(self.input_fq1, self.input_fq2):
            raise ValueError('input fastq files failed')
        
        
    def init_name(self):
        """
        ip_name, input_name, smp_name
        ip_dir, input_dir
        """
        if not isinstance(self.ip_name, str):
            self.ip_name = fx_name(self.ip_fq1[0], fix_pe=True, fix_rep=True)
        if not isinstance(self.input_name, str):
            self.input_name = fx_name(self.input_fq1[0], fix_pe=True, fix_rep=True)
        # for ip.vs.input
        if not isinstance(self.smp_name, str):
            self.smp_name = '{}.vs.{}'.format(self.ip_name, self.input_name)
        # update ip/input dirs
        self.ip_dir = os.path.join(self.outdir, self.ip_name)
        self.input_dir = os.path.join(self.outdir, self.input_name)

    
    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).fasize()
        else:
            raise ValueError('--genome or --extra-index; required')


    def init_files(self):
        # dirs
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        default_dirs = {
            'config_dir': 'config',
            'raw_dir': 'raw_data',
            'clean_dir': 'clean_data',
            'align_dir': 'align',
            'spikein_dir': 'spikein',
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'peak_dir': 'peak',
            'motif_dir': 'motif',
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
            # symlinks - bam
            'ip_bam': self.bam_dir + '/' + self.ip_name + '.bam',
            'input_bam': self.bam_dir + '/' + self.input_name + '.bam',
            # symlinks - bigWig
            'ip_bw': self.bw_dir + '/' + self.ip_name + '.bigWig',
            'input_bw': self.bw_dir + '/' + self.input_name + '.bigWig',
            # ip over input
            'bam_raw': self.bam_dir + '/' + self.smp_name + '.raw.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            #'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',
            # qc - tss
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            # qc - genebody
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            # qc - bam:cor
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            # qc - idr
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_idr_txt': self.qc_dir + '/07.peak_idr.txt',
            # qc - peak:overlap
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'peak_overlap_tiff': self.qc_dir + '/08.peak_overlap.tiff',
            # qc - fingerprint
            'bam_fingerprint_png': self.qc_dir + '/09.fingerprint.png',
            # qc - motifs:fingerprint
        }
        self = update_obj(self, default_files, force=True) # key
        # dirs
        dir_list = [
            self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
            self.spikein_dir, self.bam_dir, self.bg_dir, self.bw_dir,
            self.peak_dir, self.motif_dir, self.qc_dir, self.report_dir,
        ]
        check_path(dir_list, create_dirs=True)

    
def get_args():
    """Parsing arguments for chipseq_rx
    """
    example = '\n'.join([
        'Examples:',
        '1. Run pipeline for design.yaml, with different parameters',
        '$ python chipseq_rx.py -d design.yaml -g dm6 -o results',
        '2. Run pipeline with different parameters',
        '$ python chipseq_rx.py -d design.yaml -g dm6 -o results --extra-index te',
    ])
    parser = argparse.ArgumentParser(
        prog='chipseq_rx',
        description='chipseq_rx: for multiple groups of PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config; generated by chipseq_rd.py')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: dm6')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index (bowtie2)')
    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
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
    ChipseqRx(**args).run()


if __name__ == '__main__':
    main()

#