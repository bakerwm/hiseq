#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
CnRseq pipeline: level-3 (run _rn)

fastq files from: fq1/fq2

analysis-module:
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.cnr.cnr_r1 import CnrR1
from hiseq.cnr.cnr_rp import CnrRp
from hiseq.cnr.utils import cnr_trim, cnr_align_genome, cnr_merge_bam, \
    cnr_align_spikein, cnr_call_peak, cnr_bam_to_bw, \
    qc_trim_summary, qc_align_summary, qc_lendist, qc_frip, \
    qc_bam_cor, qc_peak_idr, qc_tss_enrich, qc_genebody_enrich, \
     qc_peak_overlap, qc_bam_fingerprint
from hiseq.utils.file import check_path, check_fx_paired, symlink_file, \
    file_exists, file_abspath, file_prefix, fx_name, Genome, list_dir
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu, \
    read_hiseq, list_hiseq_file, run_shell_cmd
from hiseq.utils.bam import Bam
from hiseq.align.align_index import AlignIndex


class CnrRn(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        args_local = CnrRnConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        self.init_files()


    def init_files(self): # from rep_list
        self.bam_list = list_hiseq_file(self.project_dir, 'bam', 'r1')
        self.bw_list = list_hiseq_file(self.project_dir, 'bw', 'r1')
        self.peak_list = list_hiseq_file(self.project_dir, 'peak', 'r1')


    def run_pipe_rn(self): # for rep_list > 1
        cnr_merge_bam(self.project_dir, 'rn')
        cnr_call_peak(self.project_dir, 'rn')
        cnr_bam_to_bw(self.project_dir, 'rn')
        cnr_call_peak(self.project_dir, 'rn')
        # qc_trim(self.project_dir, 'rn')
        # qc_align(self.project_dir, 'rn')
        qc_lendist(self.project_dir, 'rn')
        qc_frip(self.project_dir, 'rn')
        qc_tss_enrich(self.project_dir, 'rn')
        qc_genebody_enrich(self.project_dir, 'rn')
        # specific for rn
        qc_bam_cor(self.project_dir, 'rn')
        qc_peak_idr(self.project_dir, 'rn')
        qc_peak_overlap(self.project_dir, 'rn')
        qc_bam_fingerprint(self.project_dir, 'rn')


    def run_pipe_r1(self): # for rep_list == 1
        log.warning('merge() skipped, Only 1 replicate detected')
        k_list = [
            'bam', 'bw', 'peak', 'peak_seacr', 'peak_seacr_top001',
            'align_scale_json', 'trim_json', 'align_json'
        ]
        # get the replist
        rep_dir = self.rep_list[0] # first one
        for k in k_list:
            k_from = list_hiseq_file(rep_dir, k, 'r1')
            # k_to = list_hiseq_file(self.project_dir, k, 'rn')
            k_to = getattr(self, k)
            symlink_file(k_from[0], k_to)
        # copy all files in qc dir
        rep_qc_dir = list_hiseq_file(rep_dir, 'qc_dir', 'r1')
        rep_qc_files = list_dir(rep_qc_dir[0], include_dir=True)
        for f in rep_qc_files:
            symlink_file(f, self.qc_dir) # to qc_dir
        # update: bam index
        Bam(self.bam).index()


    def run_single_fx(self, i):
        """
        Parameters
        ----------
        i:  int
            The index of fq1 files in self.fq1
        """
        # required args
        args_required = ['aligner', 'fq1', 'fq2', 'genome', 'gene_bed',
            'genome_size', 'is_ip', 'trimmed', 'outdir', 'overwrite',
            'parallel_jobs', 'threads', 'spikein', 'spikein_index',
            'extra_index', 'cut_to_length', 'recursive']
        args = dict((k, getattr(self, k)) for k in args_required \
            if k in self.__dict__)
        # update fq1, fq2, rep_list, ...
        args.update({
            'fq1': self.fq1[i],
            'fq2': self.fq2[i],
            'rep_list': None,
            'build_design': None,
            'design': None,
        })
        CnrR1(**args).run()

            
    def run_multi_fx(self):
        i_list = range(len(self.fq1))
        # in parallel
        if self.parallel_jobs > 1 and len(self.fq1) > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single_fx, i_list)
        else:
            for i in i_list:
                self.run_single_fx(i)


    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_yaml)
        # 2. run CnrR1
        self.run_multi_fx()
        # 3. run CnrRn, merge
        if len(self.rep_list) == 1:
            self.run_pipe_r1()
        else:
            self.run_pipe_rn()
        # 4. generate reprot
        CnrRp(**self.__dict__).run()


class CnrRnConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'bowtie2',
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index': None,
            'smp_name': None,
            'rmdup': True, # key
            'genome': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': 0,
            'keep_tmp': None,
            'trimmed': False,
            'cut': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_rn'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.init_cut()
        self.init_fx()
        self.init_files()
        self.init_index()
        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)


    def init_cut(self):
        """
        Cut the reads to specific length
        """
        if self.cut:
            self.cut_to_length = 50
            self.recursive = True
        
        
    def init_fx(self):
        """
        required:
        1. fq1:list, fq2:list
        2. paired
        """
        flag_err = True
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]
        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]
        if not all([
            isinstance(self.fq1, list),
            isinstance(self.fq2, list),
            check_fx_paired(self.fq1, self.fq2)
        ]):
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        # auto: sample names
        if self.smp_name is None:
            fname = fx_name(self.fq1, fix_pe=True, fix_rep=True)
            fname = list(set(fname))
            if len(fname) > 1:
                log.warning('fq1 name not identical')
            self.smp_name = fname[0] # first one
        # update rep_list
        self.rep_list = [
            os.path.join(self.outdir, i) for i in fx_name(self.fq1, fix_pe=True)
        ]


    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).get_fasize()
        else:
            raise ValueError('--genome or --extra-index; required')


    def init_files(self):
        # dirs
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')
        default_dirs = {
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
        trim_prefix = os.path.join(self.clean_dir, self.smp_name)
        align_prefix = os.path.join(self.align_dir, self.smp_name)
        spikein_prefix = os.path.join(self.spikein_dir, self.smp_name)
        peak_prefix = os.path.join(self.peak_dir, self.smp_name)
        default_files = {
            # basic files
            'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
            'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
            'bam_raw': self.bam_dir + '/' + self.project_name + '.raw.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'peak': peak_prefix+'_peaks.narrowPeak',
            'peak_seacr': peak_prefix+'.stringent.bed',
            'peak_seacr_top001': peak_prefix+'.top0.01.stringent.bed',
            
            # trimming
            'trim_stat': trim_prefix+'.trim.stat',
            'trim_json': trim_prefix+'.trim.json',
            
            # align files
            'align_scale_json': align_prefix+'.scale.json',
            'align_stat': align_prefix+'.align.stat',
            'align_json': align_prefix+'.align.json',
            'align_flagstat': align_prefix+'.align.flagstat',

            # spikein files
            'spikein_scale_json': spikein_prefix+'.scale.json',
            'spikein_stat': spikein_prefix+'.align.stat',
            'spikein_json': spikein_prefix+'.align.json',
            'spikein_flagstat': spikein_prefix+'.align.flagstat',
            'align_scale_json': self.bam_dir + '/' + 'scale.json',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.bowtie2.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.bowtie2.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',

            # qc files
            'trim_summary_json':self.qc_dir +  '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
            'lendist_csv': self.qc_dir + '/02.length_distribution.csv',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.pdf',
            'frip_json': self.qc_dir + '/03.FRiP.json',
            'tss_enrich_matrix': self.qc_dir + '/04.tss_enrich.mat.gz',
            'tss_enrich_matrix_log': self.qc_dir + '/04.tss_enrich.log',
            'tss_enrich_png': self.qc_dir + '/04.tss_enrich.png',
            'tss_enrich_cmd': self.qc_dir + '/04.tss_enrich.cmd.sh',
            'genebody_enrich_matrix': self.qc_dir + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png',
            'peak_idr_png': self.qc_dir + '/07.peak_idr.png',
            'peak_idr_txt': self.qc_dir + '/07.peak_idr.txt',
            'peak_overlap_png': self.qc_dir + '/08.peak_overlap.png',
            'peak_overlap_tiff': self.qc_dir + '/08.peak_overlap.tiff',
            'bam_fingerprint_png': self.qc_dir + '/09.fingerprint.png'
        }
        self = update_obj(self, default_files, force=True) # key
        dir_list = [
            self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
            self.spikein_dir, self.bam_dir, self.bg_dir, self.bw_dir,
            self.peak_dir, self.motif_dir, self.qc_dir, self.report_dir,
        ]
        check_path(dir_list, create_dirs=True)


def get_args():
    """Parsing arguments for cnr: rn
    """
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python cnr_rn.py --fq1 *1.fq.gz --fq2 *2.fq.gz -o results -g dm6',
        '2. for specific index',
        '$ python cnr_rn.py --fq1 *1.fq.gz --fq2 *2.fq.gz -o results -x bowtie2_index/te',
    ])
    parser = argparse.ArgumentParser(
        prog='cnr_rn',
        description='cnr_rn: for multiple PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='Fastq file, read1 of PE')
    parser.add_argument('-2', '--fq2', nargs='+', required=True,
        help='Fastq file, read2 of PE')
    parser.add_argument('-o', '--outdir', default=None,
        help='Directory saving results, default: [pwd]')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: dm6')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index (bowtie2)')
    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes, for TSS enrichment analysis')
    parser.add_argument('--is-ip', dest='is_ip', action='store_true',
        help='Is the IP sample')
    parser.add_argument('--trimmed', action='store_true',
        help='Skip trimming, input reads are already trimmed')
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
    CnrRn(**args).run()


if __name__ == '__main__':
    main()

#