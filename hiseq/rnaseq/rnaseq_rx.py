#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
RNAseq pipeline: level-3 (run rx)


loading fastq from: mut_fq1, mut_fq2, wt_fq1, wt_fq2
run pipeline, with specific parameters


analysis-module:
"""


import os
import sys
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.rnaseq.rnaseq_r1 import RnaseqR1
from hiseq.rnaseq.rnaseq_rn import RnaseqRn
from hiseq.rnaseq.rnaseq_rp import RnaseqRp
from hiseq.rnaseq.utils import rnaseq_trim, rnaseq_align_spikein, \
    rnaseq_align_rRNA, rnaseq_align_genome, rnaseq_quant, rnaseq_bam2bw, \
    rnaseq_deseq, qc_trim_summary, qc_align_summary, qc_bam_cor, \
    qc_genebody_enrich  # , salmon_align, salmon_deseq
from hiseq.rnaseq.rnaseq_salmon import RnaseqSalmon

from hiseq.utils.file import check_path, check_fx_paired, symlink_file, \
    file_exists, file_abspath, file_prefix, fx_name, Genome
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu, \
    read_hiseq, list_hiseq_file, run_shell_cmd
from hiseq.utils.bam import Bam
from hiseq.align.align_index import AlignIndex



class RnaseqRx(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RnaseqRxConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)

        
    def run_rn(self):
        args_required = [
            'aligner', 'outdir', 'gene_bed', 'gene_gtf', 
            'extra_index', 'genome', 'genome_index', 'spikein', 'spikein_index',
            'rRNA_index', 'to_rRNA', 'genome_size', 'genome_size_file',
            'threads', 'parallel_jobs', 'overwrite', 'binsize', 
            'trimmed', 'cut_to_length', 'recursive'
        ]
        args_global = dict((k, getattr(self, k)) for k in args_required \
            if k in self.__dict__)
        # for mut
        args_mut = args_global.copy()
        args_mut.update({
            'is_mut': True,
            'fq1': self.mut_fq1,
            'fq2': self.mut_fq2,
            'smp_name': self.mut_name,
            'parallel_jobs': 1,
        })
        RnaseqRn(**args_mut).run()
        # for wt
        args_wt = args_global.copy()
        if self.wt_fq1 is not None:
            args_wt.update({
                'is_mut': False,
                'fq1': self.wt_fq1,
                'fq2': self.wt_fq2,
                'smp_name': self.wt_name,
            })
        RnaseqRn(**args_wt).run()

    
    def copy_r1_files(self):
        # bam, bw, count
        for i in [self.mut_dir, self.wt_dir]:
            for j in ['r1', 'rn']:
                b = list_hiseq_file(i, 'bam', j)
                bw1 = list_hiseq_file(i, 'bw_fwd', j)
                bw2 = list_hiseq_file(i, 'bw_rev', j)
                c1 = list_hiseq_file(i, 'count_sens', j)
                c2 = list_hiseq_file(i, 'count_anti', j)
                [symlink_file(bam, self.bam_dir) for bam in b]
                [symlink_file(bw, self.bw_dir) for bw in bw1+bw2]
                [symlink_file(c, self.count_dir) for c in c1+c2]


    def run_rx(self):
        """
        Run mut over wt, for quality control
        """
        # DEseq
        rnaseq_deseq(self.project_dir, 'rx')
        # qc
        qc_genebody_enrich(self.project_dir, 'rx')
        qc_bam_cor(self.project_dir, 'rx')
        
        
    def run_salmon(self):
        """
        Run RNAseq with Salmon+DESeq2
        """
        args = {
            'mut_fq1': self.mut_fq1,
            'mut_fq2': self.mut_fq2,
            'wt_fq1': self.wt_fq1,
            'wt_fq2': self.wt_fq2,
            'extra_index': self.salmon_index,
            'outdir': self.salmon_dir,
            'threads': self.threads,
            'genome': self.genome,
        }
        if AlignIndex(self.salmon_index).is_valid():
            RnaseqSalmon(**args).run()
        else:
            log.warning('--salmon-index not valid, {}'.foramt(self.salmon_index))
        
        
    def run(self):
        # 1. save config
        Config().dump(self.__dict__, self.config_yaml)
        # 2. run RnaseqRn
        self.run_rn()
        # 3. run RnaseqRx
        self.copy_r1_files()
        self.run_rx()
        # 4. generate report
        RnaseqRp(self.project_dir).run()
        # 5. run salmon
        self.run_salmon()
        

class RnaseqRxConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'aligner': 'STAR',
            'mut_fq1': None,
            'mut_fq2': None,
            'wt_fq1': None,
            'wt_fq2': None,
            'mut_name': None,
            'wt_name': None,
            'smp_name': None, # mut.vs.wt
            'outdir': None,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': None,
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'to_rRNA': False,
            'rRNA_index': None,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'norm_project': None,
            'trimmed': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rx'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(os.path.expanduser(self.outdir))
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        self.mut_fq1, self.mut_fq2, self.mut_is_paired = self.init_fq(self.mut_fq1, self.mut_fq2)
        self.wt_fq1, self.wt_fq2, self.wt_is_paired = self.init_fq(self.wt_fq1, self.wt_fq2)
        # paired-end
        self.is_paired = self.mut_is_paired and self.wt_is_paired
        self.init_index()
        self.init_name()
        self.init_files()


    def init_fq(self, fq1, fq2):
        # convert to list
        if isinstance(fq1, str):
            fq1 = [fq1]
        if isinstance(fq2, str):
            fq2 = [fq2]
        if not isinstance(fq1, list):
            raise ValueError('fq1 require list, got {}'.format(
                type(fq1).__name__))
        # check message
        c1 = isinstance(fq1, list)
        c2 = isinstance(fq2, list)
        c1e = all(file_exists(fq1))
        c1x = all([c1, c1e])
        if c1:
            if isinstance(fq2, list):
                c2e = all(file_exists(fq2))
                c2p = all(check_fx_paired(fq1, fq2))
                c2x = all([c2, c2e, c2p])
            elif fq2 is None:
                fq2 = None
                c2x = c2p = c2e = True # skipped
            else:
                c2x = c2e = c2p = False # force                
        else:
            c1x = c2x = c2e = c2p = False # force
        # final
        c = all([c1x, c2x])
        if not c:
            msg = '\n'.join([
                '='*80,
                'Check fastq:',
                '{:>14} : {}'.format('fq1', fq1),
                '{:>14} : {}'.format('fq2', fq2),
                '-'*40,
                'Status',
                '{:>14} : {}'.format('fq1 is list', c1),
                '{:>14} : {}'.format('fq1 exists', c1e),
                '{:>14} : {}'.format('fq2 is list', c2),
                '{:>14} : {}'.format('fq2 is exists', c2e),
                '{:>14} : {}'.format('fq is paired', c2p),
                '-'*40,
                'Status: {}'.format(c),
                '='*80                
            ])
            print(msg)
            raise ValueError('fq1, fq2 not valid')
        return (file_abspath(fq1), file_abspath(fq2), c2p)
                
    
    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).get_fasize()
        else:
            raise ValueError('--genome or --extra-index; required')


    def init_name(self):
        """
        mut_name, wt_name, smp_name
        mut_dir, wt_dir
        """
        if not isinstance(self.mut_name, str):
            # fix mut_name
            self.mut_name = fx_name(self.mut_fq1[0], fix_pe=True, fix_rep=True)
        if not isinstance(self.wt_name, str):
            # fix wt_name
            self.wt_name = fx_name(self.wt_fq1[0], fix_pe=True, fix_rep=True)
        if not isinstance(self.smp_name, str):
            # fix smp_name
            self.smp_name = '{}.vs.{}'.format(self.mut_name, self.wt_name)
        # update mut/wt dirs
        self.mut_dir = os.path.join(self.outdir, self.mut_name)
        self.wt_dir = os.path.join(self.outdir, self.wt_name)
    
    
    def init_files(self):
        # dirs
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        default_dirs = {
            'config_dir': 'config',
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'count_dir': 'count',
            'deseq_dir': 'deseq',
            'enrich_dir': 'enrich',
            'qc_dir': 'qc',
            'report_dir': 'report',
            'salmon_dir': 'salmon'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key
        # files
        mut_bam_prefix = os.path.join(self.bam_dir, self.mut_name)
        mut_bw_prefix = os.path.join(self.bw_dir, self.mut_name)
        mut_count_prefix = os.path.join(self.count_dir, self.mut_name)
        wt_bam_prefix = os.path.join(self.bam_dir, self.wt_name)
        wt_bw_prefix = os.path.join(self.bw_dir, self.wt_name)
        wt_count_prefix = os.path.join(self.count_dir, self.wt_name)
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml',
            'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
            
            # mut files
            'mut_bam': mut_bam_prefix + '.bam',
            'mut_bw_fwd': mut_bw_prefix + '.fwd.bigWig',
            'mut_bw_rev': mut_bw_prefix + '.rev.bigWig',
            'mut_count_sens': mut_count_prefix + '.sens.txt',
            'mut_count_anti': mut_count_prefix + '.anti.txt',
            
            # wt files
            'wt_bam': wt_bam_prefix + '.bam',
            'wt_bw_fwd': wt_bw_prefix + '.fwd.bigWig',
            'wt_bw_rev': wt_bw_prefix + '.rev.bigWig',
            'wt_count_sens': wt_count_prefix + '.sens.txt',
            'wt_count_anti': wt_count_prefix + '.anti.txt',
            
            # qc files
            'genebody_enrich_matrix': os.path.join(self.qc_dir, '05.genebody_enrich.mat.gz'),
            'genebody_enrich_matrix_log': os.path.join(self.qc_dir, '05.genebody_enrich.log'),
            'genebody_enrich_png': os.path.join(self.qc_dir, '05.genebody_enrich.png'),
            'genebody_enrich_cmd': os.path.join(self.qc_dir, '05.genebody_enrich.cmd.sh'),
            'bam_cor_npz': os.path.join(self.qc_dir, '06.bam_cor.npz'),
            'bam_cor_counts': os.path.join(self.qc_dir, '06.bam_cor.counts.tab'),
            'bam_cor_heatmap_png': os.path.join(self.qc_dir, '06.bam_cor.cor_heatmap.png'),
            'bam_cor_pca_png': os.path.join(self.qc_dir, '06.bam_cor.cor_PCA.png'),
            
            # deseq
            'deseq_fix_xls': os.path.join(self.deseq_dir, 'transcripts_deseq2.fix.xls')
        }
        self = update_obj(self, default_files, force=True) # key
        # dirs
        dir_list = [
            self.config_dir, self.bam_dir, self.bw_dir, self.count_dir,
            self.deseq_dir, self.enrich_dir, self.qc_dir, self.report_dir
        ]
        check_path(dir_list, create_dirs=True)

         
def get_args():
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python rnaseq_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6',
        '2. for specific index',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te',
    ])    
    parser = argparse.ArgumentParser(
        prog='rnaseq_rx',
        description='run rnaseq_rx',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m1', '--mut-fq1', nargs='+', required=True,
        dest='mut_fq1',
        help='read1 files, (or read1 of PE reads) of treatment/mutant samples')
    parser.add_argument('-m2', '--mut-fq2', nargs='+', required=False,
        dest='mut_fq2',
        help='read2 files, (or read2 of PE reads) of treatment/mutant samples')
    parser.add_argument('-w1', '--wt-fq1', nargs='+', required=True,
        dest='wt_fq1',
        help='read1 files, (or read1 of PE reads) of control/wt samples')
    parser.add_argument('-w2', '--wt-fq2', nargs='+', required=False,
        dest='wt_fq2',
        help='read2 files, (or read2 of PE reads) of control/wt samples')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: hg38')
    parser.add_argument('--gtf', dest='gene_gtf', default=None,
        help='The gtf file for quantification, defaut: genome.gtf (None)')
    parser.add_argument('--bed', dest='gene_bed', default=None,
        help='The BED of genes')
    # optional arguments - 0
    parser.add_argument('-mn', '--mut-name', dest='mut_name', default=None,
        help='Name of mutant samples')
    parser.add_argument('-wn', '--wt-name', dest='wt_name', default=None,
        help='Name of wildtype/control samples')
    parser.add_argument('-n', '--smp-name', dest='smp_name', default=None,
        help='Name of the samples, default: from fq1 prefix')
    # optional arguments - 0
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    # optional arguments - 1
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    ## extra: index
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--genome-index', dest="genome_index", default=None,
        help='align index of genome')
    parser.add_argument('-k', '--spikein', default=None,
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--spikein-index', dest="spikein_index", default=None,
        help='align index of spikein')    
    parser.add_argument('--to-rRNA', dest='to_rRNA', action='store_true',
        help='Align to rRNA')
    parser.add_argument('--rRNA-index', dest="rRNA_index", default=None,
        help='align index of rRNA')   
    parser.add_argument('--aligner', default='STAR',
        choices=['STAR', 'bowtie', 'bowtie2', 'bwa', 'hisat2', 'kallisto',
            'salmon'],
        help='Aligner option: [STAR, bowtie, bowtie2, bwa], default: [STAR]')
    parser.add_argument('--salmon-index', dest='salmon_index', required=False,
        default=None,help='The path to salmon index')
    
    ## extra:
    parser.add_argument('-bs', '--bin-size', default=10, type=int,
        help='bin size of the bigWig file, default [10]')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. \
        default: [None]')
    parser.add_argument('--norm-project', dest='norm_project', default=None,
        help='The RNAseq_Rx project, for parseing norm scale. eg: \
        RNAseq_gene/wt.vs.mut for RNAseq_te, default: [None]')
    return parser


def main():
    args = vars(get_args().parse_args())
    RnaseqRn(**args).run()


if __name__ == '__main__':
    main()

#
