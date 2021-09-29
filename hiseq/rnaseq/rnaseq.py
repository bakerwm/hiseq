#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
RNAseq pipeline: level-1 (run all)

loading fastq from: design.yaml
run pipeline, with specific parameters

mission-1: generate design.yaml

mission-2: run_pipe, parsing config from design.yaml

Run Salmon alone:
outdir/salmon
  - config
  - index 
  - quant
  - deseq
  - report
  - ...
"""

import os
import pathlib
import argparse
from multiprocessing import Pool
from hiseq.rnaseq.rnaseq_rd import RnaseqRd
from hiseq.rnaseq.rnaseq_rx import RnaseqRx
from hiseq.utils.file import check_path, symlink_file, file_abspath
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu


class Rnaseq(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RnaseqConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def run_rd(self):
        RnaseqRd(**self.__dict__).run()
        

    def run_single_rx(self, i):
        i = {
            k:v for k,v in i.items() if k in \
            ['mut_fq1', 'mut_fq2', 'wt_fq1', 'wt_fq2']
        }
        i.update({
        'build_design': False,
        'design': None,
        })
#         if len(self.fq_groups) > 1:
#             i['parallel_jobs'] = 1 # force
        args_local = self.__dict__.copy()
        args_local.update(i)
        RnaseqRx(**args_local).run()


    def run_multiple_rx(self): # !!!! rx -> parallel_jobs=1
        # load fq groups
        self.fq_groups = Config().load(self.design)
        if len(self.fq_groups) == 0:
            raise ValueError('no data in design: {}'.format(self.design))
        # run multiple in parallel
        for i in list(self.fq_groups.values()):
            self.run_single_rx(i)
#         if self.parallel_jobs > 1 and len(self.fq_groups) > 1:
#             with Pool(processes=self.parallel_jobs) as pool:
#                 pool.map(self.run_single_rx, self.fq_groups.values())
#         else:
#             for i in list(self.fq_groups.values()):
#                 self.run_single_rx(i)
        
        
    def run(self):
        if self.build_design:
            self.run_rd()
        else:
            self.run_multiple_rx()
            

class RnaseqConfig(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'build_design': False,
            'design': None, # required
            'fq_dir': None,
            'mut': None,
            'wt': None,
            'mut_fq1': None,
            'mut_fq2': None,
            'wt_fq1': None,
            'wt_fq2': None,
            'threads': 4,
            'parallel_jobs': 1,
            'overwrite': False,
            'verbose': False,
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_ra'
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)
        if not isinstance(self.design, str):
            raise ValueError('--design, expect str, got {}'.format(
                type(self.design).__name__))


def get_args():
    example = '\n'.join([
        'Examples:',
        '1. Generate design.yaml， with --fq-dir',
        '$ python rnaseq.py -b -d design.yaml -r data --mut mutant --wt wildtype',
        '2. Generate design.yaml, with -1 and -2',
        '$ python rnaseq_rd.py -d design.yaml -1 *1.fq.gz -2 *2.fq.gz',
        '3. run pipeline',
        '$ python rnaseq.py -d design.yaml -o results -g dm6 --gtf g.gtf --bed g.bed -p 8',
    ])    
    parser = argparse.ArgumentParser(
        prog='rnaseq',
        description='RNAseq pipeline',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='generate design.yaml, with --fq-dir, or fq1/fq2')
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config; generated by rnaseq_rd.py')
    # for design
    parser.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('--mut', nargs='+', default=None,
        help='keyword of mut fastq file, auto-find read1/2')
    parser.add_argument('--wt', nargs='+', dest='wt', default=None,
        help='keyword of wt fastq file, auto-find read1/2')
    parser.add_argument('--se', dest='as_se', action='store_true',
        help='choose only fastq1 for PE reads')
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
    
    # for pipeline
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        # choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: hg38')
    parser.add_argument('--gtf', '--gene-gtf', dest='gene_gtf', default=None,
        help='The gtf file for quantification, defaut: genome.gtf (None)')
    parser.add_argument('--bed', '--gene-bed', dest='gene_bed', default=None,
        help='The BED of genes')

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
        # choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
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
    
    ## extra:    
    parser.add_argument('--salmon-index', dest='salmon_index', required=False,
        default=None,help='The path to salmon index')
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
    Rnaseq(**args).run()


if __name__ == '__main__':
    main()

#
