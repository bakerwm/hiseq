#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Generate argument parser for rnaseq pipeline

1. build design.yaml
2. run pipeline
3. arguments in details
"""

import os
import pathlib
import argparse
from hiseq.utils.utils import log, update_obj


class RnaseqArgs(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_parser() # main, sub_parsers

            
    def init_parser(self):
        example_main = '\n'.join([
            'Examples:',
            '$ hiseq rnaseq build -d design.json -r data --wt ctrl --mut treatment',
            '$ hiseq rnaseq run -d design.json -o results -g dm6'
        ])
        example_build = '\n'.join([
            'Examples:',
            '$ hiseq rnaseq build -d design.json -r data --wt ctrl --mut treatment'
        ])
        example_run = '\n'.join([
            'Examples:',
            '$ hiseq rnaseq run -d design.json -o results -g dm6'
        ])
        example_salmon = '\n'.join([
            'Examples:',
            '$ hiseq rnaseq salmon -d design.json -o results -g dm6'
        ])
        parser = argparse.ArgumentParser(
            prog='rnaseq',
            description='RNAseq pipeline',
            epilog=example_main,
            formatter_class=argparse.RawTextHelpFormatter
        )        
        subparsers = parser.add_subparsers(title='sub-commands', 
                                           help='choose one command')
        self.build = subparsers.add_parser('build', epilog=example_build,
                                          formatter_class=argparse.RawTextHelpFormatter)
        self.run = subparsers.add_parser('run', epilog=example_run,
                                        formatter_class=argparse.RawTextHelpFormatter)
        self.salmon = subparsers.add_parser('salmon', epilog=example_salmon,
                                           formatter_class=argparse.RawTextHelpFormatter)
        self.parser = parser # main
        # return (parser, build, run, salmon)

    
    def add_build(self, p):
        #-- arguments for build 
        p.add_argument('-d', '--design', required=True,
            help='The file saving fastq files config.yaml')
        p.add_argument('-r', '--fq-dir', dest='fq_dir', default=None,
            help='directory of fastq files, for --design')
        p.add_argument('--mut', nargs='+', default=None,
            help='keyword of mut fastq file, auto-find read1/2')
        p.add_argument('--wt', nargs='+', dest='wt', default=None,
            help='keyword of wt fastq file, auto-find read1/2')
        p.add_argument('--se', dest='as_se', action='store_true',
            help='choose only fastq1 for PE reads')
        return p


    def add_main(self, p):
        """
        Add arguments for 'run' and 'salmon'
        """
        p.add_argument('-d', '--design', required=True,
            help='The file saving fastq files config.yaml')
        p.add_argument('-o', '--outdir', default=None,
            help='The directory to save results, default, [pwd]')
        p.add_argument('-g', '--genome', default=None,
            help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: [None]')
        p.add_argument('-p', '--threads', default=1, type=int,
            help='Number of threads for each job, default [1]')
        p.add_argument('-j', '--parallel-jobs', default=1, type=int,
            dest='parallel_jobs',
            help='Number of jobs run in parallel, default: [1]')
        p.add_argument('--trimmed', action='store_true',
            help='specify if input files are trimmed')
        return p
        

    def add_extra(self, p):
        """
        Extra arguments for rnaseq
        """
        p.add_argument('-x', '--extra-index', dest="extra_index",
            help='Provide alignment index(es) for alignment, support multiple\
            indexes. if specified, ignore -g, -k')
        p.add_argument('--genome-index', dest="genome_index", default=None,
            help='align index of genome')
        p.add_argument('-k', '--spikein', default=None,
            # choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
            help='Spike-in genome : dm3, hg19, hg38, mm10, default: [None]')
        p.add_argument('--spikein-index', dest="spikein_index", default=None,
            help='align index of spikein')    
        p.add_argument('--to-rRNA', dest='to_rRNA', action='store_true',
            help='Align to rRNA')
        p.add_argument('--rRNA-index', dest="rRNA_index", default=None,
            help='align index of rRNA')
        p.add_argument('-bs', '--bin-size', default=10, type=int,
            help='bin size of the bigWig file, default [10]')
        p.add_argument('--overwrite', action='store_true',
            help='if spcified, overwrite exists file')
        p.add_argument('--gtf', '--gene-gtf', dest='gene_gtf', default=None,
            help='The gtf file for quantification, defaut: [auto]')
        p.add_argument('--bed', '--gene-bed', dest='gene_bed', default=None,
            help='The BED of genes, default: [auto]')
        p.add_argument('--aligner', default='STAR',
            choices=['STAR', 'bowtie', 'bowtie2', 'bwa', 'hisat2', 'kallisto',
                'salmon'],
            help='Aligner option: [STAR, bowtie, bowtie2, bwa], default: [STAR]')
        p.add_argument('-mn', '--mut-name', dest='mut_name', default=None,
            help='Name of mutant samples')
        p.add_argument('-wn', '--wt-name', dest='wt_name', default=None,
            help='Name of wildtype/control samples')
        return p
    
    
    def add_extra2(self, p):
        p.add_argument('-x', dest='salmon_index', required=False,
            help='The path to salmon index')
        return p
    
    
    def get_args(self):
        self.build = self.add_build(self.build)
        self.run = self.add_main(self.run)
        self.run = self.add_extra(self.run)
        self.salmon = self.add_main(self.salmon)
        self.salmon = self.add_extra2(self.salmon)
        return self.parser

    
def main():
    p = RnaseqArgs().get_args()
    p.parse_args()


if __name__ == '__main__':
    main()

#