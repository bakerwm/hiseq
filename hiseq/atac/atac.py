#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
ATAC-seq pipeline: level-1 (main port)

mission-1: generate design.toml

mission-2: run_pipe, parsing config from design.toml
"""

import os
import pathlib
import argparse
from hiseq.atac.atac_rd import AtacRd
from hiseq.atac.atac_rx import AtacRx
from hiseq.utils.utils import log, update_obj, Config, get_date
from hiseq.utils.file import file_abspath


class Atac(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'build_design': False,
            'design': None,
            'fq_dir': None,
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'genome': None,
            'gene_bed': None,
        }
        self = update_obj(self, args_init, force=False)


    def run(self):
        if self.build_design:
            AtacRd(**self.__dict__).run()
        else:
            AtacRx(**self.__dict__).run() # Rx->Rn->R1


def get_args():
    """Parsing arguments for atac
    """
    example = '\n'.join([
        'Examples:',
        '1. Generate design.toml for fastq files',
        '$ python atac.py -b -d design.toml',
        '2. Run pipeline for design.toml, with different parameters',
        '$ python atac.py -d design.toml -g dm6 -o results',
    ])
    parser = argparse.ArgumentParser(
        prog='atac',
        description='atac: ATACseq analysis pipeline',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='generate design.toml, with --fq-dir, or fq1/fq2')
    parser.add_argument('-d', '--design', required=True,
        help='The file saving fastq files config. eg: design.toml')
    parser.add_argument('-r', '--fq-dir', dest='fq_dir',
        help='Path to the directory, contains fastq files, eg: _rep1_1.fq.gz')
    parser.add_argument('-1', '--fq1', nargs='+', required=False,
        help='Fastq file, read1 of PE')
    parser.add_argument('-2', '--fq2', nargs='+', required=False,
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
    Atac(**args).run()


if __name__ == '__main__':
    main()

#




# to-be-removed
"""
1. cutadapt, recursive trimming

2. Alignment
--very-sensitive -X 2000

3. Mito%

pysam.idxstats(), chrM

4. Remove PCR duplicates

$ sambamba markdup

5. unique alignments

$ samtools view -q 30

6. proper mapped

$ samtools view -F 1804

0x4 Read unmapped
0x8 Mate unmapped
0x100 Not primary alignment
0x200 Read fails platform/vendor quality checks
0x400 Read is PCR or optical duplicate

7. merge replicate

$ samtools merge

8. Shifting reads (For TF footprinting)
For Tn5 insertion: +4, -5

$ samtools sort -n
$ bedtools bamtobed -i in.bam -bedpe | awk -v OFS="\t" '{($9=="+"){print $1,$2+4,$6+4} \
  ($9=="-"){print $1,$2-5,$6-5}}' > fragments.bed


9. Peak calling (MACS2)

$ macs2 callpeak --nomodel -f BAMPE --keep-dup all --cutoff-analysis
$ macs2 callpeak

10. bigWig

$ bamCoverage --binSize 10 --normalizeUsing RPGC


# use --ATACshift
alignmentSieve --numberOfProcessors 8 --ATACshift --bam sample1.bam -o sample1.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ 8 -O bam -o sample1.shifted.bam sample1.tmp.bam
samtools index -@ 8 sample1.shifted.bam
rm sample1.tmp.bam

Working mode:
1. ATACseqSingle(): for single sample/fastq
2. ATACseqMultiple(): for multiple sample/fastq files
3. ATACseqSample(): for single sample (could have multiple replicates)
4. ATACseq():

bam -> rmdup -> proper_paired -> peak/bw/...

##
AtacSingle()
AtacMultiple()
AtacReplicate()
Atac()

Input config file:

genome    outdir    name    fq1    fq2

Output:
- raw_data
- clean_data
- align
- bam_files
- bw_files
- peak
- motif
- report
  - reads num
  - peaks num
  - peaks annotation
  - qc

####

compare replicates:
- bam_files
- bw_files
- peak
- motif
- report

Mission:
atac_r1: single
atac_rn: merge replicates
atac_rx: multiple groups
"""




# class AtacConfig(object):
#     def __init__(self, **kwargs):
#         self = update_obj(self, kwargs, force=True)
#         self.init_args()


#     def init_args(self):
#         args_init = {
#             'aligner': 'bowtie2',
#             'build_design': False,
#             'design': None,
#             'fq_dir': None,
#             'rep_list': None,
#             'fq1': None,
#             'fq2': None,
#             'outdir': None,
#             'genome': None,
#             'genome_index': None,
#             'extra_index': None,
#             'spikein': None,
#             'spikein_index': None,
#             'threads': 1,
#             'parallel_jobs': 1,
#             'overwrite': False,
#             'binsize': 10,
#             'genome_size': 0,
#             'genome_size_file': 0,
#             'gene_bed': None,
#             'keep_tmp': None,
#             'trimmed': False,
#             'cut_to_length': 0,
#             'recursive': False
#         }
#         self = update_obj(self, args_init, force=False)
#         # outdir
#         if self.outdir is None:
#             self.outdir = str(pathlib.Path.cwd())
#         self.outdir = file_abspath(self.outdir)
#         # aligner
#         if self.aligner is None:
#             self.aligner = 'bowtie2'

