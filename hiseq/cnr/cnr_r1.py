#!/usr/bin/env python

"""
Run CnR for single fastq file (PE or SE)

main:

1. trimming
double-trimming
a. full length adapter 
b. <6bp adapters

2. Alignment
dovetail alignment using bowtie2

--dovetail 

--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700

## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt

optional: 
spike-in
chrM

3. remove duplicates
remove or not

picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam 
O=$projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam 
REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt


4. peak calling
SEACR
Macs2

5. cut metrix


Quality control

1. fragment size (expected < 120bp)
2. adapter content pct (10-15%)
3. read duplication (10-15%)
4. alignment pct (>90%)
5. peaks ()
6. enrichment of motif 
7. 

"""


import os
import sys
from hiseq.cnr.cnr_rp import CnrRp
from hiseq.cnr.utils import cnr_trim, cnr_align_genome, \
    cnr_align_spikein, cnr_call_peak, cal_norm_scale, cnr_bam_to_bw, \
    qc_trim, qc_align, qc_lendist, qc_frip, qc_tss_enrich, qc_genebody_enrich
from hiseq.utils.file import check_path, check_fx_paired, symlink_file, \
    file_abspath, file_prefix, fx_name, Genome
from hiseq.utils.utils import log, update_obj, Config, get_date, init_cpu, \
    read_hiseq
from hiseq.align.align_index import AlignIndex


class CnrR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = CnrR1Config(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_toml)


    def run_pipe(self):
        # symlink raw data
        raw_fq1, raw_fq2 = self.raw_fq_list
        symlink_file(self.fq1, raw_fq1, absolute_path=True)
        symlink_file(self.fq2, raw_fq2, absolute_path=True)
        cnr_trim(self.project_dir, '_r1')
        cnr_align_genome(self.project_dir, '_r1')
        cnr_call_peak(self.project_dir, '_r1')
        cnr_bam_to_bw(self.project_dir, '_r1')
        if isinstance(self.spikein_index, str):
            cnr_align_spikein(self.project_dir, '_r1')
        # qc
        qc_trim(self.project_dir, '_r1')
        qc_align(self.project_dir, '_r1')
        qc_lendist(self.project_dir, '_r1')
        qc_frip(self.project_dir, '_r1')
        qc_tss_enrich(self.project_dir, '_r1')
        qc_genebody_enrich(self.project_dir, '_r1')


    def run(self):
        # 1. run pipe: CnrR1
        self.run_pipe()
        # 2. generate report
        CnrRp(**self.__dict__).run()


class CnrR1Config(object):
    """
    Prepare directories for R1 single replicates

    require: fq1/fq2, outdir, ...
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'is_ip': True,
            'aligner': 'bowtie2',
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'index': None,
            'smp_name': None,
            'rmdup': True,
            'genome': None,
            'gene_bed': None,
            'genome_index': None,
            'extra_index': None,
            'spikein': None,
            'spikein_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 10,
            'genome_size': 0,
            'genome_size_file': None,
            'keep_tmp': None,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_r1'
        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.init_fx()
        self.init_files()
        self.init_index()
        # threads
        self.threads, self.parallel_jobs = init_cpu(
            self.threads,
            self.parallel_jobs)


    def init_fx(self):
        """
        required:
        1. fq1:str, fq2:str
        2. paired
        """
        # check message
        c1 = isinstance(self.fq1, str)
        c2 = isinstance(self.fq2, str)
        c3 = check_fx_paired(self.fq1, self.fq2)
        if not all([c1, c2, c3]):
            # get message
            msg = '\n'.join([
                '='*80,
                'Input',
                '{:>14} : {}'.format('fq1', self.fq1),
                '{:>14} : {}'.format('fq2', self.fq2),
                '-'*40,
                'Status',
                '{:>14} : {}'.format('fq1 is str', c1),
                '{:>14} : {}'.format('fq1 is str', c2),
                '{:>14} : {}'.format('fq1,fq2 paired', c3),
                '-'*40,
                'Output: {}'.format(all([c1, c2, c3])),
                '='*80                
            ])
            print(msg)
            raise ValueError('fq1, fq2 not valid')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        # auto: sample names
        if self.smp_name is None:
            self.smp_name = fx_name(self.fq1, fix_pe=True)


    # update: genome_size_file    
    def init_index(self):
        # get data from: genome, extra_index
        if isinstance(self.extra_index, str):
            self.genome_size_file = AlignIndex(self.extra_index).index_size(out_file=True)
        elif isinstance(self.genome, str):
            self.genome_size_file = Genome(self.genome).get_fasize()
        else:
            raise ValueError('--genome or --extra-index; required')


    def init_files(self, create_dirs=True):
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
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True)
        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml', # updated
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
            'raw_fq1': self.raw_dir + '/' + os.path.basename(self.fq1),
            'raw_fq2': self.raw_dir + '/' + os.path.basename(self.fq2),
            'clean_fq1': self.clean_dir + '/' + os.path.basename(self.fq1),
            'clean_fq2': self.clean_dir + '/' + os.path.basename(self.fq2),

            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bam_rmdup': self.bam_dir + '/' + self.smp_name + '.rmdup.bam',
            'bam_proper_pair': self.bam_dir + '/' + self.smp_name + '.proper_pair.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'peak': self.peak_dir + '/' + self.project_name + '_peaks.narrowPeak',
            'peak_seacr': self.peak_dir + '/' + self.project_name + '.stringent.bed',
            'peak_seacr_top001': self.peak_dir + '/' + self.project_name + '.top0.01.stringent.bed',

            'align_scale_json': self.bam_dir + '/scale.json',
            'trim_stat_txt': self.clean_dir + '/' + self.project_name + '/' + self.project_name + '.trim.stat',
            'trim_stat_json': self.clean_dir + '/' + self.project_name + '/' + self.project_name + '.trim.json',
            'align_flagstat': self.align_dir + '/' + self.smp_name + '.flagstat',
            'align_stat': self.align_dir + '/' + self.smp_name + '.align.stat',
            'align_json': self.align_dir + '/' + self.smp_name + '.align.json',
            'spikein_bam': self.spikein_dir + '/spikein.bam',
            'spikein_scale_json': self.spikein_dir + '/scale.json',
            'spikein_flagstat': self.spikein_dir + '/spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/spikein.align.stat',
            'spikein_json': self.spikein_dir + '/spikein.align.json',

            'trim_summary_json': self.qc_dir + '/00.trim_summary.json',
            'align_summary_json': self.qc_dir + '/01.alignment_summary.json',
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
        }
        self = update_obj(self, default_files, force=True) # key
        self.raw_fq_list = [self.raw_fq1, self.raw_fq2]
        self.clean_fq_list = [self.clean_fq1, self.clean_fq2]
        # dirs
        dir_list = [
            self.config_dir, self.raw_dir, self.clean_dir, self.align_dir,
            self.spikein_dir, self.bam_dir, self.bg_dir, self.bw_dir,
            self.peak_dir, self.motif_dir, self.qc_dir, self.report_dir,
        ]
        check_path(dir_list, create_dirs=True)


def get_args():
    """
    Parsing arguments for cnr: r1
    """
    example = '\n'.join([
        'Examples:',
        '1. support fastq input',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -g dm6',
        '2. for specific index',
        '$ python cnr_r1.py --fq1 f1.fq.gz --fq2 f2.fq.gz -o results -x bowtie2_index/te',
    ])
    parser = argparse.ArgumentParser(
        prog='cnr_r1',
        description='cnr_r1: for single PE reads',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--fq1', required=True,
        help='Fastq file, read1 of PE')
    parser.add_argument('-2', '--fq2', required=True,
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
    CnrR1(**args).run()


if __name__ == '__main__':
    main()

#
