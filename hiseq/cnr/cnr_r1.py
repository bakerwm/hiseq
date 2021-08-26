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
"""


import os
import sys
from hiseq.cnr.cnr_rp import CnrRp
from hiseq.cnr.utils import (
    hiseq_bam2bw, hiseq_pcr_dup,
    cnr_trim, cnr_align_genome,  cnr_align_spikein, cnr_call_peak, 
    qc_trim_summary, qc_align_summary, qc_lendist, 
    qc_frip, qc_tss_enrich, qc_genebody_enrich
)
from hiseq.utils.file import (
    check_fx_args, check_path, check_fx_paired, symlink_file, file_abspath,
    file_prefix, fx_name, Genome
)
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, init_cpu, read_hiseq
)
from hiseq.align.align_index import AlignIndex, check_index_args


class CnrR1(object):
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = CnrR1Config(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Config().dump(self.__dict__, self.config_yaml)


    def run_pipe(self):
        # symlink raw data
        raw_fq1, raw_fq2 = self.raw_fq_list
        symlink_file(self.fq1, raw_fq1, absolute_path=True)
        symlink_file(self.fq2, raw_fq2, absolute_path=True)
        cnr_trim(self.project_dir, '_r1')
        if isinstance(self.spikein_index, str):
            cnr_align_spikein(self.project_dir, '_r1')
        cnr_align_genome(self.project_dir, '_r1')
        hiseq_pcr_dup(self.project_dir, '_r1')
        cnr_call_peak(self.project_dir, '_r1')
        hiseq_bam2bw(self.project_dir, '_r1')
        if isinstance(self.spikein_index, str):
            cnr_align_spikein(self.project_dir, '_r1')
        # qc
        qc_trim_summary(self.project_dir, '_r1')
        qc_align_summary(self.project_dir, '_r1')
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
            'genome_path': None,
            'keep_tmp': None,
            'trimmed': False,
            'cut': False,
            'cut_to_length': 0,
            'recursive': False
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'cnr_r1'
        # output
        if self.outdir is None:
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)
        self.init_cut()
        self.init_fq()
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


    def init_fq(self):
        """
        required:
        1. fq1:str, fq2:str/None
        2. paired
        """
        if not check_fx_args(self.fq1, self.fq2):
            raise ValueError('fq1, fq2 not valide, check above message')
        self.fq1 = file_abspath(self.fq1)
        self.fq2 = file_abspath(self.fq2)
        self.is_paired = check_fx_paired(self.fq1, self.fq2)
        # auto: sample names
        if not isinstance(self.smp_name, str):
            self.smp_name = fx_name(self.fq1, fix_pe=self.is_paired, fix_unmap=True)


    # update: genome_size_file    
    def init_index(self):
        index_list = check_index_args(**self.__dict__)
        if len(index_list) == 0:
            raise ValueError('no index found')
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
        trim_prefix = os.path.join(self.clean_dir, self.smp_name)
        spikein_prefix = os.path.join(self.spikein_dir, self.smp_name)
        align_prefix = os.path.join(self.align_dir, self.smp_name)
        peak_prefix = os.path.join(self.peak_dir, self.smp_name)
        default_files = {
            'config_yaml': self.config_dir + '/config.yaml', # updated
            'report_log': self.report_dir + '/report.log',
            'report_html': self.report_dir + '/HiSeq_report.html',
            'raw_fq1': self.raw_dir + '/' + os.path.basename(self.fq1),
            'raw_fq2': self.raw_dir + '/' + os.path.basename(self.fq2),
            'clean_fq1': self.clean_dir + '/' + os.path.basename(self.fq1),
            'clean_fq2': self.clean_dir + '/' + os.path.basename(self.fq2),

            # basic files
            'config_yaml': os.path.join(self.config_dir, 'config.yaml'),
            'report_html': os.path.join(self.report_dir, 'HiSeq_report.html'),
            'bam_raw': self.bam_dir + '/' + self.project_name + '.raw.bam',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'peak': peak_prefix+'_peaks.narrowPeak',
            'peak_seacr': peak_prefix+'.stringent.bed',
            'peak_seacr_top001': peak_prefix+'.top0.01.stringent.bed',

            # trimming
            'trim_stat': trim_prefix+'.trim.stat',
            'trim_json': trim_prefix+'.trim.json',
            
            # align files (genome)
            'align_scale_json': self.bam_dir + '/' + 'scale.json',
            'pcr_dup_json': self.bam_dir + '/' + 'pcr_dup.json',
            'align_stat': align_prefix+'.align.stat',
            'align_json': align_prefix+'.align.json',
            'align_flagstat': align_prefix+'.align.flagstat',
            'unmap1': align_prefix+'.unmap.1.fastq',
            'unmap2': align_prefix+'.unmap.2.fastq',
            
            # spikein files
            'spikein_bam': spikein_prefix+'.bam',
            'spikein_scale_json': spikein_prefix+'.scale.json',
            'spikein_stat': spikein_prefix+'.align.stat',
            'spikein_json': spikein_prefix+'.align.json',
            'spikein_flagstat': spikein_prefix+'.align.flagstat',
            'spikein_unmap1': spikein_prefix+'.unmap.1.fastq',
            'spikein_unmap2': spikein_prefix+'.unmap.2.fastq',
            
            # qc
            'trim_summary_json': os.path.join(self.qc_dir, '00.trim_summary.json'),
            'align_summary_json': os.path.join(self.qc_dir, '01.alignment_summary.json'),
            'dup_summary_json': os.path.join(self.qc_dir, '01.pcr_dup_summary.json'),
            'lendist_csv': self.qc_dir + '/02.length_distribution.fragsize.csv',
            'lendist_txt': self.qc_dir + '/02.length_distribution.txt',
            'lendist_pdf': self.qc_dir + '/02.length_distribution.fragsize.pdf',
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
    CnrR1(**args).run()


if __name__ == '__main__':
    main()

#
