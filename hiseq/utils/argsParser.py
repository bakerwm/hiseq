
# -*- coding: utf-8 -*-


"""
Parse arguments from command line

- trimmer
- align

"""

import argparse
import pathlib




def add_sheet_args():
    """
    Prepare sample sheet for Demx
    output:
    1. sample_name,i7,i5,barcode
    2. i7_name,i7,reads
    3. sample_name,NULL,NULL,barcode (bc only)
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Prepare sample sheet for demx/demx2',
        epilog='''Description:
YY00.xlsx : required columns, ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']

Output:
1. sample_name,i7,i5,barcode
2. i7_name,i7,reads
3. sample_name,NULL,NULL,barcode

Example:
hiseq sheet -s YY00.xlsx -o data'''
    )
    parser.add_argument('-s', '--xlsx-table', dest='x', required=True,
        help='sample table in xlsx format, eg: YY00.xlsx')
    parser.add_argument('-o', '--outdir', dest='outdir',
        help='directory to save the reulsts')
    return parser


def add_demx_args():
    """
    Demultiplexing
    """
    parser = argparse.ArgumentParser(description='hiseq demx')
    parser.add_argument('-1', '--fq1', required=True,
        help='read1 in fastq format, gzipped')
    parser.add_argument('-2', '--fq2',
        help='read2 in fastq format, gzipped, (optional)')
    parser.add_argument('-o', '--outdir', required=True,
        help='directory to save the reulsts')
    parser.add_argument('-s', '--index-table', dest='index_table', required=True,
        help='index list in csv format, [filename,i7,i5,barcode]')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demostration, default: off')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-x', '--barcode-in-read2', action='store_true',
        help='barcode in read2')
    parser.add_argument('-l', '--barcode-n-left', type=int, dest='barcode_n_left',
        default=0, help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int, dest='barcode_n_right',
        default=0, help='bases locate on the right of barcode')
    parser.add_argument('-p', '--threads', type=int, default=1,
        help='number of threads, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', type=int, dest='parallel_jobs',
        default=1, help='number of josb run in parallel, default: [1]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def add_demx2_args():
    """
    Demultiplexing, multi barcode files
    """
    parser = argparse.ArgumentParser(description='hiseq demx2')
    parser.add_argument('-s', '--xlsx-table', dest='x', required=True,
        help="Sample table in (xlsx|csv) format; xlsx: require the columns\
        ['Sample_name*', 'P7_index_id*', 'Barcode_id*', 'Reads, M']; \
        csv: require the columns: ['name', 'i7', 'i5', 'bc', 'reads'] \
        the csv file could be `hiseq sheet -s a.xlsx -o data` output: *.demx.csv")
    parser.add_argument('-d', '--datadir', dest='datadir', required=True,
        help='Directory saving the fastq files')
    parser.add_argument('-o', '--outdir', dest='outdir',
        help='directory to save the reulsts')
    parser.add_argument('--demo', action='store_true',
        help='run demo (1M reads) for demostration, default: off')
    parser.add_argument('-m', '--mismatch', type=int, default=0,
        help='mismatches allowed to search index, default: [0]')
    parser.add_argument('-x', '--barcode-in-read2', dest='barcode_in_read2',
        action='store_true', help='barcode in read2')
    parser.add_argument('-l', '--barcode-n-left', type=int, dest='barcode_n_left',
        default=0, help='bases locate on the left of barcode')
    parser.add_argument('-r', '--barcode-n-right', type=int, dest='barcode_n_right',
        default=0, help='bases locate on the right of barcode')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite exists files, default: off')
    return parser


def add_qc_args():
    """
    utils:
      - fastqc
    """
    parser = argparse.ArgumentParser(
        description='hiseq qc, fastqc')
    parser.add_argument('-i', '--fq', nargs='+', required=True,
        help='reads in FASTQ files, or directory contains fastq files')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')
    parser.add_argument('--fastqc', default='fastqc',
        help='The path to the fastqc command, default: [fastqc]')

    parser.add_argument('-f', '--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, only for multiple fastq files, default: [1]')
    return parser


def add_p7_args():
    """
    utils:
      - fastqc
    """
    parser = argparse.ArgumentParser(
        description='hiseq p7')
    parser.add_argument('-i', '--fq', nargs='+', required=True,
        help='reads in FASTQ files, or directory contains fastq files')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')
    parser.add_argument('-f', '--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('-s', '--save-seq', dest='save_seq', 
        action='store_true',
        help='Save the i7 sequence to file')
    return parser


def add_trim_args():
    """
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser = argparse.ArgumentParser(
        description='hiseq qc, trim adapters and qc')
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='The read2 of pair-end reads')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results.')

    parser.add_argument('--library-type', dest='library_type', default=None,
        type=str, choices=['TruSeq', 'Nextera', 'smallRNA'],
        help='Type of the library structure, \
        TruSeq, TruSeq standard library \
        Nextera, Tn5 standard library, \
        smallRNA, small RNA library')
    parser.add_argument('-m', '--len_min', default=15, metavar='len_min',
        type=int, help='Minimum length of reads after trimming, defualt [15]')

    parser.add_argument('--cut-to-length', default=0, dest='cut_to_length',
        type=int,
        help='cut reads to from right, default: [0], full length')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')

    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1,
        type=int, help='Number of jobs to run in parallel, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    ## global arguments
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')

    ## specific
    parser.add_argument('--rm-untrim', action='store_true',
        dest='rm_untrim',
        help='discard reads without adapter')
    parser.add_argument('--save-untrim', action='store_true',
        dest='save_untrim',
        help='Save untrim reads to file')
    parser.add_argument('--save-too-short', action='store_true',
        dest='save_too_short',
        help='Save too short reads to file')
    parser.add_argument('--save-too-long', action='store_true',
        dest='save_too_long',
        help='Save too short reads to file')
    parser.add_argument('--cut-before-trim', default='0',
        dest='cut_before_trim',
        help='cut n-bases before trimming adapter; positive value, \
        cut from left; minus value, cut from right, eg: 3 or -4 or 3,-4, \
        default [0]')
    parser.add_argument('--cut-after-trim', default='0',
        dest='cut_after_trim',
        help='cut n-bases after trimming adapter; positive value, \
        cut from left; minus value, cut from right, eg: 3 or -4 or 3,-4, \
        default [0]')

    parser.add_argument('-a', '--adapter3', default=None,
        help='3-Adapter sequence, default [].')
    parser.add_argument('-g', '--adapter5', default='',
        help='5-Adapter, default: None')

    ## PE arguments
    parser.add_argument('-A', '--Adapter3', default=None,
        help='The 3 adapter of read2, default []')
    parser.add_argument('-G', '--Adapter5', default=None,
        help='The 5 adapter of read1, default: None')
    return parser


def add_align_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(
        description='Align short reads to reference sequence')
    parser.add_argument('-1', '--fq1', nargs='+', required=True,
        help='path to HiSeq reads in FASTQ format, support multipe \
        files, separated by white spaces.')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='path to HiSeq read2 of pair-end reads, optional, support \
        multiple files separated by white spaces.')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', required=True, default=None,
        choices=[None, 'dm6', 'dm3', 'hg38', 'hg19', 'mm10', 'mm9'],
        help='Reference genome : dm6, dm3, hg38, hg19, mm10, mm9, default: dm6')
    parser.add_argument('-k', '--spikein', default=None,
        choices=[None, 'dm6', 'dm3', 'hg38', 'hg19', 'mm10', 'mm9'],
        help='Spike-in genome : dm6, dm3, hg38, hg19, mm10, mm9, default: None')
    parser.add_argument('--aligner', default='bowtie',
        choices=['bowtie', 'bowtie2', 'STAR', 'hisat2', 'bwa', 'kalisto', 'salmon'],
        help='Choose which aligner to use. default: bowtie')

    ## extra: index
    parser.add_argument('--index-list', nargs='+', dest='index_list', default=None,
        help='ignore genome/spikein, add index directly, default: []')
    parser.add_argument('--index-name', nargs='+', dest='index_name', default=None,
        help='names for the input index list')
    parser.add_argument('-x', '--extra-index', nargs='+', dest="extra_index", default=None,
        help='Extra index for alignment, default: []')
    parser.add_argument('-n', '--smp-name', dest='smp_name', required=False,
        help='Name of the experiment')
    parser.add_argument('--index-list-equal', action='store_true',
        help='Align reads to each index list in parallel, if specified')

    ## extra: para
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--extra-para', dest='extra_para', default=None, type=str,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. default: [None]')
    parser.add_argument('--n-map', dest='n_map', type=int, default=0,
        help='Report up to N alignments per read. use -k for bowtie and \
        bowtie2 (default 1), --outFilterMultimapNmax for STAR \
        (default 20).')
    parser.add_argument('--repeat-masked-genome', dest='repeat_masked_genome',
        action='store_true',
        help='map to repeat masked reference genome, data from EnsEMBL')
    parser.add_argument('--path_data',
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')

    ## extra rRNA
    parser.add_argument('--align-to-chrM', dest='align_to_chrM',
        action='store_true',
        help='if specified, align to Mitochondrial DNA before genome, for supported genomes')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome, for supported genomes')
    parser.add_argument('--align-to-MT-trRNA', dest='align_to_MT_trRNA',
        action='store_true',
        help='if specified, align to Mito, tRNA and rRNA before genome, for supported genomes')
    parser.add_argument('--genomeLoad', dest='genomeLoad',
        default='LoadAndRemove',
        choices=['NoSharedMemory', 'LoadAndKeep', 'LoadAndRemove', 'LoadAndExit', 'Remove', 'NoSharedMemory'],
        help='--genomeLoad for STAR, default: [LoadAndRemove]'),

    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads for each job, default: [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, only for multiple fastq files, default: [1]')
    return parser


def add_quant_args():
    """
    quantify features using featureCounts
    support input file: BAM + GTF/BED/GFF ...
    """
    parser = argparse.ArgumentParser(
        description='quant reads')
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='BAM files, support multiple files separated by \
        white spaces')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=8, type=int,
        help='Number of threads to launch, default: 8.')
    return parser


def add_peak_args():
    """
    quantify features using featureCounts
    support input file: BAM + GTF/BED/GFF ...
    """
    parser = argparse.ArgumentParser(
        description='call peaks')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files, from IP sample')
    parser.add_argument('-c', '--control', nargs='+', required=False,
        help='BAM files for control sample, optional')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-n', '--name', default=None,
        help='The prefix of output files, default: None')
    parser.add_argument('-g', '--genome', required=False, default=None,
        choices=[None, 'dm6', 'dm3', 'hg38', 'hg19', 'mm10', 'mm9'],
        help='Reference genome : dm6, dm3, hg38, hg19, mm10, mm9, default: dm6')
    parser.add_argument('-gs', '--genome-size', required=False, type=int, default=0,
        help='The effective genome size, auto get from --genome, default')
    parser.add_argument('--is-atac', action='store_true',
        help='Change parameters for ATACseq alignment')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=8, type=int,
        help='Number of threads to launch, default: 8.')
    return parser


def add_motif_args():
    """
    quantify features using featureCounts
    support input file: BAM + GTF/BED/GFF ...
    """
    parser = argparse.ArgumentParser(
        description='call motif')
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='BAM files, support multiple files separated by \
        white spaces')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=8, type=int,
        help='Number of threads to launch, default: 8.')
    return parser


def add_rnaseq_args():
    """
    Arguments for RNAseq pipeline
    """
    parser = argparse.ArgumentParser(
        description='RNA-seq pipeline')
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='Create design for fastq files')
    parser.add_argument('-d', '--design', default=None,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('--wt', nargs='+', default=None, dest='wildtype',
        help='read1 fq files for control sample')
    parser.add_argument('--wt-fq2', nargs='+', dest='wildtype_fq2',
        default=None, help='read2 fq files for control sample')
    parser.add_argument('--mut', nargs='+', default=None, dest='mutant',
        help='read2 fq files for treatment sample')
    parser.add_argument('--mut-fq2', nargs='+', dest='mutant_fq2',
        default=None, help='read2 fq files for treatment sample')
    parser.add_argument('-1', '--fq1', nargs='+', default=None,
        help='read1 files, (or read1 of PE reads)')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='read2 of PE reads')
    parser.add_argument('-c', '--wt-dir', nargs='+', dest='wildtype_dir',
        default=None, help='path to the dirs of control samples')
    parser.add_argument('-t', '--mut-dir', nargs='+', dest='mutant_dir',
        default=None, help='path to the dirs of experiment samples')
    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg38, mm10, default: hg38')
    parser.add_argument('--gtf', default=None,
        help='The gtf file for quantification, defaut: genome.gtf (None)')
    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes')
    # optional arguments - 0
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--cut-to-length', dest='cut_to_length', default=0,
        type=int,
        help='cut the read to specific length, from right, default: [0], \
        not cut')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')
    # optional arguments - 1
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    ## extra: index
    parser.add_argument('-k', '--spikein', default=None,
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--align-to-rRNA', action='store_true',
        help='Align to rRNA')
    parser.add_argument('--aligner', default='STAR',
        choices=['STAR', 'bowtie', 'bowtie2', 'bwa', 'hisat2', 'kallisto',
            'salmon'],
        help='Aligner option: [STAR, bowtie, bowtie2, bwa], default: [STAR]')
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


def add_rnaseq_args2():
    """
    Arguments for RNAseq pipeline: packed
    """
    parser = argparse.ArgumentParser(
        description='RNA-seq pipeline')
    parser.add_argument('-d', '--design', required=True,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('-o', '--outdir', required=True,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', required=True, default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: hg19')
    parser.add_argument('--threads', default=1, type=int,
        help='Number of threads for each job, default [1]')
    parser.add_argument('-m', '--mode', default='gtp',
        choices=['g', 't', 'p', 'gt', 'gp', 'tp', 'gtp'],
        help='Run for g:gene, t:te, p:piRNA_cluster, default: [gtp]')
    return parser


def add_atac_args():
    """
    Arguments for ATAC-seq pipeline
    """
    parser = argparse.ArgumentParser(
        description='ATACseq pipeline')
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='Create design for fastq files')
    parser.add_argument('-d', '--design', default=None,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('--fq-dir', dest='fq_dir', default=None,
        help='directory of fastq files, for --build-design')
    parser.add_argument('-1', '--fq1', nargs='+', default=None,
        help='read1 files, (or read1 of PE reads)')
    parser.add_argument('-2', '--fq2', nargs='+', default=None,
        help='read2 of PE reads')

    parser.add_argument('-o', '--outdir', default=None,
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: hg19')
    parser.add_argument('--spikein', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: hg19')
    parser.add_argument('--spikein-index', dest='spikein_index', default=None,
        help='Index for Spikein')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')

    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes')

    # optional arguments - 1
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    parser.add_argument('--cut-to-length', dest='cut_to_length', default=0, type=int,
        help='cut the read to specific length, from right, default: [0], not cut')
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')

    return parser


def add_chipseq_args():
    """
    Arguments for ChIP-seq pipeline
    """
    parser = argparse.ArgumentParser(
        description='ChIPseq pipeline')
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='Create design for fastq files')
    parser.add_argument('-d', '--design', default=None,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('--ip', nargs='+', default=None,
        help='fastq for IP of ChIPseq, read1 of PE')
    parser.add_argument('--input', nargs='+', default=None,
        help='fastq for Input of ChIPseq, read1 of PE')
    parser.add_argument('--ip-fq2', nargs='+', dest='ip_fq2', default=None,
        help='fastq for IP of ChIPseq, read2 of PE, optional')
    parser.add_argument('--input-fq2', nargs='+', dest='input_fq2', default=None,
        help='fastq for Input of ChIPseq, read2 of PE, optional')

    parser.add_argument('-o', '--outdir', default=str(pathlib.Path.cwd()),
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: hg19')

    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')

    # optional arguments - 1
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    ## extra: index
    parser.add_argument('-x', '--extra-index', nargs='+', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')

    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. default: [None]')
    parser.add_argument('--copy-raw-fq', dest='copy_raw_data',
        action='store_true',
        help='whether copy the raw fastq files to output')
    return parser


def add_cnr_args():
    """
    Arguments for CnR pipeline
    """
    parser = argparse.ArgumentParser(
        description='CUT&RUN pipeline')
    parser.add_argument('-b', '--build-design', dest='build_design',
        action='store_true',
        help='Create design for fastq files')
    parser.add_argument('-d', '--design', default=None,
        help='design for RNAseq, json format, ignore fq1, fq2')
    parser.add_argument('--ip', nargs='+', default=None,
        help='fastq for IP of ChIPseq, read1 of PE')
    parser.add_argument('--input', nargs='+', default=None,
        help='fastq for Input of ChIPseq, read1 of PE')
    parser.add_argument('--ip-fq2', nargs='+', dest='ip_fq2', default=None,
        help='fastq for IP of ChIPseq, read2 of PE, optional')
    parser.add_argument('--input-fq2', nargs='+', dest='input_fq2',
        default=None,
        help='fastq for Input of ChIPseq, read2 of PE, optional')

    parser.add_argument('-o', '--outdir', default=str(pathlib.Path.cwd()),
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, \
        default: hg19')
    parser.add_argument('--spikein', default=None,
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, \
        default: hg19')
    parser.add_argument('--spikein-index', dest='spikein_index', default=None,
        help='Index for Spikein')

    parser.add_argument('--cut-to-length', dest='cut_to_length', default=0,
        type=int,
        help='cut the read to specific length, from right, default: [0], \
        not cut')
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--gene-bed', dest='gene_bed', default=None,
        help='The BED or GTF of genes')

    # optional arguments - 1
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('-j', '--parallel-jobs', default=1, type=int,
        dest='parallel_jobs',
        help='Number of jobs run in parallel, default: [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    ## extra: index
    parser.add_argument('--aligner', default='bowtie2',
        help='The alignment tool for CnR pipeline, default: [bowtie2]')
    parser.add_argument('-x', '--extra-index', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--recursive', action='store_true',
        help='trim adapter recursively')

    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. \
        default: [None]')

    return parser


##################################
## Utils
def add_trackhub_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog='get_trackhub',
        description='Generate trackhub for bigWig and bigBed files',
        epilog='Example: \n\
               python get_trackhub.py --config config.yaml')
    parser.add_argument('-d', '--demo', action='store_true',
        help='Show the tutorial')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
        help='Generate trackhub files, Do not copy track files to remote_dir')
    parser.add_argument('-c', '--config', default=None,
        help='Config file, run -d, generate the template')
    # parser.add_argument('-i', '--data-dir', required=False, dest='data_dir',
    #     help='The directory of bigWig and bigBed files')
    # parser.add_argument('-s', '--subgroup-config', dest='subgroups_yaml',
    #     help='Define the subgroups, subgroups.yaml')
    # parser.add_argument('-n', '--hub-name', metavar='hub_name', required=False,
    #     default=None, help='hub name')
    # parser.add_argument('-g', '--genome', metavar='GENOME', required=False,
    #     default='dm6', help='genome for the trackhub, UCSC genome build, \
    #     [hg19, hg38, mm9, mm10, dm3, dm6]')
    # parser.add_argument('-u', '--user', default='UCSC',
    #     help='Who maintain the trackhub')
    # parser.add_argument('-e', '--email', default='abc@abc.com',
    #     help='email of the maintainer')
    # parser.add_argument('--dry-run', dest='dry_run', action='store_true',
    #     help='Do not copy the files')
    # parser.add_argument('-l', '--short-label', default=None, dest='short_label',
    #     help='short label for the hub, default: [--hub-name]')
    # parser.add_argument('-L', '--long-label', default=None, dest='long_label',
    #     help='long label for the hub, default: [--hub]')
    # parser.add_argument('-d', '--description-url', default='',
    #     help='URL for description')
    # parser.add_argument('-r', '--recursive', action='store_true',
    #     help='search files in data_dir Recursively')
    return parser


def add_deseq_pair_args():
    """
    Run RNAseq compare
    """
    parser = argparse.ArgumentParser(description='hiseq deseq_pair -a dirA -b dirB -f gene -o outdir')
    parser.add_argument('-a', '--dirA', required=True,
        help='deseq_dir (a.vs.b) for groupA')
    parser.add_argument('-b', '--dirB', required=True,
        help='deseq_dir (a.vs.b) for groupB')
    parser.add_argument('-f', '--feature', default='gene',
        help='feqture of the RNAseq, gene|te|..., ')
    parser.add_argument('-o', '--outdir', default=None,
        help='directory to save the results')
    return parser


def add_go_args():
    """
    Run GO analysis
    """
    parser = argparse.ArgumentParser(description='hiseq go -i gene.xls -g dm6 -o output')
    parser.add_argument('-a', '--all',
        help='for all siginificantly changed genes, require "sig" in header, ignore: -i')
    parser.add_argument('-i', '--input',
        help='directory of deseq_dir (a.vs.b), or path to file, contain genes')
    parser.add_argument('-o', '--outdir',
        help='directory to save the GO results, only works if -i is gene_list')
    parser.add_argument('-g', '--genome',
        help='genome name, scientific_name prefer, eg: Drosophila melanogaster, also support: dm3/hg19/mm10')
    parser.add_argument('-f', '--foldChange',
        help='path to file, contains log2FoldChange, Gene')
    parser.add_argument('-t', '--feature', default='gene',
        help='feature of the analysis, only works if -i is deseq_dir')
    parser.add_argument('-c', '--ctl-vs-exp', default='1',
        choices=['1', '2'],
        help='1=exp/ctl, 2=ctl/exp; default:1')
    return parser


def add_fragsize_args():
    """
    required:
    bam
    outdir
    labels
    threads
    """
    parser = argparse.ArgumentParser(description='hiseq fragsize -i bam -o outdir')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-l', '--labels', nargs='+', default=None,
        help='label of the bam files')
    parser.add_argument('-p', '--threads', default=4, type=int,
        help='number of processes, default: [4]')
    return parser


def add_bam2bw_args():
    """
    required:
    bam
    outdir
    binsize
    """
    parser = argparse.ArgumentParser(description='hiseq bam2bw -i bam -g dm6 -o outdir')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-s', '--strandness', default=0, type=int,
        help='the strandness, 0=no, 1=fwd, 2=rev, 12=fwd&rev, default: [0]')
    parser.add_argument('-b', '--binsize', default=50, type=int,
        help='set binSize for bigWig, default: [50]')
    parser.add_argument('-g', '--genome', default=None,
        help='choose genome for the bam file, default: [None]')
    parser.add_argument('--scaleFactor', nargs='+', type=float, default=[1.0],
        help='The scaling factor, default: [1.0]')
    parser.add_argument('--normalizeUsing', default='None',
        choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
        help='Possible choices: RPKM, CPM, BPM, RPGC, None, default: [None] \
        see https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html, \
        for details from deeptools documentation.')
    parser.add_argument('-gs', '--genome-size', dest='genome_size', default=None,
        help='set the genome size for input genome, default: [None]' )
    parser.add_argument('-r', '--reference', default=None,
        help='the reference genome in fasta format, default: [None]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads, default: [4]')
    return parser


def add_bam2cor_args():
    """
    required:
    bam
    outdir
    """
    parser = argparse.ArgumentParser(description='hiseq bam2cor -i bam -o outdir')
    parser.add_argument('-i', '--bam', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-m', '--cor-method', default='pearson',
        choices=['pearson', 'spearman'],
        help='method to calculate correlation, default: [pearson]')
    parser.add_argument('-np', '--no-plot', dest='no_plot', action='store_false',
        help='do not make plots')
    parser.add_argument('-n', '--prefix', default=None,
        help='set the prefix for output files, default: [multibam]')
    parser.add_argument('-p', '--threads', default=1, type=int,
        help='number of threads, default: [1]')
    parser.add_argument('-b', '--binsize', default=50, type=int,
        help='set binSize for bigWig, default: [50]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    return parser


def add_peak2idr_args():
    """
    required:
    bam
    outdir
    """
    parser = argparse.ArgumentParser(description='hiseq peak2idr -i peak -o outdir')
    parser.add_argument('-i', '--peak', nargs='+', required=True,
        help='peak files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-m', '--cor-method', default='pearson',
        choices=['pearson', 'spearman'],
        help='method to calculate correlation, default: [pearson]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    return parser


def add_bed2overlap_args():
    """
    required:
            'peak': None,
            'outdir': str(pathlib.Path.cwd()),
            'flag': False,
            'prefix': None,
            'overwrite': False
    """
    parser = argparse.ArgumentParser(description='hiseq peak2idr -i peak -o outdir')
    parser.add_argument('-i', '--peak', nargs='+', required=True,
        help='peak files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-n', '--prefix', default=None,
        help='set the prefix for output files, default: [multibam]')
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Whether overwrite exists files')
    return parser


def add_sample_args():
    """
    required:
        'input':
        'outdir':
        'sample_size':
    """
    parser = argparse.ArgumentParser(description='hiseq sample')
    parser.add_argument('-i', '--fx', nargs='+', required=True,
        help='fastx files')
    parser.add_argument('-o', '--outdir', default=None,
        help='output directory to save results')
    parser.add_argument('-n', '--number', type=int, default=1000,
        help='Number of records, default: 1000')
    parser.add_argument('-r', '--random', action='store_true',
        help='Get random subset records'),
    parser.add_argument('-w', '--overwrite', action='store_true',
        help='Overwrite the exists files')
    parser.add_argument('-j', '--parallel-jobs', dest='parallel_jobs',
        default=1, type=int,
        help='Number of threads run in parallel, default [1]')
    return parser
