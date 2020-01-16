
# -*- coding: utf-8 -*-


"""
Parse arguments from command line

- trimmer
- align

"""

import argparse
       

def add_qc_args():
    """
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser = argparse.ArgumentParser(
        description='hiseq qc, trim adapters and qc')    
    parser.add_argument('-i', '--fq1', nargs='+', required=True, 
        help='reads in FASTQ files, support (*.gz), 1-4 files.')
    parser.add_argument('-o', '--outdir', default=None, 
        help='The directory to save results.')

    parser.add_argument('--library-type', dest='library_type', default='unknown',
        type=str, choices=['TruSeq', 'Nextera', 'smRNA', 'NSR', 'ChIPseq', 'iCLIP', 'eCLIP', 
            'CLIP_NSR', 'smRNA_NSR', 'unknown'],
        help='Type of the library structure, \
        TruSeq, TruSeq standard libraries \
        Nextera, Tn5 standard libraries, \
        NSR, (TruSeq), cut 7-nt at the left-end of both reads, \
        eCLIP, (TruSeq), cut 10-nt at left-end, 7-nt at-right end of read1, \
                 cut 7-nt at left-end and 10-nt at right-end of read2, \
                 This is Yulab version eCLIP, random barcode (N10) at P5, \
                 and barcode (6-nt + 1A) at P7 end. \
        iCLIP, (TruSeq), cut 9-nt at left-end of read1 (barcode) \
        determine the way to trim the raw reads, default: [unknown] \
        ignore --cut-after-trim \
        the following arguments are masked, if --library-type is not unknown: \
        --adapter3, --adapter5, --AD3, --AD5, --len-min, --rmdup, --cut-after-trim, ...')

    parser.add_argument('-m', '--len_min', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-a', '--adapter3', metavar='adapter', type=str,
    	default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',         
        help='3-Adapter, default: [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC].')    
    parser.add_argument('-g', '--adapter5', default='',
        help='5-Adapter, default: None')
    parser.add_argument('--read12', type=int, default=1,
        help='which one of PE reads, 1=read1, 2=read2, default: 1')

    parser.add_argument('--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    ## global arguments    
    parser.add_argument('-q', '--qual-min', default=20, type=int,
        dest='qual_min',
        help='The cutoff of base quality, default [20]')    
    parser.add_argument('-e', '--error-rate', default=0.1, type=float,
        dest='error_rate',
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', '--overlap', default=3, type=int,
        help='Required N bases overlap between reads and adapter, default [3]')
    parser.add_argument('-p', '--percent', default=80, type=int,
        help='minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('--rm-untrim', action='store_true', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--keep-name', action='store_true', dest='keep_name',
        help='if specified, do not change file names')
    
    ## extra arguments
    parser.add_argument('--adapter-sliding', dest='adapter_sliding', 
        action='store_true',
        help='Trim reads by sliding windows on adapter')
    parser.add_argument('--trim-times', dest='trim_times', type=int,
        default=1, help='Trim adapter from reads by N times, default:1')
    parser.add_argument('--double-trim', action='store_true', 
        dest='double_trim', help='if specified, trim adapters twice')
    parser.add_argument('--rmdup', action='store_true', dest='rmdup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--trim-to-length', default=0, metavar='max-length',
        dest='trim_to_length', type=int,
        help='trim reads from right, save the specific length of reads. \
              default: [0], 0=the full length')

    ## PE arguments
    parser.add_argument('--fq2', nargs='+', default=None, 
        help='The read2 of pair-end reads')
    parser.add_argument('-A', '--AD3', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        help='The 3 adapter of read2, default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-G', '--AD5', default=None,
        help='The 5 adapter of read1, default: None')
    return parser


def add_align_args():
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser = argparse.ArgumentParser(
        description='Align short reads to reference sequence')
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='path to HiSeq reads in FASTQ format, support multipe \
        files, separated by white spaces.')
    parser.add_argument('--fq2', nargs='+', default=None, 
        help='path to HiSeq read2 of pair-end reads, optional, support \
        multiple files separated by white spaces.')
    parser.add_argument('-o', '--outdir', default=None, 
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', required=True, default='dm6', 
        choices=[None, 'dm6', 'dm3', 'hg38', 'hg19', 'mm10', 'mm9'],
        help='Reference genome : dm6, dm3, hg38, hg19, mm10, mm9, default: dm6')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')

    ## extra: index
    parser.add_argument('-k', '--spikein', default=None, 
        choices=[None, 'dm6', 'dm3', 'hg38', 'hg19', 'mm10', 'mm9'],
        help='Spike-in genome : dm6, dm3, hg38, hg19, mm10, mm9, default: None')
    parser.add_argument('-x', '--extra-index', nargs='+', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('-n', '--smp_name', required=False,
        help='Name of the experiment')

    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
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
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')

    ## extra rRNA
    parser.add_argument('--align-to-chrM', dest='align_to_chrM',
        action='store_true',
        help='if specified, align to Mitochondrial DNA before genome')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome')
    parser.add_argument('--align-to-MT-trRNA', dest='align_to_MT_trRNA',
        action='store_true',
        help='if specified, align to Mito, tRNA and rRNA before genome')

    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=8, type=int, 
        help='Number of threads to launch, default: 8.')
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
    parser.add_argument('-i', '--fq1', nargs='+', required=True,
        help='BAM files, support multiple files separated by \
        white spaces')   
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


def add_atac_args():
    """
    Arguments for ATAC-seq pipeline
    """
    parser = argparse.ArgumentParser(
        description='ATAC-seq pipeline')
    parser.add_argument('-c', '--config', nargs='+', default=None,
        help='config, a tab-separated file; config for multiple ATACseq samples \
        in this format: [fq1  fq2  genome  outdir] \
        support multiple config.txt files \
        if specified, ignore [--design], [--fq1, --fq2, -o, --genome]')
    parser.add_argument('-d', '--design', nargs='+', default=None,
        help='design, a json file, config for one ATACseq sample; \
        supprot multiple design.json files \
        (ignore --fq1, --fq2, --genome, --outdir)')
    parser.add_argument('-1', '--fq1', nargs='+', default=None,
        help='read1 files, (or read1 of PE reads)')
    parser.add_argument('-2', '--fq2', nargs='+', default=None, 
        help='read2 of PE reads')
    parser.add_argument('-o', '--outdir', default=None, 
        help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-g', '--genome', default='dm6', 
        choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm9, mm10, default: hg19')

    # optional arguments - 0    
    parser.add_argument('-n', '--smp_name', required=False,
        help='Name of the experiment, works for only one input fastq file')
    parser.add_argument('--trimmed', action='store_true',
        help='specify if input files are trimmed')
    parser.add_argument('--copy-raw-fq', dest='copy_raw_data',
        action='store_true',
        help='whether copy the raw fastq files to output')
    # parser.add_argument('--skip-rmdup', dest='rmdup', # odd!!!!
    #     action='store_false',
    #     help='do not remove PCR duplicates')
    # conflict with trim()

    # optional arguments - 1    
    parser.add_argument('--threads', default=1, type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')

    # optional arguments - 2
    parser.add_argument('-m', '--len_min', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-a', '--adapter3', metavar='adapter', type=str,
        default='CTGTCTCTTATACACATCT',         
        help='3-Adapter, default: [CTGTCTCTTATACACATCT].')

    ## extra arguments - 3
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')

    ## extra: index
    parser.add_argument('-k', '--spikein', default=None, 
        choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', '--extra-index', nargs='+', dest="extra_index",
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')

    ## extra: para
    parser.add_argument('--extra-para', dest='extra_para', default=None,
        help='Extra parameters for aligner, eg: -X 2000 for bowtie2. default: [None]')

    ## extra: call-peak
    parser.add_argument('--genome-size', dest='genome_size', default=0,
        type=int, help='The genome size for the genome. default: [0]; use --genome')

    return parser





