#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
## RNAseq-pipeline

+ 1. qc (fastqc report)

  - a. fastqc + rmarkdown report (html)

+ 2. trimming (optional)

  - a. NSR library (-6,6)  

  - b. standard illumina RNAseq (0)

+ 3. mapping 

  - 0. remove MT_trRNA (optional)

  - a. map to genome (STAR)  

  - b. map to transposon (STAR, consensus sequence, optional)  

  - c. map to extra sequence (STAR, firefly, ..., sequence + annotation)

+ 4. quantification (featureCounts) 

  - a. count.txt

+ 5. report

  - a. fastqc report

  - b. trimming report  

  - c. mapping report (.txt)  

  - d. count table (.stat)


## step 2

+ 6. de analysis

DESeq2, require >= 2 replicates for each sample

  - a. DESeq2 output

  - b. report (high-quality plots *.pdf, *.rda)

  - c. bigWig (compare between samples)

  
directory structure:

mapping
├── gene
│   ├── qc
│   ├── mapping
│   ├── bigWig
│   ├── count
│   └── report
├── transposon
│   ├── qc
│   ├── mapping
│   ├── bigWig
│   ├── count
│   └── report
└── extra_genes
    ├── mapping
    ├── bigWig
    ├── count
    └── report

de_analysis (A vs B)
├── gene
│   ├── bigWig
│   ├── de_analysis (DESeq2 output)
│   ├── GO_analysis (plots, tables, clusterProfiler, Reactome)
│   └── report (DE gene list)
├── transposon
│   ├── bigWig
│   ├── de_analysis
│   └── report
└── extra_genes
    ├── bigWig
    ├── de_analysis
    └── report

"""

import os
import sys
import argparse
import collections
import shlex
import subprocess
from arguments import args_init
from helper import *
from alignment import Alignment, Alignment_log, Alignment_stat
from hipipe_reporter import QC_reporter, Alignment_reporter
from utils_parser import DesignReader


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def get_args():
    parser = argparse.ArgumentParser(prog='RNAseq-pipeline', 
                                     description='DEseq analysis')
    parser.add_argument('--design', default=None,
        help='TAB-delimited file saving arguments, 6-column required: \
        <control-prefix> \
        <treatment-prefix> <fq_dir> <output_dir> <genome> <spikein>')
    parser.add_argument('-c1', nargs='+', metavar='control_fq', default=None,
        help='FASTQ files of Control sample replicates')
    parser.add_argument('-t1', nargs='+', metavar='treatment_fq', default=None,
        help='FASTQ files of Treatment sample, replicates')
    parser.add_argument('-c2', nargs='+', required=False, metavar='control_fq_2',
        default=None,
        help='the mate of paired-end reads for Control sample')
    parser.add_argument('-t2', nargs='+', required=False, metavar='treatment_fq_2',
        default=None,
        help='the mate of paired-end reads for Treatment sample')
    parser.add_argument('-o', '--path-out', default=None, dest='path_out',
        help='The directory to save results, default: [cwd]')
    parser.add_argument('-C', metavar='Control_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-T', metavar='Treatment_NAME', default=None,
        help='Name control samples')
    parser.add_argument('-g', '--genome', default='dm6',
        # choices=['dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10', 'GRCh38'],
        help='Reference genome : dm3, dm6, hg19, hg39, mm10, default: dm6')
    parser.add_argument('-k', '--spikein', default=None, 
        # choices=[None, 'dm3', 'dm6', 'hg19', 'hg38', 'mm9', 'mm10', 'GRCh38'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--gtf', default=None,
        help='genome annotation file in GTF format, from ensembl')
    parser.add_argument('--align-to-te', action='store_true',
        dest='align_to_te', help='if spcified, align reads to TE')
    parser.add_argument('--te-index', default=None, dest='te_index',
        help='path to the align index for TE, default: None')
    parser.add_argument('--te-gtf', default=None, dest='te_gtf',
        help='required if use --te-index, the path to GTF file')
    parser.add_argument('-x', '--extra-index', nargs='+', dest='extra_index',
        help='Provide extra alignment index(es) for alignment, support multiple\
        indexes. eg. Transposon, tRNA, rRNA and so on.')
    parser.add_argument('--extra-gtf', default=None, dest='extra_gtf',
        help='genome annotation file in GTF format, from ensembl')
    parser.add_argument('--include-multi-reads', action='store_true', 
        dest='include_multi_reads',
        help='if specified, including multi-mapped reads for further analysis')
    parser.add_argument('-p', '--threads', type=int, default=8,
        help='Number of threads to use, default: 8')
    parser.add_argument('-s', metavar='strandness', type=int,
        default=2, choices=[0, 1, 2],
        help='strandness, for NSR, dUTP: 1=anti, 2=sens, 0=ignore, default:2\
        This is for featureCounts, dUTP strand-specific mode, \
        read2 is sense strand.')
    parser.add_argument('--library-type', default='fr-firststrand',
        help='fr-unstranded: Standard Illumina Reads from the left-most end of \
        the fragment (in transcript coordinates) map to the transcript strand, \
        and the right-most end maps to the opposite strand. \
        fr-firststrand:dUTP, NSR, NNSR Same as above except we enforce the rule \
        that the right-most end of the fragment (in transcript coordinates) is \
        the first sequenced (or only sequenced for single-end reads). Equivalently, \
        it is assumed that only the strand generated during first strand synthesis \
        is sequenced. fr-secondstrand: Ligation, Standard SOLiD Same as above except \
        we enforce the rule that the left-most end of the fragment (in transcript \
        coordinates) is the first sequenced (or only sequenced for single-end reads). \
        Equivalently, it is assumed that only the strand generated during second \
        strand synthesis is sequenced.')
    parser.add_argument('--aligner', default='STAR', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: STAR')
    parser.add_argument('--bin-size', dest='bin_size', metavar='binsize', 
        type=int, default=50,
        help='binsize of bigWig, default: 50')
    parser.add_argument('--genome-path', dest='genome_path',
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--log-level', default='INFO', dest='log_level',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)

    return args


def init_rnaseq_project(x, analysis_type=1):
    """Create directories for RNAseq project
    :args analysis_type, 1 mapping
    :args analysis_tyep, 2 de_analysis
    x str, path to directory of project
    analysis_type int, 1=mapping, 2=de_analysis
    group 1, 1=gene, 2=te, 3=other
    """
    assert isinstance(x, str)
    path_dict = {
        1: ['mapping', 'bigWig', 'count', 'report'],
        2: ['bigWig', 'de_analysis', 'report']
    }

    group_dict = {
        1: 'gene',
        2: 'transposon',
        3: 'extra'
    }

    # choose A or B
    sub_path = path_dict.get(analysis_type, None)
    if sub_path is None:
        raise Exception('rnaseq-group, expect 1 or 2, not None')

    # create directories
    project_path = []
    out_dict = collections.defaultdict(dict)
    for group in [1, 2, 3]:
        # group
        group_name = group_dict.get(group, None)
        if group_name is None:
            raise Exception('[1,2,3] for group supported')
        # path
        for p in sub_path:
            p_path = os.path.join(x, group_name, p)
            is_path(p_path)
            project_path.append(p_path)
            # save to dict
            out_dict[group_name][p] = p_path

    # return values
    return out_dict


class fc(object):
    """Parsing the output of featureCounts
    :strandness
    :assigned%, warning if <60%
    """

    def __init__(self, x):
        """Input the .txt output of featureCounts
        Including *.summary in the same directory
        """
        self.txt = x
        self.summary = x + '.summary'
        if not os.path.exists(self.summary):
            log.warning('summary file not exists: {}'.format(self.summary))


    def assign(self):
        """Parsing the summary file for assign%
        warning if assign% < 60%
        """
        df = pd.read_csv(self.summary, '\t', index_col=0)
        df.columns = [os.path.basename(i) for i in df.columns.tolist()]
        ## total
        total = df.sum(axis=0, skipna=True)
        ## assign
        assign = df.loc['Assigned', ] / total * 100
        ## pct
        assign_df = assign.to_frame('assigned%')
        ## minimum
        min_pct = assign.min()
        if assign.min() < 60:
            log.warning('!!!! most of the reads are not assigned, Visual checking recommended !!!!')

        return assign_df

    def command(self):
        """Parsing the command line from .txt output
        :1-st line in file

        ## example
        # Program:featureCounts v1.6.1; Command:"featureCounts" "-M" "-O" 
        "--fraction" "-a" "dm3.ensembl.gtf" "-s" "2" "-p" "-C" "-B" 
        "-o" "count.txt" "bam1" "bam2"
        """
        with open(self.x, 'rt') as fi:
            line = next(fi) # the first line

        version, cmd_line = line.strip().split(';')
        version = version.split(' ')[2]
        cmd_line = re.sub('"', '', cmd_line.strip())

        return [version, cmd_line]


    def bam_files(self):
        """Parsing the bam files from summary
        the first line
        """
        with open(self.summary, 'rt') as fi:
            line = next(fi) # the first line

        bam_files = line.strip().split(' ')
        bam_files = [os.path.basename(i) for i in bam_files]
        # remove the 1-st item (Status)
        bam_files.pop(0)

        return bam_files


def run_featureCounts(gtf, bam_files, out, strandness=0, 
    threads=8, overwrite=False):
    """Count reads on each genes using featureCounts
    gtf, annotation file in GTF format
    bam, one or multiple alignemnt in BAM format
    out, count.txt
    featureCounts -M -O --fraction -T <cpu> -a <gtf> -s <0|1|2> <out.txt> bam.files

    ## paired-end options
    -p -C
    """
    ## init
    assert os.path.exists(gtf)
    assert isinstance(bam_files, list)
    assert isinstance(threads, int)
    log.info('running featureCounts')
    fc_exe = which('featureCounts')
    if fc_exe is None:
        raise Exception('featureCounts, not found in your $PATH')
    
    ## check BAM files (should be indexed)
    bam_check_flag = [BAM(b).index() for b in bam_files]
    if all(bam_check_flag) is False:
        raise Exception('failed to index bam files')

    ## save BAM prefix to file
    bam_prefix = [file_prefix(b)[0] for b in bam_files]
    bam_list = os.path.join(os.path.dirname(out), 'bam_ids.txt')
    with open(bam_list, 'wt') as fo:
        for b in bam_prefix:
            fo.write(b + '\n')

    ## check BAM is PE or SE
    # print(bam_files)
    pe_num = pysam.view('-c', '-f', '1', bam_files[0], catch_stdout=True)
    pe_num = int(pe_num)
    if pe_num > 0:
        pe_flag = True
    else:
        pe_flag = False

    ## prepare command
    bam_line = ' '.join(bam_files)
    fc_log = os.path.splitext(out)[0] + '.featureCounts.log'
    cmd = '{} -s {} -T {} '
    if pe_flag:
        cmd += '-p -C -B '
    cmd += '-M -O --fraction -g gene_id -t exon '
    cmd += '-a {} '
    cmd += '-o {} {} 2> {}'
    cmd = cmd.format(
        fc_exe,
        strandness,
        threads,
        gtf,        
        out,
        bam_line,
        fc_log)

    if os.path.exists(out) and overwrite is False:
        log.info('featureCounts skipped, file exists: {}'.format(out))
    else:
        run_shell_cmd(cmd)

        if not os.path.exists(out):
            raise Exception('featureCounts failed, file not found: {}'.format(out))

    # check output
    df_assign = fc(out).assign()
    log.info(df_assign)

    return out


def gene_aligner(fq1_files, smp_name, args, fq2_files=None):
    """Mapping reads to genome
    control or treatment
    args dict, the arguments of pipeline
    check index
    1. rRNA
    2. genome
    3. spike-in-rRNA
    4. spike-in
    """
    project_path = init_rnaseq_project(args['path_out'], analysis_type=1)
    gene_align_path = project_path['gene']

    ## qc-report
    qc_path = os.path.join(gene_align_path['report'], 'qc')
    # QC_reporter(fq1_files, qc_path).run()

    ## update args
    args['fq1'] = fq1_files
    args['fq2'] = fq2_files
    args['path_out'] = gene_align_path['mapping']
    args['smp_name'] = smp_name
    args['align_to_te'] = False

    ## run alignment
    map_bam_list = Alignment(**args).run()

    ## filt map_genome
    map_bam = []
    for i in map_bam_list:
        for k in i:
            if k.endswith('map_' + args['genome'] + '.bam'):
                map_bam.append(k)

    # # create bigWig files
    # for bam in map_bam:
    #     bam2bigwig(
    #         bam=bam, 
    #         genome=args['genome'], 
    #         path_out=gene_align_path['bigWig'],
    #         strandness=args['s'], 
    #         binsize=args['bin_size'],
    #         overwrite=args['overwrite'])                

    return map_bam


def te_aligner(fq1_files, smp_name, args, fq2_files=None):
    """Mapping reads to genome
    control or treatment
    args dict, the arguments of pipeline
    check index
    1. rRNA
    2. genome
    3. spike-in-rRNA
    4. spike-in
    """
    project_path = init_rnaseq_project(args['path_out'], analysis_type=1)
    te_align_path = project_path['transposon']

    args['extra_index'] = None # pre-build

    # ## qc-report
    # qc_path = os.path.join(te_align_path['report'], 'qc')
    # QC_reporter(fq1_files, qc_path).run() ## skip, run in gene_aligner

    ## update args
    args['fq1'] = fq1_files
    args['fq2'] = fq2_files
    args['path_out'] = te_align_path['mapping']
    args['smp_name'] = smp_name
    args['align_to_te'] = True

    # extra small genome
    small_genome = args['small_genome']
    args['small_genome'] = True

    ## run alignment
    map_bam_list = Alignment(**args).run()
    map_bam = [item for sublist in map_bam_list for item in sublist]

    # create bigWig files
    # for bam in map_bam:
    #     bam2bigwig(
    #         bam=bam, 
    #         genome=args['genome'], 
    #         path_out=te_align_path['bigWig'],
    #         strandness=args['s'], 
    #         binsize=args['bin_size'],
    #         overwrite=args['overwrite'])

    ## return
    args['small_genome'] = small_genome

    return map_bam


def extra_aligner(fq1_files, smp_name, args, fq2_files=None):
    """Mapping reads to genome
    control or treatment
    args dict, the arguments of pipeline
    check index
    1. rRNA
    2. genome
    3. spike-in-rRNA
    4. spike-in
    """
    project_path = init_rnaseq_project(args['path_out'], analysis_type=1)
    extra_align_path = project_path['extra']

    ## qc-report
    qc_path = os.path.join(extra_align_path['report'], 'qc')
    # QC_reporter(fq1_files, qc_path).run()

    ## update args
    args['fq1'] = fq1_files
    args['fq2'] = fq2_files
    args['path_out'] = extra_align_path['mapping']
    args['smp_name'] = smp_name
    args['align_to_te'] = False

    # extra small genome, for STAR
    small_genome = args['small_genome']
    args['small_genome'] = True

    ## run alignment
    map_bam = Alignment(**args).run()

    ## return
    args['small_genome'] = small_genome

    ## return
    return map_bam


def deseq2_exe(control, treatment, path_out, genome, nameA=None, 
    nameB=None, group='gene', pvalue=0.05, path_suffix=None):
    """Run DESeq2 witb count.txt matrix input
    control str, featureCounts output
    treatment str, featureCounts output
    path_out str, path to output
    group str, gene|transposon|extra_genes
    using custom R script
    """
    log.info('running DESeq2 analysis')
    project_path = init_rnaseq_project(path_out, analysis_type=2) # 1=mapping, 2=de
    de_path = project_path[group] # gene|transposon|extra_genes
    
    Rscript_exe = which('Rscript')
    if Rscript_exe is None:
        raise Exception('Rscript, command not found in your $PATH')

    deseq2_script = os.path.join(sys.path[0], 'run_deseq2.R')
    if not os.path.exists(deseq2_script):
        raise Exception('R script not found: %s' % deseq2_script)

    deseq2_path = de_path['de_analysis']
    if isinstance(path_suffix, str):
        deseq2_path += '_' + path_suffix
    is_path(deseq2_path) # create path
    deseq2_log = os.path.join(deseq2_path, 'run_DESeq2.log')

    cmd = ' '.join([Rscript_exe, deseq2_script, 
        control, treatment, genome, deseq2_path, nameA, nameB, str(pvalue)])
    cmd += ' > {}'.format(deseq2_log)
    print(cmd)
    # print(cmd)
    run_shell_cmd(cmd)


def gene_rnaseq(args):
    """RNAseq pipeline analysis for genes
    require, gtf, bam, ...
    """
    log.info('running for genes')

    group = 'gene' # !!!!

    ###########################
    ## sense strand analysis ##
    ###########################
    ## control, args['c1']
    ctl_args = args.copy()
    ctl_args['align_to_te'] = False # required !!!!
    ctl_args['extra_index'] = None # required !!!!
    ctl_args['path_out'] = os.path.join(args['path_out'], args['C'])
    ctl_bam = gene_aligner(args['c1'], args['C'], ctl_args, args['c2'])
    ## count reads
    ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.sens.txt')
    run_featureCounts(
        gtf=args['gtf'],
        bam_files=ctl_bam,
        out=ctl_count,
        strandness=args['s'],
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## treatment, args['t1']
    tre_args = args.copy()
    tre_args['align_to_te'] = False # required !!!!
    tre_args['extra_index'] = None # required !!!!
    tre_args['path_out'] = os.path.join(args['path_out'], args['T'])
    tre_bam = gene_aligner(args['t1'], args['T'], tre_args, args['t2'])
    ## count reads
    tre_count = os.path.join(args['path_out'], args['T'], 'gene', 'count', 'count.sens.txt')
    run_featureCounts(
        gtf=args['gtf'], 
        bam_files=tre_bam, 
        out=tre_count, 
        strandness=args['s'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## de analysis using DESeq2
    de_path = os.path.join(args['path_out'], args['C'] + '.vs.' + args['T'])
    deseq2_exe(
        control=ctl_count, 
        treatment=tre_count, 
        path_out=de_path, 
        genome=args['genome'], 
        nameA=args['C'], 
        nameB=args['T'], 
        group=group,
        path_suffix='sense')

    ###############################
    ## antisense strand analysis ##
    ###############################
    # determine the strandness
    if args['s'] == 2:
        args['anti_strand'] = 1
    elif args['s'] == 1:
        args['anti_strand'] = 2
    else:
        args['anti_strand'] = 0

    ## count reads
    ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=args['gtf'], 
        bam_files=ctl_bam, 
        out=ctl_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## count reads
    tre_count = os.path.join(args['path_out'], args['T'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=args['gtf'], 
        bam_files=tre_bam, 
        out=tre_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## de analysis using DESeq2
    de_path = os.path.join(args['path_out'], args['C'] + '.vs.' + args['T'])
    deseq2_exe(
        control=ctl_count, 
        treatment=tre_count, 
        path_out=de_path, 
        genome=args['genome'], 
        nameA=args['C'], 
        nameB=args['T'], 
        group=group,
        path_suffix='antisense')    


def te_rnaseq(args, gtf):
    """RNAseq pipeline analysis for Transposon, (only support dm3)
    require, gtf, bam, ...
    """
    log.info('running for transposon')
    group = 'transposon' # !!!! only for dm3

    if not os.path.exists(gtf):
        log.warning('transposon analysis skipped, gtf not exists: {}'.format(gtf))
        return None
    # sys.exit(gtf)


    ###########################
    ## sense strand analysis ##
    ###########################
    ## control, args['c1']
    ctl_args = args.copy()
    ctl_args['align_to_te'] = True # required !!!!
    ctl_args['extra_index'] = None # required !!!!
    ctl_args['path_out'] = os.path.join(args['path_out'], args['C'])    
    ctl_bam = te_aligner(args['c1'], args['C'], ctl_args, args['c2'])
    ## count reads
    ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.sens.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=ctl_bam, 
        out=ctl_count, 
        strandness=args['s'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## treatment, args['t1']
    tre_args = args.copy()
    ctl_args['align_to_te'] = True # required !!!!
    ctl_args['extra_index'] = None # required !!!!
    tre_args['path_out'] = os.path.join(args['path_out'], args['T'])
    tre_bam = te_aligner(args['t1'], args['T'], tre_args, args['t2'])
    ## count reads
    tre_count = os.path.join(args['path_out'], args['T'], group, 'count', 'count.sens.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=tre_bam, 
        out=tre_count, 
        strandness=args['s'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## de analysis using DESeq2
    de_path = os.path.join(args['path_out'], args['C'] + '.vs.' + args['T'])
    deseq2_exe(
        control=ctl_count, 
        treatment=tre_count, 
        path_out=de_path, 
        genome=args['genome'], 
        nameA=args['C'], 
        nameB=args['T'], 
        group=group,
        path_suffix='sense')


    ###############################
    ## antisense strand analysis ##
    ###############################
    # determine the strandness
    if args['s'] == 2:
        args['anti_strand'] = 1
    elif args['s'] == 1:
        args['anti_strand'] = 2
    else:
        args['anti_strand'] = 0

    ## count reads
    ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=ctl_bam, 
        out=ctl_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## count reads
    tre_count = os.path.join(args['path_out'], args['T'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=tre_bam, 
        out=tre_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## de analysis using DESeq2
    de_path = os.path.join(args['path_out'], args['C'] + '.vs.' + args['T'])
    deseq2_exe(
        control=ctl_count, 
        treatment=tre_count, 
        path_out=de_path, 
        genome=args['genome'], 
        nameA=args['C'], 
        nameB=args['T'], 
        group=group,
        path_suffix='antisense')


def extra_rnaseq(args, gtf):
    """RNAseq pipeline analysis for extra sequences
    require gtf, bam, ...
    """
    log.info('running for other reference')
    group = 'extra'
    ## update gtf
    # args['gtf'] = gtf

    ## check if extra_index exists
    if args['extra_index'] is None:
        log.warning('extra analysis skipped, index not exists: {}'.format(args['extra_index']))
    elif not os.path.exists(gtf):
        log.warning('extra analysis skipped, gtf file not exists: {}'.foramt(gtf))
    else:
        ## control, args['c1']
        ctl_args = args.copy()
        ctl_args['path_out'] = os.path.join(args['path_out'], args['C'])    
        ctl_bam = extra_aligner(args['c1'], args['C'], ctl_args)
        ## flaten
        ctl_bam = [item for sublist in ctl_bam for item in sublist]

        ## count reads
        ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.txt')
        run_featureCounts(
            gtf=gtf, 
            bam_files=ctl_bam, 
            out=ctl_count, 
            strandness=args['s'], 
            threads=args['threads'], 
            overwrite=args['overwrite'])

        ## treatment, args['t1']
        tre_args = args.copy()
        tre_args['path_out'] = os.path.join(args['path_out'], args['T'])
        tre_bam = extra_aligner(args['t1'], args['T'], tre_args)
        tre_bam = [item for sublist in tre_bam for item in sublist]

        ## count reads
        tre_count = os.path.join(args['path_out'], args['T'], group, 'count', 'count.txt')
        run_featureCounts(
            gtf=gtf, 
            bam_files=tre_bam, 
            out=tre_count, 
            strandness=args['s'], 
            threads=args['threads'], 
            overwrite=args['overwrite'])

        ## de analysis using DESeq2
        de_path = os.path.join(args['path_out'], args['C'] + '_vs_' + args['T'])
        deseq2_exe(
            control=ctl_count, 
            treatment=tre_count, 
            path_out=de_path, 
            genome=args['genome'], 
            nameA=args['C'], 
            nameB=args['T'], 
            group=group,
            path_suffix='sense')

    ###############################
    ## antisense strand analysis ##
    ###############################
    # determine the strandness
    if args['s'] == 2:
        args['anti_strand'] = 1
    elif args['s'] == 1:
        args['anti_strand'] = 2
    else:
        args['anti_strand'] = 0

    ## count reads
    ctl_count = os.path.join(args['path_out'], args['C'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=ctl_bam, 
        out=ctl_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## count reads
    tre_count = os.path.join(args['path_out'], args['T'], group, 'count', 'count.anti.txt')
    run_featureCounts(
        gtf=gtf, 
        bam_files=tre_bam, 
        out=tre_count, 
        strandness=args['anti_strand'], 
        threads=args['threads'], 
        overwrite=args['overwrite'])

    ## de analysis using DESeq2
    de_path = os.path.join(args['path_out'], args['C'] + '.vs.' + args['T'])
    deseq2_exe(
        control=ctl_count, 
        treatment=tre_count, 
        path_out=de_path, 
        genome=args['genome'], 
        nameA=args['C'], 
        nameB=args['T'], 
        group=group,
        path_suffix='antisense')


def map_stat(path):
    """Stat mapping reads for each replicate
    and merge files
    input: path to control/treatment
    directory structure:
    directory
      |- control
      |    |-merged
      |    |-rep1
      |    |-rep2
      |
      |-treatment
      |    |-merged
      |    |-rep1
      |    |-rep2
    """
    path_dirs = [os.path.join(path, f) for f in os.listdir(path)]
    ## exclude extra_mapping
    path_dirs = [f for f in path_dirs if os.path.isdir(f) and 
        not os.path.basename(f).startswith('extra_mapping')]
    frames = [Alignment_stat(f).stat for f in path_dirs]

    frames = [f for f in frames if isinstance(f, pd.DataFrame)]
    if len(frames) > 0:
        df = pd.concat(frames, axis=0)
        return df
    else:
        return None


def run(args):
    """Main for RNAseq analysis pipeline"""
    # args = args_init(vars(get_args()), align=True)

    log.info('running RNAseq pipeline')

    ## default arguments, for RNAseq2 only
    args['align_to_rRNA'] = True

    ## multireads should not be discarded
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/
    ## Genome Biology, 2016
    if args['include_multi_reads']:
        args['unique_only'] = False #  
    else:
        args['unique_only'] = True # default    

    # determine gtf file
    # gene gtf
    if args['gtf'] is None:
        args['gtf'] = Genome(**args).gene_gtf('refseq') # ucsc version

        # for GRCh38 genome, using GENCODE version 
        # if args['genome'] == 'GRCh38':
        if args['genome'] in ['GRCh38', 'GRCm38']:
            args['gtf'] = Genome(**args).gene_gtf('ensembl') # gencode

    # te gtf
    if args['te_gtf'] is None:
        args['te_gtf'] = Genome(**args).te_gtf()

    # print(args['gtf'])

    ## update prefix
    ctl_prefix = str_common([os.path.basename(f) for f in args['c1']])
    ctl_prefix = ctl_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args['C'] is None:
        args['C'] = ctl_prefix
    tre_prefix = str_common([os.path.basename(f) for f in args['t1']])
    tre_prefix = tre_prefix.rstrip('r|R|rep|Rep').rstrip('_|.')
    if args['T'] is None:
        args['T'] = tre_prefix
    
    ## run pipeline
    if args['extra_index']:
        extra_rnaseq(args, args['extra_gtf'])
    elif args['align_to_te']:
        te_rnaseq(args, args['te_gtf'])
    else:
        gene_rnaseq(args)

    log.info('finish')


def main():
    args = args_init(vars(get_args()), align=True)

    if args['design'] is None:
        # require -c1, -t1, -g
        if args['c1'] is None or args['t1'] is None:
            log.error('require --design, or -c1, -t1')
        else:
            if not supportedGenome(args['genome']):
                log.error('genome not supported: {}'.foramt(args['genome']))
            run(args)
    else:
        design = DesignReader(args['design'])
        # save to file
        design.to_json()
        # search dict
        d = design.to_dict()

        for k in d:
            log.info(k)
            dict_k = d[k]
            args['c1'], args['c2'] = dict_k['control']
            args['t1'], args['t2'] = dict_k['treatment']
            args['c2'] = args['c2'] if len(args['c2']) > 0 else None
            args['t2'] = args['t2'] if len(args['t2']) > 0 else None
            args['C'] = dict_k['control_name']
            args['T'] = dict_k['treatment_name']
            args['genome'] = dict_k['genome']
            args['path_out'] = dict_k['path_out']
            args['spikein'] = dict_k['spikein']
            
            ## supported
            if not supportedGenome(args['genome']):
                log.info('genome not supported: {}'.format(args['genome']))
                continue
            run(args)
        # print(d)


if __name__ == '__main__':
    main()


# EOF
