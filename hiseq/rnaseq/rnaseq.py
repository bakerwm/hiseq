#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

1. RNAseq for 1 index, and 1 transcript (gene, TE, piRNA cluster, reporter, ...)
# attention, norm_scale!?

2. Saving config in TOML format

Directory structure:

wildtype
  |- gene
  |- te
  |- piRNA
  |- reporter
  |- ...

mutant
  |- gene
  |- te
  |- piRNA
  |- reporter
  |- ...

mutant.vs.wildtype
  |- gene
  |- te
  |- piRNA
  |- reporter
  |- ...

subdirectory structure

wildtype/gene
  |- raw_data
  |- clean_data
  |- align
  |- count
  |- qc
  |- report
  |- spikein
  |    |- align
  |    |- count/norm.txt


mutant.vs.wildtype/gene
  |- count sense/antisense smp_name.txt   [sense, antisense]?
  |- deseq
  |- enrich
  |- report


## Create plots
1. QC
  - adapter content
  - mapping content
  - unique/multiple content
  - gene index

2. Reads on gene body?

3. DESeq
  - scatter plot
  - MA plot
  - Volcano
  - heatmap



2020-04-18
## 1. save arguments in self.
## 2. save default paths, files, in self.files
## 3. update all args
"""

import os
import re
import glob
import shutil
import tempfile
import collections
import pandas as pd
from multiprocessing import Pool
from Levenshtein import distance
import hiseq # for command file
from hiseq.utils.helper import *
from hiseq.trim.trimmer import Trim
from hiseq.rnaseq.alignment import Align, AlignIndex
from hiseq.atac.atac_utils import *
from hiseq.bam2bw.bam2bw import Bam2bw

def print_dict(d):
    d = collections.OrderedDict(sorted(d.items()))
    for k, v in d.items():
        print('{:>20s}: {}'.format(k, v))


def init_cpu(threads=1, parallel_jobs=1):
    """
    threads, CPUs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()

    max_jobs = int(n_cpu / 4.0)
    ## check parallel_jobs (max: 1/4 of n_cpus)
    if parallel_jobs > max_jobs: 
        log.warning('Too large, change parallel_jobs from {} to {}'.format(
            parallel_jobs, max_jobs))
        parallel_jobs = max_jobs

    ## check threads
    max_threads = int(0.8 * n_cpu / parallel_jobs)
    if threads * parallel_jobs > 0.8 * n_cpu:
        log.warning('Too large, change threads from {} to {}'.format(
            threads, max_threads))
        threads = max_threads

    return (threads, parallel_jobs)
    

def fq_paired(fq1, fq2):
    """
    Check fq1, fq2, (list)
    """
    if isinstance(fq1, str) and isinstance(fq2, str):
        chk = distance(fq1, fq2) == 1
    else:
        chk = False

    if not chk:
        log.warning('fastq pairing failed')

    return chk


class RNAseq(object):
    """
    The main port, for RNAseq analysis: N groups of RNAseq analysis
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RNAseqConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        # Toml(self.__dict__).to_toml(self.config_toml)


    def run_RNAseqRx_single(self, i):
        """
        parsing arguments as dict
        run RNAseqRx
        """
        args_local = self.__dict__
        args_local['parallel_jobs'] = 1 # force
        d = list(self.fq_groups.values())[i]
        args_local.update(d)
        RNAseqRx(**args_local).run()


    def run_RNAseqRx(self):
        """
        Run RNAseq, multiple groups (>=1)
        multiple projects
        """
        if len(self.fq_groups) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_RNAseqRx_single, range(len(self.fq_groups)))
        else:
            for i in range(len(self.fq_groups)):
                self.run_RNAseqRx_single(i)


    def run(self):
        """
        Run all
        """
        if self.hiseq_type == 'build_design':
            RNAseqDesign(**self.__dict__).run()
        elif self.hiseq_type == 'rnaseq_rx':
            self.run_RNAseqRx()
        elif self.hiseq_type == 'rnaseq_rn':
            RNAseqRn(**self.__dict__).run()
        elif self.hiseq_type == 'rnaseq_r1':
            RNAseqR1(**self.__dict__).run()
        else:
            raise ValueError('unknown hiseq_type: {}'.format(self.hiseq_type))


class RNAseqRx(object):
    """
    DEseq analysis, single pair: mutant vs wildtype
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RNAseqRxConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Toml(self.__dict__).to_toml(self.config_toml)


    def run_RNAseqRn(self):
        """
        Run RNAseqRn for wt/mut fastq files
        or copy file from RNAseqRn
        """
        c_mutant = RNAseqReader(self.mutant_dir)
        c_wildtype = RNAseqReader(self.wildtype_dir)

        if c_mutant.is_hiseq_rn and c_wildtype.is_hiseq_rn:
            log.info('run RNAseqRx for mutant_dir and wildtype_dir')
        elif isinstance(self.mutant, list) and isinstance(self.wildtype, list):
            args_local = self.__dict__.copy()
            rm_args = [
                'build_design',
                'design',
                'rep_list',
                'mutant_dir',
                'wildtype_dir',
                'mutant',
                'mutant_fq2',
                'wildtype',
                'wildtype_fq2',
                'smp_name'
            ]
            [args_local.pop(i, None) for i in rm_args] # remove

            # for mutant
            log.info('run RNAseqRn for wildtype')
            # for wildtype fastq files
            wildtype_args = args_local.copy()
            wildtype_args.update({
                'fq1': self.wildtype,
                'fq2': self.wildtype_fq2
            })
            RNAseqRn(**wildtype_args).run()

            # for mutant fastq files
            log.info('run RNAseqRn for mutant')
            mutant_args = args_local.copy()
            mutant_args.update({
                'fq1': self.mutant,
                'fq2': self.mutant_fq2
            })
            RNAseqRn(**mutant_args).run()
        else:
            log.error('RNAseqRx() failed, unknown args')

        # update
        self.mutant_bam_from = c_mutant.args.get('bam', None)
        self.mutant_bw_from = c_mutant.args.get('bw', None)
        self.wildtype_bam_from = c_wildtype.args.get('bam', None)
        self.wildtype_bw_from = c_wildtype.args.get('bw', None)


    def copy_bam_files(self):
        """
        Copy bam files
        """
        # file_symlink(self.mutant_bam_from, self.mutant_bam)
        # file_symlink(self.wildtype_bam_from, self.wildtype_bam)
        bam_list = RNAseqReader(self.project_dir).get_r1_file('bam')
        for bam in bam_list:
            file_symlink(bam, self.bam_dir)


    def copy_bw_files(self):
        """
        Copy bw files
        """
        # c_mutant = RNAseqReader(self.mutant_dir)
        # file_symlink(c_mutant.args.get('bw', None), self.mutant_bw)

        # c_wildtype = RNAseqReader(self.wildtype_dir)
        # file_symlink(c_wildtype.args.get('bw', None), self.wildtype_bw)

        bw_list = RNAseqReader(self.project_dir).get_r1_file('bw')
        for bw in bw_list:
            file_symlink(bw, self.bw_dir)


    def copy_count_files(self):
        # # sense strand
        # c_mutant = RNAseqReader(self.mutant_dir)
        # mut_sens = c_mutant.args.get('count_sens', None)
        # mut_anti = c_mutant.args.get('count_anti', None)
        # file_symlink(mut_sens, self.mutant_count_sens)
        # file_symlink(mut_anti, self.mutant_count_anti)

        # # anti strand
        # c_wildtype = RNAseqReader(self.wildtype_dir)
        # wt_sens = c_wildtype.args.get('count_sens', None)
        # wt_anti = c_wildtype.args.get('count_anti', None)
        # file_symlink(wt_sens, self.wildtype_count_sens)
        # file_symlink(wt_anti, self.wildtype_count_anti)

        c_sens = RNAseqReader(self.project_dir).get_r1_file('count_sens')
        c_anti = RNAseqReader(self.project_dir).get_r1_file('count_anti')

        try:
            for c in c_sens + c_anti:
                file_symlink(c, self.count_dir)
        except:
            log.error('copy count.txt files failed')


    def run_deseq(self):
        """DEseq analysis
        Using DESeq2, edgeR, ...
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        deseqR = os.path.join(pkg_dir, 'bin', 'run_rnaseq.R')
        deseq_stdout = os.path.join(self.deseq_dir, 'deseq.stdout')
        deseq_stderr = os.path.join(self.deseq_dir, 'deseq.stderr')
        cmd_shell = os.path.join(self.deseq_dir, 'cmd.sh')
        run_go = 1 if self.genome else 0
        cmd = 'Rscript {} {} {} 1> {} 2> {}'.format(
            deseqR,
            self.project_dir,
            run_go,
            deseq_stdout,
            deseq_stderr)
        with open(cmd_shell, 'wt') as w:
            w.write(cmd + '\n')
        if file_exists(self.deseq_fix_xls) and not self.overwrite:
            log.info('run_deseq() skipped, file exists: {}'.format(
                self.deseq_fix_xls))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('DESeq2 failed.')


    def report(self):
        """
        Create report for one sample
        html
        """
        RNAseqRp(self.project_dir).run()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        self.run_RNAseqRn()
        self.copy_bam_files()
        self.copy_bw_files()
        self.copy_count_files()
        self.run_deseq()
        self.report()


class RNAseqRn(object):
    """
    RNAseq analysis for multiple replicates, merge replicates
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = RNAseqRnConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Toml(self.__dict__).to_toml(self.config_toml)


    def get_bam_list(self):
        """
        get the proper_mapped.bam files
        """
        return [RNAseqReader(i).args.get('bam', None) for i in self.rep_list]


    def get_bw_list(self):
        """
        get the proper_mapped.bam files
        """
        return [RNAseqReader(i).args.get('bw', None) for i in self.rep_list]


    def merge_bam(self):
        """
        Merge replicates, BAM
        """
        self.bam_list = self.get_bam_list()

        cmd = ' '.join([
            'samtools merge -',
            ' '.join(self.bam_list),
            '| samtools sort -o {} -'.format(self.bam),
            '&& samtools index {}'.format(self.bam)])

        if os.path.exists(self.bam):
            log.info('merge_bam() skipped, file exists: {}'.format(
                self.bam))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.warning('merge_bam() failed.')

        ## calculate norm scale
        if file_exists(self.align_scale_txt):
            with open(self.align_scale_txt) as r:
                self.align_scale = eval(r.readline().strip())
        else:
            self.align_scale = self.cal_norm_scale(self.bam)
            with open(self.align_scale_txt, 'wt') as w:
                w.write('{:.4f}\n'.format(self.align_scale))
                
        
    def bam_to_bw(self, norm=1000000):
        """
        Create bigWig
        bam -> bigWig
        """
        args_local = {
            'bam': self.bam,
            'scaleFactor': self.align_scale,
            'outdir': self.bw_dir,
            'genome': self.genome,
            'strandness': 12,
            'binsize': self.binsize,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size
        }
        Bam2bw(**args_local).run()


    def bam_to_bg(self):
        """
        Convert bam to bedgraph

        norm
        bedtools genomecov -bg -scale $scale_factor -ibam bam > bg
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedtools')),
            'genomecov -bg -scale {}'.format(self.align_scale),
            '-ibam {}'.format(self.bam),
            '| sort -k1,1 -k2,2n > {}'.format(self.bg)
        ])

        if file_exists(self.bg) and not self.overwrite:
            log.info('bam_to_bg() skipped, file exists:{}'.format(self.bg))
        else:
            cmd_shell = self.bg_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bam_to_bg() failed')


    def bg_to_bw(self):
        """
        Create bigWig
        bedgraph to bigWig
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedGraphToBigWig')),
            '{} {} {}'.format(self.bg, self.genome_size_file, self.bw)
        ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_shell = self.bw_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor

        Bam().count_reads()
        """
        bam = Bam(bam)
        is_pe = bam.isPaired()
        n_map = bam.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    def qc_bam_cor(self, window=500):
        """
        Compute correlation (pearson) between replicates
        window = 500bp

        eg:
        multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
            --outRawCounts *counts.tab -b bam
        """
        args_local = {
            'bam': self.get_bam_list(),
            'outdir': self.qc_dir,
            'prefix': '06.bam_cor',
            'threads': self.threads,
            'overwrite': self.overwrite,
            'binsize': self.binsize
        }
        Bam2cor(**args_local).run()


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        bw_list = self.get_bw_list()
        bw_list.append(self.bw)
        bw_list_arg = ' '.join(bw_list)

        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(bw_list_arg),
            '-o {}'.format(self.genebody_enrich_matrix),
            '-b 2000 -a 2000 --regionBodyLength 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.genebody_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.genebody_enrich_matrix),
            '-o {}'.format(self.genebody_enrich_png),
            '--dpi 300',
            '--perGroup'
        ])

        if file_exists(self.genebody_enrich_png) and not self.overwrite:
            log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
                self.genebody_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.genebody_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_genebody_enrich() failed, see: {}'.format(
                        self.genebody_enrich_matrix_log))


    def report(self):
        """
        Create report for one sample
        html
        """
        RNAseqRp(self.project_dir).report()


    def run_RNAseqR1(self, i):
        """
        Run RNAseqR1 for each fastq
        """
        req_args = [
            'outdir',
            'aligner', 
            'genome', 
            'genome_index',
            'spikein',
            'spikein_index',
            'rRNA_index',
            'align_to_rRNA',
            'extra_index', 
            'gtf',
            'gene_bed',
            'trimmed', 
            'keep_tmp',
            'overwrite', 
            'threads', 
            'parallel_jobs', 
            'cut_to_length', 
            'recursive',
            'extra_para',
        ]
        args_local = dict((k, getattr(self, k, None)) for k in req_args)

        # update: fq1, fq2
        args_local.update({
            'fq1': self.fq1[i],
            'fq2': self.fq2[i] if isinstance(self.fq2, list) else None
        })
        RNAseqR1(**args_local).run()


    def pipe_rn(self):
        """Run pipeline, replicates
        data
        """
        # run
        self.merge_bam()
        if self.genome:
            self.bam_to_bw()
        else:
            self.bam_to_bg()
            self.bg_to_bw()
        self.qc_bam_cor()
        self.qc_genebody_enrich()
        self.report()


    def pipe_r1(self):
        """Run pipeline for single fastq
        Do not run merge
        """
        r1 = RNAseqReader(self.rep_list[0]).args
        bam = rep.get('bam', None)
        bg = rep.get('bg', None)
        bw = rep.get('bw', None)
        qc_dir = r1.get('qc_dir', None)

        file_symlink(bam, self.bam)
        file_symlink(bg, self.bg)
        file_symlink(bw, self.bw)
        file_symlink(qc_dir, self.qc_dir)
        self.report()

        # align_scale_txt = rep.get('align_scale_txt', None)
        # align_json = rep.get('align_json')
        # file_symlink(align_scale_txt, self.align_scale_txt)
        # file_symlink(align_json, self.align_json)
        # file_symlink(bam, self.bam)
        # Bam(self.bam).index()
        # file_symlink(bg, self.bg)
        # file_symlink(bw, self.bw)
        # # self.qc_genebody_enrich()
        # self.report()


    def run(self):
        """
        Run multiple samples in parallel
        using
        multipleprocess.Pool
        """
        print('!AAAA-6')
        # run for fastq
        if isinstance(self.fq1, list):
            if self.parallel_jobs > 1 and len(self.fq1) > 1:
                with Pool(processes=self.parallel_jobs) as pool:
                    pool.map(self.run_RNAseqR1, range(len(self.fq1)))
            else:
                for i in range(len(self.fq1)):
                    self.run_RNAseqR1(i)

        # mrege bam
        if len(self.rep_list) > 1:
            self.pipe_rn()
        elif len(self.rep_list) == 1:
            log.warning('merge() skipped, Only 1 replicate detected')
            self.pipe_r1()
        else:
            raise ValueError('merge() failed, no rep detected')


class RNAseqR1(object):
    """
    RNAseq analysis for single replicate
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        obj_local = RNAseqR1Config(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)
        Toml(self.__dict__).to_toml(self.config_toml)


    # main pipeline #
    def prep_raw(self):
        """
        Cut the reads to specific length, from the 3' end
        self.fq1, self.fq2 => raw_fq_list
        """
        raw_fq1, raw_fq2 = self.raw_fq_list
        file_symlink(self.fq1, raw_fq1, absolute_path=True)
        file_symlink(self.fq2, raw_fq2, absolute_path=True)


    def trim(self):
        """
        Trim reads:
        Trim(fq1, outdir, fq2, cut_after_trim='9,-6').run()

        if trimmed:
            do
        else:
            copy/links
        """
        fq1, fq2 = self.raw_fq_list
        clean_fq1, clean_fq2 = self.clean_fq_list

        # whether to trim or not
        if self.trimmed:
            ## copy files
            file_symlink(fq1, clean_fq1)
            file_symlink(fq2, clean_fq2)
        else:
            # update args
            args_local = self.__dict__.copy()
            args_init = {
                'fq1': fq1,
                'fq2': fq2,
                'outdir': self.clean_dir,
                'library_type': None,
                'len_min': 20,
                'cut_after_trim': '7,-7',
                'cut_to_length': self.cut_to_length,
                'recursive': self.recursive,
                'parallel_jobs': 1 # do not allowed > 1
            }
            args_local.update(args_init)
            trim = Trim(**args_local)
            trim.run()

            ## copy files
            if isinstance(fq2, str):
                file_symlink(trim.clean_fq1, clean_fq1)
                file_symlink(trim.clean_fq2, clean_fq2)
            else:
                file_symlink(trim.clean_fq, clean_fq1)


    def align_rRNA(self):
        """
        Alignment reads to rRNA
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        # !!!! force: bowtie2 !!!!
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.rRNA_dir,
            'aligner': AlignIndex().get_aligner(index=self.rRNA_index),
            'genome_index': self.rRNA_index,
            'threads': self.threads,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        # args_local.update(args_sp)
        args_local = args_sp
        align = Align(**args_local)

        # output
        tmp0 = listfile(self.rRNA_dir, '*.bam', recursive=True)
        bam = tmp[0] if len(tmp0) > 0 else None
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            align.run()

        # copy bam to bam_files
        file_symlink(bam, self.rRNA_bam)


    def align_spikein(self):
        """
        Alignment reads, spikein
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_sp = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.spikein_dir,
            'genome_index': self.spikein_index,
            'threads': self.threads,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        args_local.update(args_sp)
        align = Align(**args_local)

        # output
        tmp0 = listfile(self.spikein_dir, '*.bam', recursive=True)
        bam = tmp[0] if len(tmp0) > 0 else None
        if file_exists(bam) and not self.overwrite:
            log.info('align() skipped, file exists: {}'.format(bam))
        else:
            align.run()

        # copy bam to bam_files
        file_symlink(bam, self.spikein_bam)

        self.spikein_scale = self.cal_norm_scale(self.spikein_bam)
        with open(self.spikein_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.spikein_scale))


    def align_genome(self):
        """
        Alignment PE reads to reference genome, using bowtie2
        """
        args_local = self.__dict__.copy()
        fq1, fq2 = args_local.get('clean_fq_list', [None, None])

        # update
        args_init = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': self.align_dir,
            'aligner': self.aligner,
            'genome_index': self.genome_index,
            'extra_para': self.extra_para,
            'keep_tmp': self.keep_tmp,
            'overwrite': self.overwrite
        }
        args_local.update(args_init)
        # fix, extra_index alignment
        if args_local.get('genome_index', None) == args_local.get('extra_index', None):
            args_local['extra_index'] = None
        # align = Alignment(**args_local)
        align = Align(**args_local)
        align.run()

        # output
        tmp0 = listfile(self.align_dir, '*.bam', recursive=True)
        bam = tmp0[0] if len(tmp0) > 0 else None
        file_symlink(bam, self.bam)

        # align.toml
        tmp1 = listfile(self.align_dir, '*align.toml', recursive=True)
        if len(tmp1) > 0:
            file_symlink(tmp1[0], self.align_toml)

        self.align_scale = self.cal_norm_scale(self.bam)
        with open(self.align_scale_txt, 'wt') as w:
            w.write('{:.4f}\n'.format(self.align_scale))


    def get_bam(self, x):
        """
        Get the align bam file, from x dicectory
        align_dir/A, B, C, ...
        Add prefix/ 01, 02, ...
        """
        bamlist = listfile(x, '*.bam', recursive=True)
        return bamlist.pop() if len(bamlist) > 0 else None # last one


    def bam_to_bw(self):
        """
        Convert bam to bigWig
        """
        args = {
            'bam': self.bam,
            'outdir': self.bw_dir,
            'binsize': self.binsize,
            'strandness': 12,
            'genome': self.genome,
            'scaleFactor': self.align_scale,
            'overwrite': self.overwrite,
            'genome_size': self.genome_size,
        }
        Bam2bw(**args).run()
    
    
    def bam_to_bg(self):
        """
        Convert bam to bedgraph

        norm
        bedtools genomecov -bg -scale $scale_factor -ibam bam > bg
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedtools')),
            'genomecov -bg -scale {}'.format(self.align_scale),
            '-ibam {}'.format(self.bam),
            '| sort -k1,1 -k2,2n > {}'.format(self.bg)
            ])

        if file_exists(self.bg) and not self.overwrite:
            log.info('bam_to_bg() skipped, file exists:{}'.format(self.bg))
        else:
            cmd_shell = self.bg_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bam_to_bg() failed')


    def bg_to_bw(self):
        """
        Create bigWig
        bedgraph to bigWig
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('bedGraphToBigWig')),
            '{} {} {}'.format(self.bg, self.genome_size_file, self.bw)
            ])

        if file_exists(self.bw) and not self.overwrite:
            log.info('bg_to_bw() skipped, file exists:{}'.format(self.bw))
        else:
            cmd_shell = self.bw_dir + '/cmd.sh'
            with open(cmd_shell, 'wt') as w:
                w.write(cmd + '\n')

            try:
                run_shell_cmd(cmd)
            except:
                log.error('bg_to_bw() failed')


    def fc_count(self):
        """
        Run FeatureCounts for the bam file
        """
        args = self.__dict__.copy()

        # check library type
        s = HiSeqLibType(bam=self.bam, gtf=self.gtf).run()

        if s == 1:
            strand_sens = 1
            strand_anti = 2
        elif s == 2:
            strand_sens = 2
            strand_anti = 1
        else:
            strand_sens = strand_anti = 0

        # args for local run
        args_local = {
            'gtf': self.gtf,
            'bam_list': self.bam,
            'outdir': self.count_dir,
            'threads': self.threads,
            'overwrite': self.overwrite
        }

        # for sense
        args_sens = {
            'prefix': os.path.basename(self.count_sens),
            'strandness': strand_sens
        }
        args_sens.update(args_local)
        fc_sens = FeatureCounts(**args_sens)
        fc_sens.run()

        # for antisense
        args_anti = {
            'prefix': os.path.basename(self.count_anti),
            'strandness': strand_anti
        }
        args_anti.update(args_local)
        fc_anti = FeatureCounts(**args_anti)
        fc_anti.run()

        # save strandness to file
        msg = '\n'.join([
            'sens: {}'.format(strand_sens),
            'anti: {}'.format(strand_anti)
        ])
        with open(self.strandness_file, 'wt') as w:
            w.write(msg + '\n')


    # quality wildtype #
    def cal_norm_scale(self, bam, norm=1000000):
        """
        scale factor

        Bam().count_reads()
        """
        bam_o = Bam(bam)
        is_pe = bam_o.isPaired()
        n_map = bam_o.getNumberOfAlignments()
        if is_pe:
            n_map = n_map/2

        if n_map > 0:
            n_scale = norm/n_map
        else:
            log.error('no mapped reads detected')
            n_scale = 1

        return n_scale


    def qc_trim_summary(self):
        """
        reads trim off

        # version-1
        #sample wildtype   output  percent
        fq_rep1      2234501 2234276 99.99%

        # version-2
        #name   total   too_short       dup     too_short2      clean   percent
        """
        if file_exists(self.trim_stat_txt):
            with open(self.trim_stat_txt, 'rt') as r:
                lines = r.readlines()
            lines = [i for i in lines if not i.startswith('#')]
            line = lines.pop() # last line
            tabs = line.strip().split('\t')
            fqname, n_wildtype = tabs[:2]
            n_output, n_pct = tabs[-2:]
            d = {
                'id': fqname,
                'input': n_wildtype,
                'output': n_output,
                'out_pct': n_pct,
                'rm_pct': round(100 - float(n_pct.strip()), 2)
            }
        else:
            d = {
                'id': self.smp_name,
                'input': 0,
                'output': 0,
                'out_pct': 100,
                'rm_pct': 0
            }
        Toml(d).to_toml(self.trim_summary_toml)


    def qc_align_summary(self):
        """
        Save summary to file
        """
        # align to genome
        align = Toml().from_toml(self.align_toml)

        # spikein
        if file_exists(self.spikein_bam):
            n_spikein = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_spikein = n_spikein/2
        else:
            n_spikein = 0

        # rRNA
        if file_exists(self.rRNA_bam):
            n_rRNA = Bam(self.spikein_bam).getNumberOfAlignments()
            if is_pe:
                n_rRNA = n_rRNA/2
        else:
            n_rRNA = 0

        align['nodup'] = 0
        align['spikein'] = n_spikein
        align['chrM'] = n_rRNA

        # save to file
        Toml(align).to_toml(self.align_summary_toml)


    def qc_genebody_enrich(self):
        """
        Calculate the TSS enrichment
        """
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            'scale-regions',
            '-R {}'.format(self.gene_bed),
            '-S {}'.format(self.bw),
            '-o {}'.format(self.genebody_enrich_matrix),
            '-b 2000 -a 2000 --regionBodyLength 2000',
            '--binSize 10 --sortRegions descend --skipZeros',
            '--smartLabels',
            '-p {}'.format(self.threads),
            '2> {}'.format(self.genebody_enrich_matrix_log),
            '&& {}'.format(shutil.which('plotProfile')),
            '-m {}'.format(self.genebody_enrich_matrix),
            '-o {}'.format(self.genebody_enrich_png),
            '--dpi 300',
            '--perGroup'
            ])

        if file_exists(self.genebody_enrich_png) and not self.overwrite:
            log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
                self.genebody_enrich_png))
        else:
            if not file_exists(self.gene_bed):
                log.error('qc_tss() skipped, gene_bed not found')
            else:
                with open(self.genebody_enrich_cmd, 'wt') as w:
                    w.write(cmd + '\n')

                try:
                    run_shell_cmd(cmd)
                except:
                    log.error('qc_genebody_enrich() failed, see: ')


    def pipe_genome(self):
        """
        Run for genome
        """
        self.prep_raw()
        self.trim()

        if isinstance(self.spikein_index, str):
            self.align_spikein() # run alignment

        if isinstance(self.rRNA_index, str):
            self.align_rRNA()

        self.align_genome()
        if self.genome:
            self.bam_to_bw()
        else:
            self.bam_to_bg()
            self.bg_to_bw()

        # count reads
        if file_exists(self.count_sens) and not self.overwrite:
            log.info('fc_count() skipped, file exists: {}'.format(
                self.count_sens))
        else:
            self.fc_count()


    def pipe_qc(self):
        """
        Quality control
        """
        self.qc_trim_summary()
        self.qc_align_summary()
        self.qc_genebody_enrich()


    def report(self):
        """
        Create report for one sample
        html
        """
        RNAseqRp(self.project_dir).report()


    def run(self):
        """
        Run all steps for RNAseq pipeline
        """
        self.pipe_genome()
        self.pipe_qc()
        self.report()


class RNAseqRp(object):
    """
    Run report for RNAseq dirs: R1, Rn, Rx, ...
    input: project_dir
    output: project_dir/report/
    """
    def __init__(self, project_dir=None, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.project_dir = project_dir
        self.init_args()


    def init_args(self):
        obj_local = RNAseqRpConfig(**self.__dict__)
        self = update_obj(self, obj_local.__dict__, force=True)


    def report(self):
        """
        Create report for one sample
        html
        """
        pkg_dir = os.path.dirname(hiseq.__file__)
        qc_reportR = os.path.join(pkg_dir, 'bin', 'hiseq_report.R')
        cmd = 'Rscript {} {} {} 1> {} 2> {}'.format(
            qc_reportR,
            self.project_dir,
            self.report_dir,
            self.report_stdout,
            self.report_stderr
        )
        cmd_shell = os.path.join(self.report_dir, 'cmd.sh')
        with open(cmd_shell, 'wt') as w:
            w.write(cmd + '\n')
        if file_exists(self.report_html):
            log.info('report() skipped, file exists: {}'.format(
                self.report_html))
        else:
            try:
                run_shell_cmd(cmd)
            except:
                log.error('report() failed, {}'.format(self.report_dir))


    def run(self):
        """
        Create report for multiple samples (n groups)
        """
        self.report()


class RNAseqConfig(object):
    """
    The main port, for RNAseq analysis

    groups = {RNAseqRx}n = {RNAseqRn, RNAseqRn}n = {RNAseqR1}n
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseq analysis
        """         
        args_init = {
            'build_design': False,
            'design': None,
            'rep_list': None,
            'mutant_dir': None,
            'wildtype_dir': None,
            'mutant': None,
            'mutant_fq2': None,
            'wildtype': None,
            'wildtype_fq2': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'aligner': 'STAR',
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'rRNA_index': None,
            'align_to_rRNA': False,
            'extra_index': None,
            'threads': 1,
            'parallel_jobs': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'genome_size_file': None,
            'gene_bed': None,
            'gtf': None,
            'trimmed': False,
            'keep_tmp': False,
            'cut_to_length': 0,
            'recursive': False,
            'extra_para': None, # for alignment
        }
        self = update_obj(self, args_init, force=False)

        # output
        if self.outdir is None:
          self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir)

        # aligner
        if self.aligner is None:
            self.aligner = 'STAR'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)

        # RNAseq groups
        self.fq_groups = Config().from_toml(self.design) if \
            file_exists(self.design) else self.init_groups()

        self.init_files()
        self.init_mission()


    def init_files(self):
        """
        Prepare directories, files
        """
        self.project_dir = self.outdir
        self.config_dir = self.outdir
        self.config_toml = self.config_dir + '/config.toml'
        check_path(self.project_dir)


    def init_groups(self):
        """
        Create RNAseq compare group:
        - wildtype, wildtype_fq2
        - mutant, mutant_fq2

        - wildtype_dir
        - mutant_dir
        """
        if isinstance(self.wildtype, list) and isinstance(self.mutant, list):
            g = {
                'wildtype_dir': None,
                'mutant_dir': None,
                'wildtype': self.wildtype,
                'wildtype_fq2': self.wildtype_fq2,
                'mutant': self.mutant,
                'mutant_fq2': self.mutant_fq2,
                'group': fq_name_rmrep(self.wildtype).pop()
            }
        elif isinstance(self.wildtype_dir, str) \
            and isinstance(self.mutant_dir, str):
            c_wildtype = RNAseqReader(self.wildtype_dir)
            c_mutant = RNAseqReader(self.mutant_dir)

            if not c_wildtype.is_hiseq_rn or not c_mutant.is_hiseq_rn:
                raise ValueError('mutant_dir, wildtype_dir failed: \
                    {}, {}'.format(self.mutant_dir, self.wildtype_dir))
            g = {
                'wildtype_dir': self.wildtype_dir,
                'mutant_dir': self.mutant_dir,
                'wildtype': None,
                'wildtype_fq2': None,
                'mutant': None,
                'wildtype_fq2': None,
                'group': c_wildtype.args.get('smp_name', 'RNAseq')
            }
        else:
            g = None

        return {'RNAseq': g} if g else {} # empty


    def init_mission(self):
        """
        Determine the type of RNAseq analysis
        1. build_design
        2. Single: mutant/wildtype pair
        3. Multiple: mutant/wildtype pair
        4. x, mutant vs wildtype
        """
        # 1st level
        if self.build_design:
            print('!AAAA-1', 'build_design')
            self.hiseq_type = 'build_design'
        elif len(self.fq_groups) >= 1:
            print('!AAAA-2', 'rnaseq_rx')
            self.hiseq_type = 'rnaseq_rx'
        elif isinstance(self.fq1, list):
            print('!AAAA-3', 'rnaseq_rn')
            self.hiseq_type = 'rnaseq_rn'
        elif isinstance(self.fq1, str):
            print('!AAAA-4', 'rnaseq_r1')
            self.hiseq_type = 'rnaseq_r1'
        else:
            raise ValueError('unknown RNAseq type')


class RNAseqRxConfig(object):
    """
    RNAseq for mutant vs wildtype: DESeq analysis

    require: mutant_dir, wildtype_dir
    or:
    mutant, mutant_fq2, wildtype, wildtype_fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseq analysis
        """
        args_init = {
            'mutant_dir': None,
            'wildtype_dir': None,
            'mutant': None,
            'mutant_fq2': None,
            'wildtype': None,
            'wildtype_fq2': None,
            'outdir': None,
        }
        # args_init = {
        #     'outdir': None,
        #     'mutant_dir': None,
        #     'wildtype_dir': None,
        #     'mutant': None,
        #     'mutant_fq2': None,
        #     'wildtype': None,
        #     'wildtype_fq2': None,
        #     'smp_name': None,
        #     'aligner': 'STAR',
        #     'genome': None,
        #     'genome_index': None,
        #     'spikein': None,
        #     'spikein_index': None,
        #     'rRNA_index': None,
        #     'align_to_rRNA': False,
        #     'extra_index': None,
        #     'threads': 1,
        #     'parallel_jobs': 1,
        #     'overwrite': False,
        #     'binsize': 50,
        #     'genome_size': 0,
        #     'genome_size_file': None,
        #     'gene_bed': None,
        #     'keep_tmp': False,
        #     'trimmed': False,
        #     'cut_to_length': 0,
        #     'recursive': False,
        #     'extra_para': None
        # }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rx' #

        # aligner
        if self.aligner is None:
            self.aligner = 'STAR'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)

        self.init_fq()
        self.init_files()


    def init_fq(self):
        """
        Support fastq files
        mutant, wildtype
        fq1, list
        fq2, list or None
        """
        # fastq files
        if isinstance(self.mutant, list) and isinstance(self.wildtype, list):
            self.mutant = file_abspath(self.mutant)
            self.mutant_fq2 = file_abspath(self.mutant_fq2)
            self.wildtype = file_abspath(self.wildtype)
            self.wildtype_fq2 = file_abspath(self.wildtype_fq2)

            # file exists
            if not all(file_exists(self.mutant)):
                raise ValueError('--mutant, file not exists: {}'.format(
                    self.mutant))

            if not all(file_exists(self.wildtype)):
                raise ValueError('--wildtype, file not exists: {}'.format(
                    self.mutant))

            # check, mutant,wildtype name
            mutant_names = set(map(fq_name_rmrep, self.mutant))
            wildtype_names = set(map(fq_name_rmrep, self.wildtype))
            if len(mutant_names) > 1:
                raise ValueError('--mutant failed, filename \
                    differ: {}'.format(mutant_names))

            if len(wildtype_names) > 1:
                raise ValueError('--wildtype failed, filename \
                    differ: {}'.format(wildtype_names))

            # check paired
            if isinstance(self.mutant_fq2, list):
                mutant_pe = all([fq_paired(i, j) for i, j \
                    in zip(self.mutant, self.mutant_fq2)])
                if not mutant_pe:
                    raise ValueError('--mutant, --mutant-fq2, not paired: \
                        {}, {}'.format(self.mutant, self.mutant_fq2))
            if isinstance(self.wildtype_fq2, list):
                wildtype_pe = all([fq_paired(i, j) for i, j \
                    in zip(self.wildtype, self.wildtype_fq2)])
                if not wildtype_pe:
                    raise ValueError('--wildtype, --wildtype-fq2, not paired: \
                        {}, {}'.format(self.wildtype, self.wildtype_fq2))

            # update mutant_dir, wildtype_dir
            self.mutant_name = mutant_names.pop()
            self.wildtype_name = wildtype_names.pop()
            self.mutant_dir = self.outdir + '/' + self.mutant_name
            self.wildtype_dir = self.outdir + '/' + self.wildtype_name

        elif isinstance(self.mutant_dir, str) \
            and isinstance(self.wildtype_dir, str):
            self.mutant_dir = file_abspath(self.mutant_dir)
            self.wildtype_dir = file_abspath(self.wildtype_dir)

            c_mutant = RNAseqReader(self.mutant_dir)
            c_wildtype = RNAseqReader(self.wildtype_dir)

            if not c_mutant.is_hiseq_rn or not c_wildtype.is_hiseq_rn:
                raise ValueError('mutant_dir, wildtype_dir failed: \
                    {}, {}'.format(self.mutant_dir, self.wildtype_dir))

            self.mutant_name = c_mutant.args.get('smp_name', None)
            self.wildtype_name = c_wildtype.args.get('smp_name', None)

        else:
            raise ValueError('--mutant, --wildtype, or --mutant-dir, \
                --wildtype-dir failed')

        # update smp_name
        self.smp_name = '{}.vs.{}'.format(self.wildtype_name, self.mutant_name)


    def init_files(self):
        """
        Prepare directories, files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = self.outdir + '/' + self.project_name
        self.config_dir = self.project_dir + '/config'

        default_dirs = {
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'count_dir': 'count',
            'deseq_dir': 'deseq',
            'enrich_dir': 'enrich',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml',
            'mutant_bam': self.bam_dir + '/' + self.mutant_name + '.bam',
            'mutant_bw': self.bw_dir + '/' + self.mutant_name + '.bigWig',
            'wildtype_bam': self.bam_dir + '/' + self.wildtype_name + '.bam',
            'wildtype_bw': self.bw_dir + '/' + self.wildtype_name + '.bigWig',
            'bw': self.bw_dir + '/' + self.mutant_name + '.bigWig',
            # 'mutant_count_sens': RNAseqReader(
            #     self.mutant_dir).get_r1_file('count_sens'),
            # 'mutant_count_anti': RNAseqReader(
            #     self.mutant_dir).get_r1_file('count_anti'),
            # 'wildtype_count_sens': RNAseqReader(
            #     self.wildtype_dir).get_r1_file('count_sens'),
            # 'wildtype_count_anti': RNAseqReader(
            #     self.wildtype_dir).get_r1_file('count_anti'),
            'deseq_fix_xls': self.deseq_dir + '/transcripts_deseq2.fix.xls'
        }
        self = update_obj(self, default_files, force=True) # key

        check_path([
            self.project_dir,
            self.config_dir,
            self.bam_dir,
            self.bw_dir,
            self.count_dir,
            self.deseq_dir,
            self.enrich_dir,
            self.qc_dir,
            self.report_dir])


class RNAseqRnConfig(object):
    """
    RNAseq for single sample, N replicates; wrap replicates
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseq analysis
        """
        args_init = {
            'rep_list': None,
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'smp_name': None
        }

        # args_init = {
        #     'rep_list': None,
        #     'fq1': None,
        #     'fq2': None,
        #     'outdir': None,
        #     'aligner': 'STAR',
        #     'smp_name': None,
        #     'genome': None,
        #     'genome_index': None,
        #     'spikein': None,
        #     'spikein_index': None,
        #     'rRNA_index': None,
        #     'align_to_rRNA': False,
        #     'extra_index': None,
        #     'threads': 1,
        #     'parallel_jobs': 1,
        #     'overwrite': False,
        #     'binsize': 50,
        #     'genome_size': 0,
        #     'genome_size_file': None,
        #     'keep_tmp': False,
        #     'trimmed': False,
        #     'cut_to_length': 0,
        #     'recursive': False
        # }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_rn' #

        # aligner
        if self.aligner is None:
            self.aligner = 'STAR'

        # threads
        self.threads, self.parallel_jobs = init_cpu(self.threads,
            self.parallel_jobs)

        self.outdir = file_abspath(self.outdir)

        self.init_fq()
        self.init_files()


    def init_fq(self):
        """
        Support fastq files
        fq1 (required), list
        fq2 (optional), list or None
        """
        # convert str to list
        if isinstance(self.fq1, str):
            self.fq1 = [self.fq1]

        if isinstance(self.fq2, str):
            self.fq2 = [self.fq2]

        if isinstance(self.fq1, list):
            self.fq1 = file_abspath(self.fq1)

            # file exists
            if not all(file_exists(self.fq1)):
                raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

            # rep list
            q_names = [fq_name(i, pe_fix=True) for i in self.fq1]
            self.rep_list = [self.outdir + '/' + i for i in q_names]

            # fq2
            if isinstance(self.fq2, list):
                self.fq2 = file_abspath(self.fq2)

                if not all(file_exists(self.fq2)):
                    raise ValueError('--fq2, file not exists: {}'.format(
                        self.fq2))

                # paired
                fq_pe = all([fq_paired(i, j) \
                    for i, j in zip(self.fq1, self.fq2)])
                
                if not fq_pe and isinstance(self.fq2, list):
                    raise ValueError('--fq1, --fq2, file not paired')
        else:
            raise ValueError('fq1, fq2 required')

        # smp_name
        if not isinstance(self.smp_name, str):
            self.smp_name = fq_name_rmrep(self.rep_list).pop()


    def init_files(self):
        """
        Prepare directories, files
        for bigWig files, track files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        default_dirs = {
            'align_dir': 'align',
            'bam_dir': 'bam_files',
            'bw_dir': 'bw_files',
            'bg_dir': 'bg_files',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bed': self.bam_dir + '/' + self.project_name + '.bed',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'bw_fwd': self.bw_dir + '/' + self.project_name + '.fwd.bigWig',
            'bw_rev': self.bw_dir + '/' + self.project_name + '.rev.bigWig',
            'align_scale_txt': self.align_dir + '/' + 'scale.txt',
            'genebody_enrich_matrix': self.qc_dir \
                + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir \
                + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png'
        }
        self = update_obj(self, default_files, force=True) # key

        check_path([
            self.config_dir,
            self.align_dir,
            self.bam_dir,
            self.bw_dir,
            self.bg_dir,
            self.qc_dir,
            self.report_dir])


class RNAseqR1Config(object):
    """
    RNAseq for single replicate, single index
    input:
    - fq1, fq2
    - genome, genome_index
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseq analysis
        """
        args_init = {
            'fq1': None,
            'fq2': None,
            'outdir': None,
            'smp_name': None,
            'aligner': 'STAR',
            'genome': None,
            'genome_index': None,
            'spikein': None,
            'spikein_index': None,
            'rRNA_index': None,
            'align_to_rRNA': False,
            'extra_index': None,
            'threads': 1,
            'overwrite': False,
            'binsize': 50,
            'genome_size': 0,
            'genome_size_file': None,
            'gene_bed': None,
            'gtf': None,
            'keep_tmp': False,
            'trimmed': False,
            'cut_to_length': 0,
            'recursive': False,
            'extra_para': None
        }
        self = update_obj(self, args_init, force=False)
        self.hiseq_type = 'rnaseq_r1' #

        # aligner
        if self.aligner is None:
            self.aligner = 'STAR'

        # smp_name
        if not isinstance(self.smp_name, str):
            self.smp_name = fq_name(self.fq1, pe_fix=True) # the first one

        # output
        self.outdir = file_abspath(self.outdir)

        self.init_fq()
        self.init_index()
        self.init_files()


    def init_fq(self):
        """
        Support fastq files
        fq1 (required), str
        fq2 (optional), str or None
        """
        # fq1
        if not isinstance(self.fq1, str):
            raise ValueError('--fq1, str expected, got {}'.format(
                type(self.fq1).__name__))

        # abs
        self.fq1 = file_abspath(self.fq1)

        # file exists
        if not file_exists(self.fq1):
            raise ValueError('--fq1, file not exists: {}'.format(self.fq1))

        # fq2
        if not self.fq2 is None:
            if not isinstance(self.fq2, str):
                raise ValueError('--fq2, None or str expected, got {}'.format(
                    type(self.fq2).__name__))

            # abs
            self.fq2 = file_abspath(self.fq2)

            if not file_exists(self.fq2):
                raise ValueError('--fq2, file not exists: {}'.format(self.fq2))

            # paired
            if not fq_paired(self.fq1, self.fq2):
                raise ValueError('--fq1, --fq2, file not paired')


    def init_index(self):
        """Prepare index for alignment
        1. genome_index: genome, extra_index
        2. rRNA_index: align_to_rRNA, ...
        3. spikein_index: spikein, spikein_index
        """
        # extra index
        if isinstance(self.extra_index, str):
            self.genome_index = self.extra_index
        else:
            # genome, spikein, rRNA, MT, ...
            # 1. genome_index
            if isinstance(self.genome_index, str):
                ai = AlignIndex(index=self.genome_index)

                # valid
                if not ai.is_index():
                    raise ValueError('genome_index failed, {}'.format(
                        self.genome_index))

                # chrsizes
                if self.genome_size < 1:
                    self.genome_size = ai.index_size()

                if not isinstance(self.genome_size_file, str): 
                    self.genome_size_file = ai.index_size(return_file=True)

            elif isinstance(self.genome, str):
                self.genome_index = AlignIndex(aligner=self.aligner).search(
                    genome=self.genome, group='genome')

                self.genome_size_file = Genome(genome=self.genome).get_fasize()

                if self.genome_size < 1:
                    with open(self.genome_size_file, 'rt') as r:
                        s = [i.strip().split('\t')[1] for i in r.readlines()]
 
                    self.genome_size = sum(map(int, s))
            
            else:
                self.genome_index = None

        # 2. rRNA_index:
        if isinstance(self.rRNA_index, str):
            ai = AlignIndex(aligner=self.aligner, index=self.rRNA_index)
            if not ai.is_index():
                raise ValueError('rRNA_index failed, {}'.format(
                    self.rRNA_index))
        elif self.align_to_rRNA:
            # !!!! force: bowtie2 for rRNA !!!!
            if isinstance(self.genome, str):
                # ri = AlignIndex(aligner=self.aligner)
                ri = AlignIndex(aligner='bowtie2')
                for rRNA in ['rRNA', 'MT_trRNA', 'trRNA']:
                    ri_index = ri.search(genome=self.genome, group=rRNA)
                    if ri_index:
                        log.info('fetch rRNA index: {}, {}'.format(
                            rRNA, ri_index))
                        self.rRNA_index = ri_index
                        break
        else:
            self.rRNA_index = None

        # 3. spikein_index
        if isinstance(self.spikein_index, str):
            ai = AlignIndex(aligner=self.aligner, index=self.spikein_index)
            if not ai.is_index():
                raise ValueError('spikein_index failed, {}'.format(
                    self.spikein_index))
        elif isinstance(self.spikein, str):
            self.spikein_index = AlignIndex(aligner=self.aligner).search(
                genome=self.spikein, group='genome')

        else:
            self.spikein_index = None

        # summarize all index
        self.spikein_index = getattr(self, 'spikein_index', None)
        self.extra_index = getattr(self, 'extra_index', [])
        self.tags_index = getattr(self, 'tags_index', [])
        self.genome_index = getattr(self, 'genome_index', [])

        # check
        ai = AlignIndex(index=self.genome_index)
        if not ai.is_index():
            raise ValueError('genome_index failed, {}'.format(
                self.genome_index))


    def init_files(self):
        """
        Prepare directories, files
        """
        # path, files
        self.project_name = self.smp_name
        self.project_dir = os.path.join(self.outdir, self.project_name)
        self.config_dir = os.path.join(self.project_dir, 'config')

        default_dirs = {
            'raw_dir': 'raw_data',
            'clean_dir': 'clean_data',
            'align_dir': 'align',
            'spikein_dir': 'spikein',
            'rRNA_dir': 'rRNA',
            'bam_dir': 'bam_files',
            'bg_dir': 'bg_files',
            'bw_dir': 'bw_files',
            'count_dir': 'count',
            'qc_dir': 'qc',
            'report_dir': 'report'
        }
        # convert to path
        for k, v in default_dirs.items():
            default_dirs[k] = os.path.join(self.project_dir, v)
        self = update_obj(self, default_dirs, force=True) # key

        # files
        default_files = {
            'config_toml': self.config_dir + '/config.toml',
            'bam': self.bam_dir + '/' + self.project_name + '.bam',
            'bg': self.bg_dir + '/' + self.project_name + '.bedGraph',
            'bw': self.bw_dir + '/' + self.project_name + '.bigWig',
            'bw_fwd': self.bw_dir + '/' + self.project_name + '.fwd.bigWig',
            'bw_rev': self.bw_dir + '/' + self.project_name + '.rev.bigWig',
            'strandness_file': self.count_dir + '/strandness.txt',
            'count_sens': self.count_dir \
                + '/' + self.project_name + '.sens.txt',
            'count_anti': self.count_dir \
                + '/' + self.project_name + '.anti.txt',
            'trim_stat_txt': self.clean_dir + '/' + self.project_name \
                + '/' + self.project_name + '.trim.stat',
            'align_scale_txt': self.align_dir + '/' + 'scale.txt',
            'align_flagstat': self.align_dir + '/' \
                + self.smp_name + '.flagstat',
            'align_stat': self.align_dir \
                + '/' + self.smp_name + '.bowtie2.stat',
            'align_toml': self.align_dir \
                + '/' + self.smp_name + '.bowtie2.toml',

            'rRNA_bam': self.rRNA_dir + '/rRNA.bam',
            'rRNA_stat': self.rRNA_dir + '/rRNA.align.stat',
            'rRNA_toml': self.rRNA_dir + '/rRNA.align.toml',
            'spikein_bam': self.spikein_dir + '/spikein.bam',
            'spikein_scale_txt': self.spikein_dir + '/scale.txt',
            'spikein_flagstat': self.spikein_dir + '/spikein.flagstat',
            'spikein_stat': self.spikein_dir + '/spikein.align.stat',
            'spikein_toml': self.spikein_dir + '/spikein.align.toml',
            'trim_summary_toml': self.qc_dir + '/00.trim_summary.toml',
            'align_summary_toml': self.qc_dir + '/01.alignment_summary.toml',
            'genebody_enrich_matrix': self.qc_dir \
                + '/05.genebody_enrich.mat.gz',
            'genebody_enrich_matrix_log': self.qc_dir \
                + '/05.genebody_enrich.log',
            'genebody_enrich_png': self.qc_dir + '/05.genebody_enrich.png',
            'genebody_enrich_cmd': self.qc_dir + '/05.genebody_enrich.cmd.sh',
            'bam_cor_npz': self.qc_dir + '/06.bam_cor.npz',
            'bam_cor_counts': self.qc_dir + '/06.bam_cor.counts.tab',
            'bam_cor_heatmap_png': self.qc_dir + '/06.bam_cor.cor_heatmap.png',
            'bam_cor_pca_png': self.qc_dir + '/06.bam_cor.cor_PCA.png'
        }
        self = update_obj(self, default_files, force=True) # key

        # raw data
        self.raw_fq_list = [self.raw_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            raw_fq2 = self.raw_dir + '/' + os.path.basename(self.fq2)
        else:
            raw_fq2 = None
        self.raw_fq_list.append(raw_fq2)

        ## clean data
        self.clean_fq_list = [self.clean_dir + '/' + os.path.basename(self.fq1)]
        if isinstance(self.fq2, str):
            clean_fq2 = self.clean_dir + '/' + os.path.basename(self.fq2)
        else:
            clean_fq2 = None
        self.clean_fq_list.append(clean_fq2)

        check_path([
            self.config_dir,
            self.raw_dir,
            self.clean_dir,
            self.align_dir,
            self.spikein_dir,
            self.rRNA_dir,
            self.bam_dir,
            self.bg_dir,
            self.bw_dir,
            self.count_dir,
            self.qc_dir,
            self.report_dir])


class RNAseqRpConfig(object):
    """
    RNAseq for report
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        """
        required arguments for RNAseqRp report
        """
        # if not os.path.exists(self.project_dir):
        #     raise Exception('project_dir, dir not exists: {}'.format(
        #         self.project_dir))
        self.hiseq_type = RNAseqReader(self.project_dir).hiseq_type #
        if self.hiseq_type == None:
            raise Exception('project_dir, not a hiseq_ dir: {}'.format(
                self.project_dir))

        # path, files
        self.report_dir = os.path.join(self.project_dir, 'report')
        self.report_html = os.path.join(self.report_dir, 'HiSeq_report.html')
        self.report_stdout = os.path.join(self.report_dir, 'report.stdout')
        self.report_stderr = os.path.join(self.report_dir, 'report.stderr')
        check_path(self.report_dir)


class RNAseqDesign(object):
    """
    Prepare mutant/wildtype samples for RNAseq analysis
    append/update

    format:
    mutant: fq1, fq2
    wildtype: fq1, fq2
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'wildtype_dir': None,
            'mutant_dir': None,
            'wildtype': None,
            'wildtype_fq2': None,
            'mutant': None,
            'mutant_fq2': None,
            'design': None
        }
        self = update_obj(self, args_init, force=False)

        ## check
        self.design = file_abspath(self.design)
        if not isinstance(self.design, str):
            raise ValueError('--design, required')

        ## check fastq files
        fq_chk = []
        log.info('Check: wildtype')
        fq1 = self.wildtype
        fq2 = self.wildtype_fq2 if isinstance(self.wildtype_fq2, list) \
            else [None]*len(fq1)
        for q1, q2 in zip(fq1, fq2):
            fq_chk.append(self.fq_check(q1))
            # paired
            if isinstance(q2, str):
                fq_chk.append(self.fq_check(q2))
                fq_chk.append(self.fq_paired(q1, q2))

        log.info('Check: mutant')
        fq1 = self.mutant
        fq2 = self.mutant_fq2 if isinstance(self.mutant_fq2, list) \
            else [None]*len(fq1)
        for q1, q2 in zip(fq1, fq2):
            fq_chk.append(self.fq_check(q1))
            # paired
            if isinstance(q2, str):
                fq_chk.append(self.fq_check(q2))                
                fq_chk.append(self.fq_paired(q1, q2))

        # check
        if not all(fq_chk):
            raise ValueError('fastq files checking failed')


    def fq_check(self, x):
        """Check input file exists
        """
        chk = file_exists(x) if file_exists(x) else False
        log.info('{} : {}'.format(chk, x))
        return chk


    def fq_paired(self, fq1, fq2):
        """
        Check fq1, fq2, (list)
        """
        if isinstance(fq1, str) and isinstance(fq2, str):
            chk = distance(fq1, fq2) == 1
        else:
            chk = False

        if not chk:
            log.warning('fastq pairing failed')

        return chk


    def run(self, design=None):
        """
        Save args in Json format
        """
        d = Toml().from_toml(self.design) if file_exists(self.design) else {}
        # d = {} # fresh start

        # local
        args_local = {
            'design': self.design,
            'mutant': self.mutant,
            'mutant_fq2': self.mutant_fq2,
            'wildtype': self.wildtype,
            'wildtype_fq2': self.wildtype_fq2}
        # remove None
        args_local = {k:v for k, v in args_local.items() if not v == None}

        # saving to dict
        key = 'rnaseq_{:03d}'.format(len(d) + 1)
        if args_local in list(d.values()):
            log.warning('design exists, skipped ...')
        else:
            d.update({key:args_local})

        Toml(d).to_toml(self.design)


class RNAseqReader(object):
    """
    Read config.toml from the local directory
    """
    def __init__(self, x):
        self.x = x
        self.args = self.read()
        self.hiseq_type = self.args.get('hiseq_type', None)
        self.is_hiseq_r1 = self.hiseq_type == 'rnaseq_r1'
        self.is_hiseq_rn = self.hiseq_type == 'rnaseq_rn'
        self.is_hiseq_rx = self.hiseq_type == 'rnaseq_rx'
        self.is_hiseq_rp = self.hiseq_type == 'rnaseq_rp'
        self.is_hiseq_rd = self.hiseq_type == 'build_design'

        if isinstance(self.hiseq_type, str):
            self.is_hiseq = self.hiseq_type.startswith('hiseq_')
        else:
            self.is_hiseq = False


    def read(self, x=None):
        """
        # hiseq
        hiseq
          |-config
          |   |-config.toml

        # alignment
        align_dir
          |- smp_nmae
          |    |- index
          |    |    |- config.toml
        """
        if x is None:
            x = self.x

        if not isinstance(x, str):
            args = {}
        else:
            if os.path.isdir(x):
                p1x = os.path.join(x, 'config', 'config.toml')
                p2x = os.path.join(x, '*', '*', 'config.toml')
                p1 = glob.glob(p1x)
                p2 = glob.glob(p2x)
                # read config
                if len(p1) == 1:
                    args = Toml().from_toml(p1[0])
                elif len(p2) == 1:
                    args = Toml().from_toml(p2[0])
                else:
                    args = {}
            else:
                args = {}

        return args


    def get_r1_dir(self):
        """
        search r1 dirs, corresponding to the x (R1, Rn, Rx)
        """
        if self.hiseq_type == 'rnaseq_r1':
            r1_list = self.x
        elif self.hiseq_type == 'rnaseq_rn':
            r1_list = self.args.get('rep_list', None)
        elif self.hiseq_type == 'rnaseq_rx':
            # mut
            mut_dir = self.args.get('mutant_dir', None)
            mut_args = self.read(mut_dir)
            mut_r1 = mut_args.get('rep_list', None)

            # wt
            wt_dir = self.args.get('wildtype_dir', None)
            wt_args = self.read(wt_dir)
            wt_r1 = wt_args.get('rep_list', None)

            # output
            r1_list = []

            if mut_r1:
                r1_list.extend(mut_r1)

            if wt_r1:
                r1_list.extend(wt_r1)

        else:
            r1_list = []

        return r1_list


    def get_rn_dir(self):
        """
        search rn dirs, corresponding to the x (R1, Rn, Rx)
        """
        if self.hiseq_type == 'rnaseq_rn':
            rn_list = self.x
        elif self.hiseq_type == 'rnaseq_rx':
            mut_dir = self.args.get('mutant_dir', None)
            wt_dir = self.args.get('wildtype_dir', None)
            rn_list = [i for i in [mut_dir, wt_dir] if i]
        else:
            rn_list = []

        return rn_list


    def get_rx_dir(self):
        """
        search rx dir
        """
        return self.x if self.hiseq_type == 'rnaseq_rx' else None


    def get_r1_file(self, name='bam'):
        """
        Parsing files from RNAseqR1 directories
        """
        if self.hiseq_type == 'rnaseq_r1':
            f_list = self.args.get(name, None)
        elif self.hiseq_type == 'rnaseq_rn':
            r1_list = self.args.get('rep_list', [])
            f_list = [self.read(i).get(name, None) for i in r1_list]
        elif self.hiseq_type == 'rnaseq_rx':
            # mut
            mut_dir = self.args.get('mutant_dir', None)
            mut_args = self.read(mut_dir)
            mut_r1 = mut_args.get('rep_list', None)

            # wt
            wt_dir = self.args.get('wildtype_dir', None)
            wt_args = self.read(wt_dir)
            wt_r1 = wt_args.get('rep_list', None)

            # output
            r1_list = []

            if mut_r1:
                r1_list.extend(mut_r1)

            if wt_r1:
                r1_list.extend(wt_r1)

            f_list = [self.read(i).get(name, None) for i in r1_list]
        else:
            f_list = []

        return f_list


class HiSeqLibType(object):
    """
    Guess the library type of the HiSeq data
    Strandness

    strandness: 1 ++, 1 --, / 2 +-, 2 -+ :
    dUTP, NSR: 1 +-, 1 -+, / 2 ++, 2 -- :

    infer_experiment.py from RSeQC package.
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'bam': None,
            'gtf': None,
            'size': 200000,
            'outdir': None,
            'clean_up': True
        }
        self = update_obj(self, args_init, force=False)

        if not isinstance(self.bam, str):
            raise Exception('bam, str and list expected, got {}'.format(
                type(self.bam).__name__))

        if not isinstance(self.gtf, str):
            raise Exception('gtf, str expected, got {}'.format(
                type(self.gtf).__name__))

        # output, tmp dir
        if self.outdir is None:
            self.outdir = self._tmp(is_dir=True)
        check_path(self.outdir)


    def _tmp(self, is_dir=False, suffix='.txt'):
        """
        Create a tmp file to save json object
        """
        if is_dir:
            tmp = tempfile.TemporaryDirectory(prefix='tmp')
        else:
            tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix=suffix,
                delete=False)
        return tmp.name


    def stat(self, strandness=1):
        """
        strandness; 1 or 2; for featureCounts
        using -s 1, -s 2
        and return the mapping reads
        """
        prefix = 'count.s_{}.txt'.format(strandness)
        args = {
            'gtf': self.gtf,
            'bam_list': self.bam_sub,
            'outdir': self.outdir,
            'prefix': prefix,
            'strandness': strandness
        }
        return FeatureCounts(**args).run()
        

    def run(self, with_status=False):
        """
        Check -s 1, -s 2:
        """
        # subset BAM
        self.bam_sub = Bam(self.bam).subset(self.size, self.outdir)

        # run fc
        dfa = self.stat(strandness=1)
        dfb = self.stat(strandness=2)

        sa = dfa.pct.to_list().pop()
        sb = dfb.pct.to_list().pop()


        # check strand:
        if sa + sb < 0.2:
            log.warning('Less than 20\% reads assigned, {}'.format(self.gtf))
            s = 0
        elif sa == sb:
            log.info('Library type (s=0): non-stranded')
            s = 0
        elif sa > sb:
            log.info('Library type (s=1): 1 ++, 1 --, / 2 +-, 2 -+')
            s = 1
        else:
            log.info('Library type (s=2): 1 +-, 1 -+, / 2 ++, 2 --')
            s = 2
        
        # clean up
        if self.clean_up:
            path_remove(self.outdir, ask=False)

        # log
        msg = '\n'.join([
            'Check BAM strandness by featureCounts',
            'Library strandness: -s=1 or -s=2',
            'bam: {}'.format(self.bam),
            'gtf: {}'.format(self.gtf),
            'forward stranded, -s=1, {:>6.2f}%'.format(sa*100),
            'reverse stranded, -s=2, {:>6.2f}%'.format(sb*100)
        ])
        log.info(msg)

        # output
        return s


