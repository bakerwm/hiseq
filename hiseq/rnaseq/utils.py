#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
General modules for RNAseq analysis

analysis-module:
"""

import os
import sys
import re
import glob
import shutil
import pysam
import hiseq
import pyfastx
import hiseq # for pkg files
from hiseq.trim.trim_r1 import TrimR1
from hiseq.align.align import Align
from hiseq.bam2bw.bam2bw import Bam2bw, bw_compare
from hiseq.fragsize.fragsize import BamFragSize
from hiseq.utils.file import (
    list_file, list_dir, check_file, check_path, copy_file, symlink_file, 
    remove_file, fx_name, file_exists, file_abspath, file_prefix, file_nrows 
)
from hiseq.utils.bam import Bam, Bam2cor, Bam2fingerprint
from hiseq.utils.bed import PeakIDR, BedOverlap, PeakFRiP
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, read_hiseq, list_hiseq_dir, 
    list_hiseq_file, run_shell_cmd, find_longest_common_str, print_dict
)
from hiseq.utils.featurecounts import LibStrand, FeatureCounts


def guess_nsr(fq1, fq2, top_n=10000, cutoff=0.9):
    """
    Guess the input file is NSR library:
    fq1: CT
    fq2: GA
    pct: >90%    
    """
    n1 = 0
    n2 = 0
    try:
        # fq1
        i = 0
        if isinstance(fq1, str):
            for _,seq,_,_ in pyfastx.Fastx(fq1):
                i += 1
                if seq.startswith('CT'):
                    n1 += 1
                if i >= top_n:
                    break
        # fq2
        j = 0
        if isinstance(fq2, str):
            for _,seq,_,_ in pyfastx.Fastx(fq2):
                j += 1
                if seq.startswith('GA'):
                    n2 += 1
                if j >= top_n:
                    break
    except TyepeError as e:
        log.error(e)
    # check
    return n1 / top_n > cutoff and n2 / top_n > cutoff


def hiseq_norm_scale(x, hiseq_type='_r1', by_spikein=False, norm=1000000):
    """
    Parameters
    ---------
    x: str
        The path to CnrRn() dir, hiseq_type=cnr_rn

    cal the norm scale: combine rep_list
    (fragments in bam files: reads, paired-reads)
    bam.count
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('hiseq_norm_scale() skipped, not a hiseq dir: {}'.format(x))
    out = None
    sj = a.align_scale_json
    # direct to r1
    if file_exists(sj):
        d = Config().load(sj)
        out = d
    else:
        bam = a.spikein_bam if by_spikein else a.bam
        m = Bam(bam).count(reads=False) # paired
        try:
            s1 = round(norm / m, 4)
        except ZeroDivisionError as e:
            log.error(e)
            s1 = 1.0
        # save to sj
        d = {
            'smp_name': a.smp_name,
            'is_spikein': by_spikein,
            'norm': norm,
            'map': m,
            'scale': s1
        }
        Config().dump(d, sj)
        out = d
    return out


################################################################################
## Raseq-sub-modules
def rnaseq_trim(x, hiseq_type='rnaseq_r1'):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1
    
    # guess: NSR
    read1: CT
    read2: GA
    
    Trimming: 
    3' adapter: guess
    cut5: 9
    cut3: 9
    min: 20

    operation:
    1. do trimming
    2. create symlink
   """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('rnaseq_trim() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    # do-the-trimming
    fq1, fq2 = a.raw_fq_list
    clean_fq1, clean_fq2 = a.clean_fq_list # output
    # whether to trim or not
    if a.trimmed:
        symlink_file(fq1, clean_fq1)
        symlink_file(fq2, clean_fq2)
    else:
        # guess NSR or not
        is_nsr = guess_nsr(fq1, fq2)
        args_local = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': a.clean_dir,
            'len_min': 20,
            'cut_after_trim': '9,-9' if is_nsr else None,
            'keep_tmp': True,
            'parallel_jobs': 1 # do not allowed > 1 !!!!
        }
        trim = TrimR1(**args_local)
        trim.run()
        ## copy files
        if a.is_paired:
            symlink_file(trim.clean_fq1, clean_fq1)
            symlink_file(trim.clean_fq2, clean_fq2)
        else:
            symlink_file(trim.clean_fq, clean_fq1)
        symlink_file(trim.trim_json, a.trim_json)


def rnaseq_align_spikein(x, hiseq_type='rnaseq_r1'):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1

    Align reads to spikein reference, using STAR
    **Update clean data** 
    non-spikein reads
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('rnaseq_align_spikein() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    if not isinstance(a.spikein_index, str):
        log.info('rnaseq_align_spikein() skipped, no spikein_index')
        return None #
    args_local = a.__dict__
    fq1, fq2 = a.clean_fq_list
    args_init = {
        'aligner': a.aligner,
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.spikein_dir,
        'smp_name': a.smp_name,
        'genome': None,
        'genome_index': None,
        'spikein': None,
        'spikein_index': a.spikein_index,
        'to_rRNA': None,
        'rRNA_index': None,
        'extra_index': None,
        'keep_tmp':a.keep_tmp,
        'unique_only': False,
        'threads': a.threads,
        'parallel_jobs': a.parallel_jobs,
        'overwrite': a.overwrite,
    }
    # args_local.update(args_init)
    args_local = args_init
    if file_exists(a.spikein_bam) and not a.overwrite:
        log.info('rnaseq_align_spikein() skipped, file exists: {}'.format(
            a.spikein_bam))
    else:
        Align(**args_local).run()
    # copy files; go to align_r1 directory
    d_list = list_dir(a.spikein_dir, include_dir=True)
    d_list = [i for i in d_list if os.path.isdir(i)]
    for i in d_list:
        t = read_hiseq(i, 'alignment_rn') # alignment_rn
        if t.is_hiseq:
            symlink_file(t.bam, a.spikein_bam)
            symlink_file(t.align_json, a.spikein_json)
            symlink_file(t.align_flagstat, a.spikein_flagstat)
            symlink_file(t.unmap, a.spikein_unmap)
            symlink_file(t.unmap1, a.spikein_unmap1)
            symlink_file(t.unmap2, a.spikein_unmap2)
            # index bam
            if not file_exists(a.spikein_bam+'.bai'):
                Bam(a.spikein_bam).index()
    # calculate norm scale
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=True)


def rnaseq_align_rRNA(x, hiseq_type='r1'):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1

    Align reads to rRNA, using STAR
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('rnaseq_align_rRNA() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    if not isinstance(a.rRNA_index, str):
        log.info('rnaseq_align_rRNA() skipped, no rRNA_index')
        return None #
    args_local = a.__dict__
    ##################################################################
    # load unmap from spikein
    sp = read_hiseq(a.spikein_dir, 'rx')
    if sp.is_hiseq:
        fq1 = a.spikein_unmap1
        fq2 = a.spikein_unmap2
    else:
        fq1, fq2 = a.clean_fq_list
    # convert str('None') -> None
    if fq2 == 'None':
        fq2 = None
    ##################################################################
    args_init = {
        'aligner': a.aligner,
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.rRNA_dir,
        'smp_name': a.smp_name,
        'genome': a.genome,
        'genome_index': None,
        'spikein': None,
        'spikein_index': None,
        'to_rRNA': False,
        'rRNA_index': None,
        'extra_index': a.rRNA_index,
        'keep_tmp': a.keep_tmp,
        'unique_only': False,
        'threads': a.threads,
        'parallel_jobs': a.parallel_jobs,
        'overwrite': a.overwrite,
    }
    # args_local.update(args_init)
    args_local = args_init
    if file_exists(a.bam) and not a.overwrite:
        log.info('align() skipped, file exists: {}'.format(a.bam))
    else:
        Align(**args_local).run()
    # copy files; go to align_r1 directory
    # align/smp_name/smp_name.bam -> bam_files/smp_name.bam
    d_list = list_dir(a.rRNA_dir, include_dir=True)
    d_list = [i for i in d_list if os.path.isdir(i)]
    for i in d_list:
        t = read_hiseq(i, 'alignment_rn') # alignment_rn
        if t.is_hiseq:
            symlink_file(t.align_json, a.rRNA_json)
            symlink_file(t.align_flagstat, a.rRNA_flagstat)
            symlink_file(t.unmap, a.rRNA_unmap)            
            symlink_file(t.unmap1, a.rRNA_unmap1)
            symlink_file(t.unmap2, a.rRNA_unmap2)


def rnaseq_align_genome(x, hiseq_type='r1'):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1

    Align reads to reference genome, using STAR
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('rnaseq_align_genome() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    args_local = a.__dict__
    ##################################################################
    # load unmap from rRNA/spikein unmap
    r = read_hiseq(a.rRNA_dir, 'rx')
    sp = read_hiseq(a.spikein_dir, 'rx')
    if r.is_hiseq:
        fq1 = a.rRNA_unmap1
        fq2 = a.rRNA_unmap2
    elif sp.is_hiseq:
        fq1 = a.spikein_unmap1
        fq2 = a.spikein_unmap2
    else:
        fq1, fq2 = a.clean_fq_list
    # convert str('None') -> None
    if fq2 == 'None':
        fq2 = None
    ##################################################################
    args_init = {
        'aligner': a.aligner,
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.align_dir,
        'smp_name': a.smp_name,
        'genome': a.genome,
        'genome_index': a.genome_index,
        'spikein': None,
        'spikein_index': None,
        'to_rRNA': False,
        'rRNA_index': None,
        'extra_index': a.extra_index,
        # overwrite args: end # 
        'keep_tmp': a.keep_tmp,
        'unique_only': True,
        'threads': a.threads,
        'parallel_jobs': a.parallel_jobs,
        'overwrite': a.overwrite,
        'verbose': a.verbose,
    }
    # args_local.update(args_init)
    args_local = args_init
    if file_exists(a.bam) and not a.overwrite:
        log.info('align() skipped, file exists: {}'.format(a.bam))
    else:
        Align(**args_local).run()
    # copy files; go to align_r1 directory
    # align/smp_name/smp_name.bam -> bam_files/smp_name.bam
    d_list = list_dir(a.align_dir, include_dir=True)
    d_list = [i for i in d_list if os.path.isdir(i)]
    for i in d_list:
        t = read_hiseq(i, 'alignment_rn') # alignment_rn
        if t.is_hiseq:
            symlink_file(t.bam, a.bam)
            symlink_file(t.align_json, a.align_json)
            symlink_file(t.align_flagstat, a.align_flagstat)
            symlink_file(t.unmap, a.unmap)
            symlink_file(t.unmap1, a.unmap1)
            symlink_file(t.unmap2, a.unmap2)
            # index bam
            if not file_exists(a.bam+'.bai'):
                Bam(a.bam).index()
    # calculate norm scale
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=False)
        

def hiseq_merge_bam(x, hiseq_type='_rn'):
    """
    Merge bam files from r1, re-calculate the norm_scale
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('cnr_merge_bam() skipped, not a cnr_rn dir: {}'.format(x))
    bam_list = list_hiseq_file(x, 'bam', 'r1')
    cmd = ' '.join([
        'samtools merge -',
        ' '.join(bam_list),
        '| samtools sort -o {} -'.format(a.bam),
        '&& samtools index {}'.format(a.bam)])
    if not all(file_exists(bam_list)):
        raise ValueError('bam file not exists: {}'.format(bam_list))
    if file_exists(a.bam) and not a.overwrite:
        log.info('merge_bam() skipped, file exists: {}'.format(a.bam))
    else:
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('merge_bam() failed.')
    # check-point
    if not file_exists(a.bam):
        raise ValueError('cnr_merge_bam() failed, see: {}'.format(a.bam_dir))
    # save align_json
    d = {'name': a.smp_name}
    for aj in list_hiseq_file(x, 'align_json', 'r1'):
        da = Config().load(aj)
        d.update({
            'index': da.get('index', None),
            'unique_only': da.get('unique_only', False),
            'total': d.get('total', 0) + da.get('total', 0),
            'map': d.get('map', 0) + da.get('map', 0),
            'unique': d.get('unique', 0) + da.get('unique', 0),
            'multi': d.get('multi', 0) + da.get('multi', 0)
        })
    Config().dump(d, a.align_json)
    # calculate norm scale
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=False) # to genome
    

        
def rnaseq_quant(x, hiseq_type='r1'):
    """
    Parameters
    ---------
    x: str
        The path to RnaseqR1() dir, hiseq_type=rnaseq_r1

    Run FeatureCounts for the bam file
    """    
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('rnaseq_quant() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    # check args
    gtf = a.gene_gtf
    bam = a.bam
    threads = getattr(a, 'threads', 1)
    overwrite = getattr(a, 'overwrite', False)
    if not all(file_exists([gtf, bam])): 
        msg = 'rnaseq_quant() skipped, file not exists, bam: {}, gtf: {}'.format(
            file_exists(bam), file_exists(gtf)
        )
        log.error(msg)
        return None
    # args_local = a.__dict__
    args_local = {
        'gtf': gtf,
        'bam_list': bam,
        'outdir': a.count_dir,
        'threads': a.threads,
        'overwrite': a.overwrite
    }
    # guess sense strand
    s = LibStrand(bam=bam, gtf=gtf, clean_up=True) # sense_strand, anti_strand
    ##################################################################
    # sense strand
    args_sense = {
        'prefix': os.path.basename(a.count_sens),
        'strandness': s.sense_strand
    }
    args_sense.update(args_local)
    fc_sense = FeatureCounts(**args_sense).run()
    ##################################################################
    # anti strand
    args_anti = {
        'prefix': os.path.basename(a.count_anti),
        'strandness': s.anti_strand
    }
    args_anti.update(args_local)
    fc_anti = FeatureCounts(**args_anti).run()
    ##################################################################
    # save strandness to file
    d = {
        'project_dir': x,
        'bam': bam,
        'gtf': gtf,
        'sense': s.sense_strand,
        'anti': s.anti_strand,
    }
    Config().dump(d, a.strandness_json)


def hiseq_bam2bw(x, hiseq_type='_r1'):
    """
    cnr_bam_to_bw()
    Check the scale only for : TE/Genome
    normalizeUsing: RPGC, CPM
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('hiseq_bam_to_bw() skipped, not a hiseq dir: {}'.format(x))
        return None
    args = {
        'bam': a.bam,
        'prefix': a.smp_name,
        'outdir': a.bw_dir,
        'binsize': a.binsize, # default: 50
        'strandness': 12, # both: fwd, rev
        'genome': a.genome,
        'scaleFactor': 1.0, # force
        'normalizeUsing': 'RPGC',
        'overwrite': a.overwrite,
        'genome_size': a.genome_size
    }
    Bam2bw(**args).run()

    
def rnaseq_deseq(x, hiseq_type='rx'):
    """
    DEseq analysis
    Using DESeq2, edgeR, ...
    """    
    a = read_hiseq(x, 'rnaseq_rx')
    if not a.is_hiseq:
        log.error('rnaseq_deseq() skipped, not a rnaseq_rx dir: {}'.format(x))
        return None
    # prepare R commands
    pkg_dir = os.path.dirname(hiseq.__file__)
    deseq_r = os.path.join(pkg_dir, 'bin', 'run_rnaseq.R')
    stdout = os.path.join(a.deseq_dir, 'deseq.stdout')
    stderr = os.path.join(a.deseq_dir, 'deseq.stderr')
    cmd_shell = os.path.join(a.deseq_dir, 'cmd.sh')
    run_go = 1 if a.genome else 0
    cmd = ' '.join([
        'Rscript {}'.format(deseq_r),
        '{} {}'.format(x, run_go),
        '1> {} 2> {}'.format(stdout, stderr),
    ])
    with open(cmd_shell, 'wt') as w:
        w.write(cmd + '\n')
    if file_exists(a.deseq_fix_xls) and not a.overwrite:
        log.info('rnaseq_deseq() skipped, file exists: {}'.format(
            a.deseq_fix_xls))
    else:
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('DESeq2 failed.')

    
def deseq_salmon(x, hiseq_type='rx'):
    """
    DEseq analysis, for salmon output
    Using DESeq2, edgeR, ...
    """    
    a = read_hiseq(x, 'rnaseq_rx')
    if not a.is_hiseq:
        log.error('rnaseq_deseq() skipped, not a rnaseq_rx dir: {}'.format(x))
        return None
    # prepare R commands
    pkg_dir = os.path.dirname(hiseq.__file__)
    deseq_r = os.path.join(pkg_dir, 'bin', 'run_deseq_salmon.R')
    stdout = os.path.join(a.deseq_dir, 'deseq.stdout')
    stderr = os.path.join(a.deseq_dir, 'deseq.stderr')
    cmd_shell = os.path.join(a.deseq_dir, 'cmd.sh')
    run_go = 1 if a.genome else 0
    cmd = ' '.join([
        'Rscript {}'.format(deseq_r),
        '{} {}'.format(x, run_go),
        '# 1> {} 2> {}'.format(stdout, stderr),
    ])
    with open(cmd_shell, 'wt') as w:
        w.write(cmd + '\n')
    if file_exists(a.deseq_fix_xls) and not a.overwrite:
        log.info('deseq_salmon() skipped, file exists: {}'.format(
            a.deseq_fix_xls))
    else:
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('DESeq2 failed.')
    

################################################################################
# quality control
# 1. trim
# 2. align
# 3. satuation
# 4. fragment
# 5. ...
def qc_trim_summary(x, hiseq_type='r1'):
    """
    # format:
    # name, input, output, out_pct, rm_pct
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_trim_summary() failed, not a hiseq_r1: {}'.format(x))
        return None
    # option-1: stat.yaml
    # option-2: stat.txt
    stat_json = getattr(a, 'trim_json', None)
    stat_txt = getattr(a, 'trim_stat', None)
    # format:
    # name, total, too_short, dup, too_short2, clean, percent
    d = {
        'name': a.smp_name,
        'input': 1,
        'output': 1,
        'out_pct': 100.0,
        'rm_pct': 0,
    }
    if file_exists(stat_json):
        df = Config().load(stat_json) # laod data
        d['input'] = int(df.get('total', 1))
        d['output'] = int(df.get('clean', 1))
        d['out_pct'] = float(df.get('percent', 100.0))
        d['rm_pct'] = 100.0 - d['out_pct']
    elif file_exists(stat_txt):
        try:
            s = None # init
            with open(stat_txt) as r:
                for line in r:
                    if line.startswith('#'):
                        continue
                    s = line.strip().split('\t')
                    break
            if isinstance(s, list):
                d = {
                    'name': s[0],
                    'input': int(s[1]),
                    'output': int(s[-2]),
                    'out_pct': float(s[-1]),
                    'rm_pct': 100.0 - float(s[-1]),
                }
        except IOError as e:
            log.error(e)
    else:
        log.error('trim.stat not exists: {}'.format(stat_txt))
    # update pct
    d['out_pct'] = float('{:.2f}'.format(d['out_pct']))
    d['rm_pct'] = float('{:.2f}'.format(d['rm_pct']))
    # save to new file
    Config().dump(d, a.trim_summary_json)
    return d


def qc_align_summary(x, hiseq_type='r1'):
    """
    Organize the alignment:
    
    output:
    name, total, map, unique, multi, spikein, rRNA, unmap
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_align_summary() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    # spikein, rRNA, genome
    n_total = 0
    if file_exists(a.spikein_json):
        sp = Config().load(a.spikein_json)
        n_sp = sp.get('map', 0) if isinstance(sp, dict) else 0
        n_total = sp.get('total', 0) if isinstance(sp, dict) else 0
    else:
        n_sp = 0
    # rRNA
    if file_exists(a.rRNA_json):
        r = Config().load(a.rRNA_json)
        n_rRNA = r.get('map', 0) if isinstance(r, dict) else 0
        n_total = r.get('total', 0) if isinstance(r, dict) else 0
    else:
        n_rRNA = 0
    # name, total, map, unique, multi, unmap
    if file_exists(a.align_json):
        df = Config().load(a.align_json)
        if n_total == 0:
            n_total = df.get('total', 1)
        df.update({
            'total': n_total,
            'spikein': n_sp,
            'rRNA': n_rRNA,
        })
        Config().dump(df, a.align_summary_json)
    else:
        log.warning('qc_align_summary() skipped, no align_json')
        return None

    
def qc_bam_cor(x, hiseq_type='rn', bam_type='r1'):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq
        
    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']
        
    bam_type:  str
        The hiseq type of bam file, options: ['r1', 'rn', 'rx']
        default: ['r1']

    Compute correlation (pearson) between replicates
    window = 500bp

    eg:
    multiBamSummary bins --binSize 500 --smartLabels -o *bam.npz \
        --outRawCounts *counts.tab -b bam
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_bam_cor() failed, not a hiseq dir: {}'.format(x))
        return None
    # r1
    bam_list = None
    if a.is_hiseq_r1: 
        pass # require > 1 bam
    elif a.is_hiseq_rn:
        bam_list = list_hiseq_file(x, 'bam', bam_type)
    elif a.is_hiseq_rx:
        if a.hiseq_type in ['chipseq_rx', 'cnr_rx', 'cnt_rx']:
            bam_list = [a.ip_bam, a.input_bam]
        elif a.hiseq_type in ['rnaseq_rx']:
            bam_list = [a.mut_bam, a.wt_bam]
        else:
            pass
    else:
        log.error('qc_bam_cor() failed, unknown hiseq dir: {}'.format(x))
    # check point
    if bam_list is None:
        log.error('qc_bam_cor() failed, require >= 2 bam files: {}'.format(x))
        return None
    # bam_list = list_hiseq_file(x, 'bam_rmdup', bam_type)
    args = {
        'bam_list': bam_list,
        'outdir': a.qc_dir,
        'prefix': '06.bam_cor',
        'threads': a.threads,
        'overwrite': a.overwrite,
        'binsize': 500, # a.binsize,
    }
    if file_exists(a.bam_cor_heatmap_png) and not a.overwrite:
        log.info('qc_bam_cor() skipped, file exists')
    else:
        print('\n'.join(bam_list))
        if all(file_exists(bam_list)):
            Bam2cor(**args).run()
        else:
            log.error('qc_bam_cor() failed, bam files not exists')

            
################################################################################
## function for tss,genebody enrich: for general hiseq purpose
##
## computematrix, plotProfile from deeptools
def qc_tss_enrich_tool(x, hiseq_type='r1', bw_type='r1',
                       subcmd='reference-point', **kwargs):
    """
    1. get bw list
    2. get labels
    3. get title
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_tss_enrich() failed, not a hiseq dir: {}'.format(x))
        return None
    bed = getattr(a, 'gene_bed', None)
    if not file_exists(bed):
        log.error('qc_tss_enrich_tool() failed, bed not exists: {}'.format(
            bed))
        return None
    arg_bed = '-R {}'.format(bed)
    # default values
    # -b 2000 -a 2000 --binSize 500
    upstream = kwargs.get('upstream', 5000)
    downstream = kwargs.get('downstream', 5000)
    regionbody = kwargs.get('regionbody', 5000)
    binsize = kwargs.get('binSize', 50)
    # -b 2000 -a 2000 --binSize 50
    arg_body = '-b {} -a {} --binSize {}'.format(
        upstream, downstream, binsize)
    # add genebody for scale-regions
    # subcmd: scale-regions, reference-point
    if subcmd == 'scale-regions':
        arg_body += ' -m {}'.format(regionbody)
    # r1
    bw = ''
    arg_bw = ''
    arg_label = ''
    arg_title = ''
    per_group = '--perGroup'
    # r1: atac/rnaseq/cnr...
    if a.is_hiseq_r1:
        arg_title = '--plotTitle {}'.format(a.smp_name)
        if a.hiseq_type.startswith('rnaseq_'):
            bw = ' '.join([a.bw_fwd, a.bw_rev])
            arg_label = '--samplesLabel {} {}'.format('fwd', 'rev')
        else:
            bw = a.bw
            arg_label = '--samplesLabel {}'.format(a.smp_name)
        arg_bw = '-S {}'.format(bw)
    elif a.is_hiseq_rn:
        arg_title = '--plotTitle {}'.format(a.smp_name)
        if a.hiseq_type.startswith('rnaseq'):
            b1 = list_hiseq_file(x, 'bw_fwd', 'r1')
            b2 = list_hiseq_file(x, 'bw_rev', 'r1')
            n1 = list_hiseq_file(x, 'smp_name', 'rn')
            n2 = list_hiseq_file(x, 'smp_name', 'rn')
            n1 = [i.replace(a.smp_name+'_', '') for i in n1]
            n2 = [i.replace(a.smp_name+'_', '') for i in n2]
            n_list = [i+k for i in n1 + n2 for k in ['_fwd', '_rev']]
            bw_list = [a.bw_fwd, a.bw_rev] + b1 + b2
            n_list = ['merge_fwd', 'merge_rev'] + n_list
            bw = ' '.join(bw_list)
            arg_label = '--samplesLabel {}'.format(' '.join(n_list))
            arg_title = '--plotTitle {}'.format(a.smp_name)
        else: # atac, cnr, chip, ...
            bw_list = list_hiseq_file(x, 'bw', 'r1')
            smp_name = list_hiseq_file(x, 'smp_name', 'r1') # multi
            smp_name = [i.replace(a.smp_name+'_', '') for i in smp_name]
            # add merge
            bw_list.insert(0, a.bw)
            smp_name.insert(0, a.smp_name)
            # prepare args
            bw = ' '.join(bw_list)
            arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        arg_bw = '-S {}'.format(bw)
    elif a.is_hiseq_rx:
        arg_title = '--plotTitle {}'.format(a.smp_name)
        if a.hiseq_type in ['rnaseq_rx']:
            # mut, wt
            bw_list = [a.mut_bw_fwd, a.mut_bw_rev, a.wt_bw_fwd, a.wt_bw_rev]
            arg_bw = ' '.join(bw_list)
            # shorter name
            s = find_longest_common_str(a.mut_name, a.wt_name)
            s1 = a.mut_name.replace(s, '')
            s2 = a.wt_name.replace(s, '')
            ss = '{}.vs.{}'.format(s1, s2)
            smp_name = [s1+'.fwd', s1+'.rev', s2+'.fwd', s2+'.rev']
            arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        elif a.hiseq_type in ['chipseq_rx', 'cnr_rx', 'cnt_rx']:
            # ip, input
            bw_list = [a.ip_bw, a.input_bw, a.bw]
            bw = ' '.join(bw_list)
            # shorter name
            s = find_longest_common_str(a.ip_name, a.input_name)
            s1 = a.ip_name.replace(s, '')
            s2 = a.input_name.replace(s, '')
            ss = '{}.vs.{}'.format(s1, s2)
            smp_name = [s1, s2, ss]
            arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        else:
            pass
        arg_bw = '-S {}'.format(bw)
    else:
        log.error('qc_genebody_enrich() failed, unknown hiseq dir: {}'.format(x))
    if bw == '':
        return None
    arg_ma = ' '.join([arg_label, arg_body, arg_bw, arg_bed])
    arg_plot = ' '.join([arg_title, per_group])
    # return arguments
    return (arg_ma, arg_plot)            
            

def qc_genebody_enrich(x, hiseq_type='r1', bw_type='r1', **kwargs):
    """
    Calculate the Genebody enrichment

    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['r1']

    bw_type:  str
        The hiseq type of bigWig file, options: ['r1', 'rn', 'rx']
        default: ['r1']

    $ computeMatrix scale-regions -R gene.bed -S f.bigWig -o mat.gz
    $ plotProfile -m mat.gz -o gene_body.png
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    # arg_bw, arg_label, arg_title, per_group
    args = qc_tss_enrich_tool(x, hiseq_type, bw_type,
        subcmd='scale-regions', **kwargs)
    if args is None:
        log.error('qc_genebody_enrich() failed')
        return None
    cmd = ' '.join([
        '{}'.format(shutil.which('computeMatrix')),
        'scale-regions',
        '-o {}'.format(a.genebody_enrich_matrix),
        '--sortRegions descend --skipZeros',
        args[0], # -R -S -b -a -m --binSize
        '--smartLabels',
        '-p {}'.format(a.threads),
        '2> {}'.format(a.genebody_enrich_matrix_log),
        '&& {}'.format(shutil.which('plotProfile')),
        '-m {}'.format(a.genebody_enrich_matrix),
        '-o {}'.format(a.genebody_enrich_png),
        args[1], # --plotTitle --perGroup
        '--dpi 300'
    ])
    if file_exists(a.genebody_enrich_png) and not a.overwrite:
        log.info('qc_genebody_enrich() skipped, file exists: {}'.format(
            a.genebody_enrich_png))
    else:
        if not file_exists(getattr(a, 'gene_bed', None)):
            log.error('qc_tss() skipped, gene_bed not found')
        else:
            with open(a.genebody_enrich_cmd, 'wt') as w:
                w.write(cmd + '\n')
            try:
                run_shell_cmd(cmd)
            except:
                log.error('qc_genebody_enrich() failed, see: {}'.format(
                    a.genebody_enrich_matrix_log))

