#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
General modules for Chipseqseq analysis

analysis-module:
"""

import os
import re
import glob
import shutil
import pysam
import hiseq
from hiseq.trim.trimmer import TrimR1
from hiseq.align.align import Align
from hiseq.bam2bw.bam2bw import Bam2bw, bw_compare
from hiseq.chipseq.callpeak import CallPeak
from hiseq.fragsize.fragsize import BamFragSize, BamFragSizeR1
from hiseq.utils.bam import Bam, Bam2cor, Bam2fingerprint
from hiseq.utils.bed import PeakIDR, BedOverlap, PeakFRiP
from hiseq.utils.file import (
    list_file, list_dir, check_file, check_path, copy_file, symlink_file, 
    remove_file, fx_name, file_exists, file_abspath, file_prefix, file_nrows
)
from hiseq.utils.utils import (
    log, update_obj, Config, get_date, read_hiseq, list_hiseq_file, run_shell_cmd, 
    find_longest_common_str
)


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
## main ##
#
# chipseq_trim 
# chipseq_align_spikein
# chipseq_align_genome
#
def chipseq_trim(x, hiseq_type='r1'):
    """
    Parameters
    ---------
    x: str
        The path to ChipseqR1() dir, hiseq_type=chipseq_r1

    option-1: cut reads from 3' end, to X-nt in length (default: 50)

    operation:
    1. do trimming
    2. create symlink

    path -> (ChipseqReader) -> args
    fq1, fq2
    clean_fq1, clean_fq2, clean_dir
    cut_to_length, recursive
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('chipseq_trim() skipped, not a chipseq_r1 dir: {}'.format(x))
        return None
    # do-the-trimming
    fq1, fq2 = a.raw_fq_list
    clean_fq1, clean_fq2 = a.clean_fq_list
    # whether to trim or not
    if a.trimmed:
        symlink_file(fq1, clean_fq1)
        symlink_file(fq2, clean_fq2)
    else:
        args_local = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': a.clean_dir,
            'smp_name': a.smp_name,
            'library_type': None,
            'len_min': 20,
            'cut_to_length': a.cut_to_length,
            'recursive': a.recursive,
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

        
def chipseq_align_spikein(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to ChipseqR1() dir, hiseq_type=chipseq_r1

    Align reads to reference genome, using bowtie2, SE|PE

    --sensitive --local -X 2000  # --devotail

    samtools view -bhS -f 2 -F 1804

    exclude: -F
        4    0x4    Read unmapped,
        8    0x8    Mate unmapped,
        256  0x100  Not primary alignment,
        512  0x200  Reads fails platform quality checks,
        1024 0x400  Read is PCR or optical duplicate
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('chipseq_align_spikein() skipped, not a chipseq_r1 \
            dir: {}'.format(x))
        return None
#     if not isinstance(a.spikein_index, str):
#         log.info('chipseq_align_spikein() skipped, no spikein_index')
#         return None #
    args_local = a.__dict__
    fq1, fq2 = getattr(a, 'clean_fq_list', [None, None])
    args_init = {
        'aligner': a.aligner,
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.spikein_dir,
        'unique_only': False,
        'spikein': a.spikein,
        'spikein_index': a.spikein_index,
        'smp_name': a.smp_name,
        'genome': None,
        'genome_index': None,
        'to_rRNA': None,
        'rRNA_index': None,
        'extra_index': None,
        'keep_tmp': a.keep_tmp,
        'max_fragment': 2000,
        'threads': a.threads,
        'parallel_jobs': a.parallel_jobs,
        'overwrite': a.overwrite,
        'verbose': False,
        'extra_para': align_extra, # specific for chipseq
    }
    args_local = args_init
    if file_exists(a.spikein_bam) and not a.overwrite:
        log.info('chipseq_align_spikein() skipped, file exists: {}'.format(
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
            symlink_file(t.unmap1, a.unmap1)
            symlink_file(t.unmap2, a.unmap2)
    # calculate norm scale
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=True)


def chipseq_align_genome(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to ChipseqR1() dir, hiseq_type=chipseq_r1

    Align reads to reference genome, using bowtie2
    --sensitive --local -X 2000  # --devotail
    samtools view -bhS -f 2 -F 1804
    ######
    chipseq: "-I 10 -X 700"
    exclude: -F
        4    0x4    Read unmapped,
        8    0x8    Mate unmapped,
        256  0x100  Not primary alignment,
        512  0x200  Reads fails platform quality checks,
        1024 0x400  Read is PCR or optical duplicate
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('chipseq_align_genome() skipped, not a chipseq_r1 dir: \
            {}'.format(x))
        return None
    args_local = a.__dict__
    fq1, fq2 = getattr(a, 'clean_fq_list', [None, None])
    align_extra = '-I 10 -X 700' if a.is_paired else ''
    align_extra = '-I 10 -X 2000' if a.hiseq_type.startswith('atac') else \
        '-I 10 -X 700' if a.hiseq_type.startswith('cn') else ''  # chipseq, cnt
    args_init = {
        'aligner': a.aligner,
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.align_dir,
        'smp_name': a.smp_name,
        'genome': a.genome,
        'extra_index': a.extra_index,
        'keep_tmp': a.keep_tmp,
        'max_fragment': 2000,
        'unique_only': True,
        'threads': a.threads,
        'parallel_jobs': a.parallel_jobs,
        'overwrite': a.overwrite,
        'verbose': False,
        'extra_para': align_extra, # specific for chipseq
    }
    args_local = args_init
    if file_exists(a.bam) and not a.overwrite:
        log.info('chipseq_align_genome() skipped, file exists: {}'.format(a.bam))
    else:
        print('!A1', args_local)
        Align(**args_local).run()
    # copy files; go to align_r1 directory
    # align/smp_name/smp_name.bam -> bam_files/smp_name.bam
    d_list = list_dir(a.align_dir, include_dir=True)
    d_list = [i for i in d_list if os.path.isdir(i)]
    for i in d_list:
        t = read_hiseq(i, 'alignment_rn') # alignment_rn
        if t.is_hiseq:
            symlink_file(t.bam, a.bam_raw)
            symlink_file(t.align_stat, a.align_stat)
            symlink_file(t.align_json, a.align_json)
            symlink_file(t.align_flagstat, a.align_flagstat)
    # remove PCRdup
    if a.rmdup:
        if file_exists(a.bam_raw):
            if check_file(a.bam, check_empty=True) and not a.overwrite:
                log.info('rmdup() skipped, file exists: {}'.format(
                    a.bam))
            else:
                Bam(a.bam_raw).rmdup(outfile=a.bam)
                # Bam(a.bam).index()
    else:
        symlink_file(a.bam_raw, a.bam)
    s = hiseq_norm_scale(x, hiseq_type=hiseq_type, by_spikein=False)


def chipseq_merge_bam(x, hiseq_type='_rn'):
    """
    Merge bam files from r1, re-calculate the norm_scale
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('chipseq_merge_bam() skipped, not a chipseq_rn dir:\
            {}'.format(x))
    bam_list = list_hiseq_file(x, 'bam', 'r1')
    cmd = ' '.join([
        'samtools merge -',
        ' '.join(bam_list),
        '| samtools sort -o {} -'.format(a.bam_raw),
        '&& samtools index {}'.format(a.bam_raw)])
    if not all(file_exists(bam_list)):
        raise ValueError('bam file not exists: {}'.format(bam_list))
    if file_exists(a.bam_raw) and not a.overwrite:
        log.info('merge_bam() skipped, file exists: {}'.format(a.bam_raw))
    else:
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('merge_bam() failed.')
    # rmdup
    if a.rmdup:
        if file_exists(a.bam_raw):
            if check_file(a.bam, check_empty=True) and not a.overwrite:
                log.info('rmdup() skipped, file exists: {}'.format(
                    a.bam))
            else:
                Bam(a.bam_raw).rmdup(a.bam)
    else:
        symlink_file(a.bam_raw, a.bam)
    # check-point
    if not file_exists(a.bam):
        raise ValueError('chipseq_merge_bam() failed, see: {}'.format(a.bam_dir))
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


def hiseq_pcr_dup(x, hiseq_type='r1'):
    """
    Organize the PCR dup rate:
    total: align_json {'map'}
    nodup: align_scale_json {'map'}
    name total nodup dup
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_dup_summary() skipped, not a rnaseq_r1 dir: {}'.format(x))
        return None
    # name total nodup dup
    out = None
    j1 = list_hiseq_file(x, 'align_json', hiseq_type)
    j2 = list_hiseq_file(x, 'align_scale_json', hiseq_type)
    if isinstance(j1, list) and isinstance(j2, list):
        d1 = Config().load(j1[0])
        d2 = Config().load(j2[0])
        out = {
            'name': d1.get('name'),
            'total': d1.get('map', 0),
            'nodup': d2.get('map', 0),
            'dup': d1.get('map', 0) - d2.get('map', 0)
        }
        out['dup_pct'] = round(out['dup']/out['total'], 4)
        # Config().dump(out, a.dup_summary_json)
        Config().dump(out, a.pcr_dup_json)
    else:
        log.warning('qc_dup_summary() skipped')
    return out
    
    
def chipseq_call_peak(x, hiseq_type='r1'):
    """
    Call peaks using MACS2 and SEACR
    -f BAMPE
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('chipseq_call_peak() failed, not a hiseq dir: {}'.format(x))
        return None
    # r1
    ip_bam = None
    input_bam = None
    if a.is_hiseq_r1 or a.is_hiseq_rn:
        if hasattr(a, 'bam_rmdup'):
            ip_bam = a.bam_rmdup
        elif hasattr(a, 'bam'):
            ip_bam = a.bam
        else:
            pass
    elif a.is_hiseq_rx:
        if a.hiseq_type in ['chipseq_rx', 'chipseq_rx', 'cnt_rx']:
            # ip, input
            ip_bam = a.ip_bam
            input_bam = a.input_bam
    else:
        log.error('chipseq_call_peak() failed, unknown hiseq dir: {}'.format(x))
    # check point
    if ip_bam is None:
        return None
    # cmd
    args = {
        'ip': ip_bam, # rm_dup
        'input': input_bam,
        'outdir': a.peak_dir,
        'prefix': a.smp_name,
        'genome': a.genome,
        'genome_size': a.genome_size,
        'genome_size_file': a.genome_size_file
    }
    CallPeak(method='macs2', **args).run()


def chipseq_bw_compare(x, hiseq_type='rx'):
    """
    Create bw, ip over input
    bwCompare()
    """
    a = read_hiseq(x, hiseq_type)
    if not (a.is_hiseq and hiseq_type.endswith('rx')):
        log.error('chipseq_bw_compare() failed, only support for rx')
        return None
    bw_compare(a.ip_bw, a.input_bw, a.bw, 'log2',
        threads=a.threads, binsize=10)
    
    
def get_mito_count(x):
    """
    Count reads on chrM (MT), ...
    """
    try:
        lines = pysam.idxstats(x).splitlines()
        lines = [i for i in lines if re.search('^(chrM|MT)', i)]
        n_mt = sum(
            map(int, [x.split("\t")[2]
                      for x in lines if not x.startswith("#")]))
    except IndexError as e:
        msg = "can't get number of reads from bamfile, msg=%s, data=%s".format(
            (e, lines))
        log.error(msg)
        n_mt = 0
    return n_mt


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
        'strandness': 0, # non-strandness
        'genome': a.genome,
        'scaleFactor': 1.0, # force
        'overwrite': a.overwrite,
        'genome_size': a.genome_size
    }
    Bam2bw(**args).run()


################################################################################
## Quality control matrix for chipseqseq analysis
#
# 1. raw_data: > 20M
# 2. clean_data: > 20M (clean_pct, too-short)
# 3. align: map 20M (chrM, map%)
# 4. peaks
# 5. FRiP
# 6. peak_overlap (multi: _rn)
# 7. frag_size
# 8. tss_enrich
# 9. genebody_enrich
# 10. bam_fingerprint
# 11. bam_cor
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
    # pcr_dup
    if file_exists(a.pcr_dup_json):
        sd = Config().load(a.pcr_dup_json)
        n_dup = sd.get('dup', 0)
        n_nodup  = sd.get('nodup', 0)
    else:
        n_dup = n_nodup = 0
    # name, total, map, unique, multi, unmap
    if file_exists(a.align_json):
        df = Config().load(a.align_json)
        if n_total == 0:
            n_total = df.get('total', 1)
        df.update({
            'total': n_total,
            'spikein': n_sp,
            'chrM': get_mito_count(a.bam),
            'dup': n_dup,
            'nodup': n_nodup,
        })
        Config().dump(df, a.align_summary_json)
    else:
        log.warning('qc_align_summary() skipped, no align_json')
        return None


def qc_lendist(x, hiseq_type='r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_lendist() failed, not a hiseq dir: {}'.format(x))
        return None
    if os.path.exists(a.lendist_csv) and a.overwrite is False:
        log.info('qc_lendist() skipped: file exists: {}'.format(
            a.lendist_csv))
    else:
        b = BamFragSizeR1(
            bam=a.bam,
            outdir=a.qc_dir,
            labels=None,
            csv_file=a.lendist_csv,
        )
        try:
            b.run()
        except:
            log.error('qc_lendist() failed, {}'.format(a.lendist_csv))


def qc_frip(x, hiseq_type='r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_frip() failed, not a hiseq dir: {}'.format(x))
        return None
    qc_frip_dir = os.path.join(a.qc_dir, 'frip_files')
    if check_file(a.frip_json, check_empty=True):
        log.info('qc_frip() skipped, file exists: {}'.format(
            a.frip_json))
    else:
        try:
            # dict: index, total, n, frip
            s = PeakFRiP(peak=a.peak, bam=a.bam,
                method='featureCounts', # bedtools
                outdir=qc_frip_dir).run()
            # number of peaks
            s.update({
                'name': a.project_name,
                'n_peaks': file_nrows(a.peak),
            })
            Config().dump(s, a.frip_json)
        except:
            log.error('PeakFRiP() failed, see: {}'.format(a.frip_json))


################################################################################
## function for tss,genebody enrich: for general hiseq purpose
##
## computematrix, plotProfile from deeptools
##
def qc_tss_enrich_tool(x, hiseq_type='r1', bw_type='r1',
                       subcmd='reference-point', **kwargs):
    """
    1. get bw list
    2. get labels
    3. get title
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_tss_enrich_tool() failed, not a hiseq dir: {}'.format(x))
        return None
    bed = getattr(a, 'gene_bed', None)
    if not file_exists(bed):
        log.error('qc_tss_enrich_tool() failed, bed not exists: {}'.format(
            bed))
        return None
    arg_bed = '-R {}'.format(bed)
    # default values
    # -b 2000 -a 2000 --binSize 10
    upstream = kwargs.get('upstream', 1000)
    downstream = kwargs.get('downstream', 1000)
    regionbody = kwargs.get('regionbody', 2000)
    binsize = kwargs.get('binSize', 10)
    # -b 2000 -a 2000 --binSize 10
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
    # r1: atac/rnaseq/chipseq...
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
        else: # atac, chipseq, chip, ...
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
            smp_name = [s1+'.fwd', s1+'.rev', s2+'fwd', s2+'.rev']
            arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        elif a.hiseq_type in ['chipseq_rx', 'chipseq_rx', 'cnt_rx']:
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


def qc_tss_enrich(x, hiseq_type='r1', bw_type='r1', **kwargs):
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

    $ computeMatrix referencepoint -b -R gene.bed -S in.bw -o mat.gz
    $ plotProfile -m mat.gz -o tss.png
    """
    a = read_hiseq(x, hiseq_type) # for general usage    
    if not a.is_hiseq:
        log.error('qc_tss_enrich() failed, not a hiseq dir: {}'.format(x))
        return None
    # arg_bw, arg_label, arg_title, per_group
    args = qc_tss_enrich_tool(x, hiseq_type, bw_type,
        subcmd='reference-point', **kwargs)
    if args is None:
        log.error('qc_tss_enrich() failed')
        return None
    cmd = ' '.join([
        '{}'.format(shutil.which('computeMatrix')),
        'reference-point',
        '--referencePoint TSS',
        '--sortRegions descend --skipZeros',
        args[0], # -R -S -b -a --binSize
        '-o {}'.format(a.tss_enrich_matrix),
        '-p {}'.format(a.threads),
        '2> {}'.format(a.tss_enrich_matrix_log),
        '&& {}'.format(shutil.which('plotProfile')),
        '-m {}'.format(a.tss_enrich_matrix),
        '-o {}'.format(a.tss_enrich_png),
        args[1], # --plotTitle --perGroup
        '--dpi 300',
        ])
    if file_exists(a.tss_enrich_png) and not a.overwrite:
        log.info('qc_tss() skipped, file exists: {}'.format(a.tss_enrich_png))
    else:
        if not file_exists(getattr(a, 'gene_bed', None)):
            log.error('qc_tss() skipped, gene_bed not found')
        else:
            with open(a.tss_enrich_cmd, 'wt') as w:
                w.write(cmd + '\n')
            try:
                run_shell_cmd(cmd)
            except:
                log.error('failed to run computeMatrix, see: {}'.format(
                    a.tss_enrich_matrix_log))


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
    if not a.is_hiseq:
        log.error('qc_genebody_enrich() failed, not a hiseq dir: {}'.format(x))
        return None
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
################################################################################


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
    bam_list = list_hiseq_file(x, 'bam', bam_type)
    args = {
        'bam_list': bam_list,
        'outdir': a.qc_dir,
        'prefix': '06.bam_cor',
        'threads': a.threads,
        'overwrite': a.overwrite,
        'binsize': a.binsize,
    }
    if file_exists(a.bam_cor_heatmap_png) and not a.overwrite:
        log.info('qc_bam_cor() skipped, file exists')
    else:
        if all(file_exists(bam_list)):
            Bam2cor(**args).run()
        else:
            log.error('qc_bam_cor() failed, bam files not exists')


def qc_peak_idr(x, hiseq_type='rn', peak_type='r1'):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    peak_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_peak_idr() failed, not a hiseq dir: {}'.format(x))
        return None
    peak_list = list_hiseq_file(x, 'peak', peak_type)
    args = {
        'peak_list': peak_list,
        'outdir': a.qc_dir,
        'prefix': '07.peak_idr',
        'overwrite': a.overwrite
    }
    if file_exists(a.peak_idr_png) and not a.overwrite:
        log.info('qc_peak_idr() skipped, file exists')
    else:
        if all(file_exists(peak_list)):
            PeakIDR(**args).run()
        else:
            log.error('qc_peak_idr() failed, peak files not exists')


def qc_peak_overlap(x, hiseq_type='rn', peak_type='r1'):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    peak_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_peak_overlap() failed, not a hiseq dir: {}'.format(x))
        return None
    peak_list = list_hiseq_file(x, 'peak', peak_type)
    args = {
        'peak_list': peak_list,
        'outdir': a.qc_dir,
        'prefix': '08.peak_overlap',
        'overwrite': a.overwrite
    }
    if file_exists(a.peak_overlap_png) and not a.overwrite:
        log.info('qc_peak_overwrite() skipped, file exists')
    else:
        if all(file_exists(peak_list)):
            # takes too-long for peaks > 20k
            n_peaks = [file_nrows(i) > 20000 for i in peak_list]
            if all([i > 20000 for i in n_peaks]):
                log.warning('qc_peak_overlap() skipped, too many peaks (>20000): [{}]'.format(
                    '\,'.join(list(map(str, n_peaks)))))
            else:
                BedOverlap(**args).run()
        else:
            log.error('qc_peak_overwrite() failed, peak files not exists')


def qc_bam_fingerprint(x, hiseq_type='rn', bam_type='r1'):
    """
    Parameters
    ----------
    x:  str
        The project dir of hiseq

    hiseq_type:  str
        The hiseq type of `x`, options: ['r1', 'rn', 'rx']
        default: ['rn']

    bam_type:  str
        The hiseq type of peaks, options: ['r1', 'rn', 'rx']
        default: ['r1']
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    if not a.is_hiseq:
        log.error('qc_bam_fingerprint() failed, not a hiseq dir: {}'.format(x))
        return None
    bam_list = list_hiseq_file(x, 'bam', bam_type)
    args = {
        'bam_list': bam_list,
        'outdir': a.qc_dir,
        'prefix': '09.fingerprint',
        'threads': a.threads

    }
    if file_exists(a.bam_fingerprint_png) and not a.overwrite:
        log.info('qc_bam_fingerprint() skipped, file exists')
    else:
        if all(file_exists(bam_list)):
            Bam2fingerprint(**args).run()
        else:
            log.error('qc_bam_fingerprint() failed, peak files not exists')

