#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
General modules for CnRseq analysis

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
from hiseq.bam2bw.bam2bw import Bam2bw
from hiseq.cnr.callpeak import CallPeak
from hiseq.fragsize.fragsize import BamPEFragSize
from hiseq.utils.file import list_file, list_dir, check_file, \
    check_path, copy_file, symlink_file, remove_file, fx_name, \
    file_exists, file_abspath, file_prefix, file_nrows
from hiseq.utils.bam import Bam, Bam2cor, Bam2fingerprint
from hiseq.utils.bed import PeakIDR, BedOverlap, PeakFRiP
from hiseq.utils.utils import log, update_obj, Config, get_date, \
    read_hiseq, list_hiseq_file, run_shell_cmd


def cnr_trim(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to CnrR1() dir, cnrseq_type=cnrseq_r1

    option-1: cut reads from 3' end, to X-nt in length (default: 50)

    operation:
    1. do trimming
    2. create symlink

    path -> (CnrReader) -> args
    fq1, fq2
    clean_fq1, clean_fq2, clean_dir
    cut_to_length, recursive
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    # do-the-trimming
    fq1, fq2 = a.raw_fq_list
    clean_fq1, clean_fq2 = a.clean_fq_list
    # whether to trim or not
    if a.trimmed:
        symlink(fq1, clean_fq1)
        symlink(fq2, clean_fq2)
    else:
        args_local = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': a.clean_dir,
            'library_type': None,
            'len_min': 20,
            'cut_to_length': a.cut_to_length,
            'recursive': a.recursive,
            'parallel_jobs': 1 # do not allowed > 1 !!!!
        }
        trim = TrimR1(**args_local)
        trim.run()
        ## copy files
        symlink_file(trim.clean_fq1, clean_fq1)
        symlink_file(trim.clean_fq2, clean_fq2)


def cnr_align_genome(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to CnrR1() dir, cnrseq_type=cnrseq_r1

    Align reads to reference genome, using bowtie2

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
    args_local = a.__dict__
    fq1, fq2 = getattr(a, 'clean_fq_list', [None, None])
    args_init = {
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.align_dir,
        'genome': a.genome,
        'extra_index': getattr(a, 'extra_index', None),
        'aligner': a.aligner,
        'keep_tmp': getattr(a, 'keep_tmp', False),
        'max_fragment': 2000,
        'unique_only': True,
        'extra_para': '-I 10 -X 2000', # specific for CnR
    }
    args_local.update(args_init)
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
            symlink_file(t.align_stat, a.align_stat)
            symlink_file(t.align_json, a.align_json)
            symlink_file(t.align_flagstat, a.align_flagstat)
    # remove PCRdup
    if a.rmdup:
        if file_exists(a.bam):
            if file_exists(a.bam_rmdup) and not a.overwrite:
                log.info('rmdup() skipped, file exists: {}'.format(
                    a.bam_rmdup))
            else:
                Bam(a.bam).rmdup(a.bam_rmdup)
                Bam(a.bam_rmdup).index()
    else:
        symlink_file(a.bam, a.bam_rmdup)
    # calculate norm scale
    if not file_exists(a.align_scale_json):
        df = {
            'bam': a.bam_rmdup,
            'norm': 1e6,
            'scale': cal_norm_scale(a.bam_rmdup),
        }
        Config().dump(df, a.align_scale_json)
#         a.align_scale = cal_norm_scale(a.bam_rmdup)        
#         with open(a.align_scale_txt, 'wt') as w:
#             w.write('{:.4f}\n'.format(a.align_scale))


# to-do: update required
def cnr_align_spikein(x, hiseq_type='_r1'):
    """
    Parameters
    ---------
    x: str
        The path to CnrR1() dir, cnrseq_type=cnrseq_r1

    Align reads to reference genome, using bowtie2

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
    args_local = a.__dict__
    fq1, fq2 = getattr(a, 'clean_fq_list', [None, None])
    args_init = {
        'fq1': fq1,
        'fq2': fq2,
        'outdir': a.spikein_dir,
        'genome': None,
        'index': getattr(a, 'spikein_index', None),
        'aligner': a.aligner,
        'keep_tmp': getattr(a, 'keep_tmp', False),
        'max_fragment': 2000,
        'extra_para': '-X {}'.format(2000),
    }
    args_local.update(args_init)
    if file_exists(a.bam) and not a.overwrite:
        log.info('align() skipped, file exists: {}'.format(a.bam))
    else:
        Align(**args_local).run()
    # copy files
    t = HiseqReader(a.align_dir) #
    symlink_file(t.args.bam, a.spikein_bam)
    symlink_file(t.args.align_stat, a.spikein_stat)
    symlink_file(t.args.align_json, a.spikein_json)
    symlink_file(t.args.align_flagstat, a.spikein_flagstat)
    ## calculate norm scale
    if not file_exists(a.spikein_scale_json):
        df = {
            'bam': a.spikein_bam,
            'norm': 1e4,
            'scale': cal_norm_scale(a.spikein_bam),
        }
        Config().dump(df, a.spikein_scale_json)
#         spikein_scale = cal_norm_scale(a.spikein_bam, 10000) # to 1e4 ?
#         with open(a.spikein_scale_txt, 'wt') as w:
#             w.write('{:.4f}\n'.format(spikein_scale))


def cnr_merge_bam(x, hiseq_type='_r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    # get bam_list
    bam_list = list_hiseq_file(x, 'bam_rmdup', 'r1')
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
    # rmdup
    if a.rmdup:
        if file_exists(a.bam):
            if file_exists(a.bam_rmdup) and not a.overwrite:
                log.info('rmdup() skipped, file exists: {}'.format(
                    a.bam_rmdup))
            else:
                Bam(a.bam).rmdup(a.bam_rmdup)
                Bam(a.bam_rmdup).index()
    else:
        symlink_file(a.bam, a.bam_rmdup)
    # check-point
    if not file_exists(a.bam_rmdup):
        raise ValueError('cnr_merge_bam() failed, see: {}'.format(a.bam_dir))
    # calculate norm scale
    if not file_exists(a.align_scale_json):
        df = {
            'bam': a.bam_rmdup,
            'norm': 1e6,
            'scale': cal_norm_scale(a.bam_rmdup),
        }
        Config().dump(df, a.align_scale_json)


def cnr_call_peak(x, hiseq_type='_r1'):
    """
    Call peaks using MACS2 and SEACR
    -f BAMPE
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    args = {
        'ip': a.bam_rmdup, # rm_dup
        'input': None,
        'outdir': a.peak_dir,
        'prefix': a.smp_name,
        'genome': a.genome,
        'genome_size': a.genome_size,
        'genome_size_file': a.genome_size_file
    }
    CallPeak(method='macs2', **args).run()
    CallPeak(method='seacr', **args).run()


# quality control #
def cal_norm_scale(bam, norm=1000000): # to 1e6
    """
    scale factor
    Bam().count_reads()
    """
    bam_o = Bam(bam)
    is_pe = bam_o.is_paired()
    n_map = bam_o.getNumberOfAlignments()
    if is_pe:
        n_map = n_map/2
    if n_map > 0:
        n_scale = norm/n_map
    else:
        log.error('no mapped reads detected')
        n_scale = 1
    return n_scale


def cnr_bam_to_bw(x, hiseq_type='_r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    args = {
        'bam': a.bam_rmdup if getattr(a, 'rmdup', False) else a.bam,
        'prefix': a.smp_name,
        'outdir': a.bw_dir,
        'binsize': a.binsize,
        'strandness': 0, # non-strandness
        'genome': a.genome,
        'scaleFactor': load_scale(a.align_scale_json),
        'overwrite': a.overwrite,
        'genome_size': a.genome_size,
    }
    Bam2bw(**args).run()


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


def load_scale(x):
    """
    read norm scale from file, convert to float
    default: 1.0

    $ cat scale.json
    {
        'bam': 'a.bam',
        'norm': 1000000,
        'scale': 0.7345,
    }
    """
    try:
        df = Config().load(x)
        if not isinstance(df, dict):
            df = {}
    except:
        log.error('Could not read file: {}'.format(x))
        df = {}
    return df.get('scale', 1.0)

    # get the first line
#     try:
#         with open(x) as r:
#             s = next(r).strip() # 1st line
#     except IOError as e:
#         log.error(e)
#         s = '1.0'
#     # convert to float
#     try:
#         f = float(s)
#     except ValueError as v:
#         log.error(v)
#         f = 1.0
#     return f


################################################################################
## Quality control matrix for CnRseq analysis
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
def qc_trim(x, hiseq_type='r1'):
    """
    # format:
    # name, total, too_short, dup, too_short2, clean, percent
    """
    a = read_hiseq(x, hiseq_type) # for general usage
    # option-1: stat.toml
    # option-2: stat.txt
    stat_json = getattr(a, 'trim_stat_json', None)
    stat_txt = getattr(a, 'trim_stat_txt', None)
    # format:
    # name, total, too_short, dup, too_short2, clean, percent
    d = {
        'id': a.smp_name,
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
                    'id': s[0],
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


def qc_align(x, hiseq_type='r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    # load data
    df = {
        'name': a.smp_name,
        'index': 'genome',
        'total': 1,
        'unique': 1,
        'multi': 0,
        'map': 1,
        'unmap': 0,
        'unique_only': False
    }
    align_json = getattr(a, 'align_json', None)
    align_toml = getattr(a, 'align_toml', None)
    f = align_json if file_exists(align_json) else \
        align_toml if file_exists(align_toml) else None
    if file_exists(f):
        df = Config().load(f)
    # update: nodup, spikein, chrM, ...
    nodup = Bam(a.bam_rmdup)
    is_pe = nodup.is_paired()
    n_nodup = nodup.getNumberOfAlignments()
    if is_pe:
        n_nodup = n_nodup/2
    # spikein
    if file_exists(a.spikein_bam):
        n_spikein = Bam(a.spikein_bam).getNumberOfAlignments()
        if is_pe:
            n_spikein = n_spikein/2
    else:
        n_spikein = 0
    # chrM
    n_mt = get_mito_count(a.bam_rmdup)
    df['nodup'] = n_nodup
    df['spikein'] = n_spikein
    df['chrM'] = n_mt
    Config().dump(df, a.align_summary_json)
    return df


def qc_lendist(x, hiseq_type='_r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    if os.path.exists(a.lendist_txt) and a.overwrite is False:
        log.info('lendist() skipped: file exists: {}'.format(
            a.lendist_txt))
    else:
        BamPEFragSize(a.bam_rmdup).saveas(a.lendist_txt)
        
        
def qc_frip(x, hiseq_type='_r1'):
    a = read_hiseq(x, hiseq_type) # for general usage
    if check_file(a.frip_json, check_empty=True):
        log.info('qc_frip() skipped, file exists: {}'.format(
            a.frip_json))
    else:
        try:
            # dict: index, total, n, frip
            s = PeakFRiP(peak=a.peak, bam=a.bam_rmdup,
                method='featureCounts').run()
            # number of peaks
            s.update({
                'id': a.project_name,
                'n_peaks': file_nrows(a.peak),
            })
            Config().dump(s, a.frip_json) 
        except:
            log.error('PeakFRiP() failed, see: {}'.format(a.frip_json))


def qc_tss_enrich(x, hiseq_type='r1', bw_type='r1'):
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
    # r1
    if a.is_hiseq_r1:
        arg_bw = a.bw
        arg_label = '--samplesLabel {}'.format(a.smp_name)
        per_group = '' # plotProfile
    elif a.is_hiseq_rn:
        bw_list = list_hiseq_file(x, 'bw', 'r1')
        arg_bw = ' '.join(bw_list)
        smp_name = list_hiseq_file(x, 'smp_name', 'r1') # multi
        arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        per_group = '--perGroup' # plotProfile
    else:
        arg_bw = ''
        arg_label = ''
        per_group = ''
    cmd = ' '.join([
        '{}'.format(shutil.which('computeMatrix')),
        'reference-point',
        '-R {}'.format(getattr(a, 'gene_bed', None)),
        '-S {}'.format(arg_bw),
        arg_label,
        '-o {}'.format(a.tss_enrich_matrix),
        '--referencePoint TSS',
        '-b 2000 -a 2000 --binSize 10 --sortRegions descend --skipZeros',
        '-p {}'.format(a.threads),
        '2> {}'.format(a.tss_enrich_matrix_log),
        '&& {}'.format(shutil.which('plotProfile')),
        '-m {}'.format(a.tss_enrich_matrix),
        '-o {}'.format(a.tss_enrich_png),
        '--dpi 300',
        per_group,
        ])
    if file_exists(a.tss_enrich_png) and not a.overwrite:
        log.info('qc_tss() skipped, file exists')
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


def qc_genebody_enrich(x, hiseq_type='r1', bw_type='r1'):
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
    # r1
    if a.is_hiseq_r1:
        arg_bw = a.bw
        arg_label = '--samplesLabel {}'.format(a.smp_name)
        per_group = '' # plotProfile
    elif a.is_hiseq_rn:
        bw_list = list_hiseq_file(x, 'bw', 'r1')
        arg_bw = ' '.join(bw_list)
        smp_name = list_hiseq_file(x, 'smp_name', 'r1') # multi
        arg_label = '--samplesLabel {}'.format(' '.join(smp_name))
        per_group = '--perGroup' # plotProfile
    else:
        arg_bw = ''
        arg_label = ''
        per_group = ''
    cmd = ' '.join([
        '{}'.format(shutil.which('computeMatrix')),
        'scale-regions',
        '-R {}'.format(getattr(a, 'gene_bed', None)),
        '-S {}'.format(arg_bw),
        arg_label,
        '-o {}'.format(a.genebody_enrich_matrix),
        '-b 2000 -a 2000 --regionBodyLength 2000',
        '--binSize 10 --sortRegions descend --skipZeros',
        '--smartLabels',
        '-p {}'.format(a.threads),
        '2> {}'.format(a.genebody_enrich_matrix_log),
        '&& {}'.format(shutil.which('plotProfile')),
        '-m {}'.format(a.genebody_enrich_matrix),
        '-o {}'.format(a.genebody_enrich_png),
        '--dpi 300',
        per_group,
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
    bam_list = list_hiseq_file(x, 'bam_rmdup', bam_type)
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
    bam_list = list_hiseq_file(x, 'bam_rmdup', bam_type)
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
