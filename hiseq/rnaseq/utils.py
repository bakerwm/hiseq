#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
General modules for RNAseq analysis

analysis-module:
"""

import os
import re
import glob
import shutil
import pysam
import hiseq
import pyfastx
from hiseq.trim.trimmer import TrimR1
from hiseq.align.align import Align
from hiseq.bam2bw.bam2bw import Bam2bw, bw_compare
from hiseq.fragsize.fragsize import BamFragSize
from hiseq.utils.file import list_file, list_dir, check_file, \
    check_path, copy_file, symlink_file, remove_file, fx_name, \
    file_exists, file_abspath, file_prefix, file_nrows
from hiseq.utils.bam import Bam, Bam2cor, Bam2fingerprint
from hiseq.utils.bed import PeakIDR, BedOverlap, PeakFRiP
from hiseq.utils.utils import log, update_obj, Config, get_date, \
    read_hiseq, list_hiseq_file, run_shell_cmd, \
    find_longest_common_str


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
        for _,seq,_,_ in pyfastx.Fastx(fq1):
            i += 1
            if seq.startswith('CT'):
                n1 += 1
            if i >= top_n:
                break
        # fq2
        j = 0
        for _,seq,_,_ in pyfastx.Fastx(fq2):
            j += 1
            if seq.startswith('GA'):
                n2 += 1
            if j >= top_n:
                break
    except IOError as e:
        log.error(e)
    # check
    return n1 / top_n > cutoff and n2 / top_n > cutoff




def rnaseq_trim(x, hiseq_type='r1'):
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
    # do-the-trimming
    fq1, fq2 = a.raw_fq_list
    clean_fq1, clean_fq2 = a.clean_fq_list
    # whether to trim or not
    if a.trimmed:
        symlink(fq1, clean_fq1)
        symlink(fq2, clean_fq2)
    else:
        # guess NSR or not
        is_nsr = guess_nsr(fq1, fq2)
        args_local = {
            'fq1': fq1,
            'fq2': fq2,
            'outdir': a.clean_dir,
            'len_min': 20,
            # 'cut_after_trim': '9,-9',
            'keep_tmp': True,
            'parallel_jobs': 1 # do not allowed > 1 !!!!
        }
        trim = TrimR1(**args_local)
        trim.run()
        ## copy files
        symlink_file(trim.clean_fq1, clean_fq1)
l