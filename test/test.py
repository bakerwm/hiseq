#!/usr/bin/env python3

## test

import hiseq
from alignment import *
from helper import *

## qc
fq1 = 'data/pe_rep1_1.fq.gz'
fq2 = 'data/pe_rep1_2.fq.gz'
outdir = 'results/qc'

## cutadapt
# hiseq.qc.trimmer.Cutadapt(fq1, outdir).run()
# hiseq.qc.trimmer.Cutadapt(fq1, outdir, fq2, cut_after_trim='9,-6').run()

## trimmer
# hiseq.qc.trimmer.Trimmer(fq1, outdir).run()
hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

# outdir = 'results/align'
# 
# args = {
#     'fq': fq1,
#     'fq2': fq2,
#     'outdir': outdir,
#     'aligner': 'bowtie',
#     'genome': 'dm3',
#     'spikein': 'dm6',
#     'smp_name': None,
#     'unique_only': False,
#     'index_parallel': False}
# 
# # hiseq.align.alignment.Alignment(**args).run()
# Alignment(**args).run()
# 
# 
# 
# ## Json
# f = 'pe_rep1.json'
# # print(Json(f).dict)
