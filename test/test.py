#!/usr/bin/env python3

## test

import hiseq
from hiseq.utils.helper import *
import trimmer
# from hiseq.utils import *
# from call_peak import *
# from rep_cor import *
# from atac import Atac2
import atac

# # input
fq1 = 'data/pe_rep1_1.fq.gz'
fq2 = 'data/pe_rep1_2.fq.gz'
# outdir = 'results/atac'
outdir = 'aaa'
genome = 'dm6'

args = {'fq1': fq1, 'fq2': fq2, 'outdir': outdir}
trimmer.Trimmer(**args).run()



# # config = AtacConfig(fq1, fq2, genome, outdir)

# config = "config_atac_template.json"
# Atac(config).run()



# args = {
# 	'fq1': ['data/pe_rep1_1.fq.gz'], 
# 	'fq2': ['data/pe_rep1_2.fq.gz'],
# 	'genome': 'dm6',
# 	'outdir': 'results/atac3',
# 	'overwrite': True,
# 	'len_min': 18,
# 	'smp_name': 'aaaaaa'}

# args['config'] = 'config.txt' 

# a = AtacBatch(**args).run()
# 
# bb = '/data/yulab/wangming/work/yu_2019/projects/20191204_lxh_ATACseq/results/20191230_hiseq_pipeline/results/'
# args = {}
# rep_list = [
#     bb + '/ATACseq_DaGal4Xsh3893_3h_rep1_1',
#     bb + '/ATACseq_DaGal4Xsh3893_3h_rep2_1']
# args['smpname'] = 'demo'
# # args['outdir'] = None
# atac.Aatc2(rep_list)
 

# ## qc
# fq1 = 'data/pe_rep1_1.fq.gz'
# fq2 = 'data/pe_rep1_2.fq.gz'
# outdir = 'results/qc'

## cutadapt
# hiseq.qc.trimmer.Cutadapt(fq1, outdir).run()
# hiseq.qc.trimmer.Cutadapt(fq1, outdir, fq2, cut_after_trim='9,-6').run()

## trimmer
# args = {'fq1': fq1, 'outdir': outdir}
# args = {'fq1': fq1, 'fq2': fq2, 'outdir': outdir}
# args = {'fq1': fq1, 'fq2': fq2, 'outdir': outdir, 'cut_after_trim': '9,-6'}

# hiseq.qc.trimmer.Trimmer(**args).run()
# trimmer.Trimmer(**args).run()
# hiseq.qc.trimmer.Trimmer(fq1, outdir, fq2, cut_after_trim='9,-6').run()

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

# ## Peak
# bam = 'results/align/demo/2.genome/demo.bam'
# genome = "dm6"
# output = "results/peak"
# prefix = "demo"
# atac = True

# # Macs2(bam, genome, output, prefix, atac=True).callpeak()
# m = Macs2(bam, genome, output, prefix, atac=True)
# p = m.callpeak()
# print(p)

# p = m.annotation()
# print(p)

# Bam(bam).index()

# peak="results/peak/demo_peaks.narrowPeak"
# f = cal_FRiP(peak, bam)
# print(f)


# bam = 'ATACseq_DaGal4Xsh3893_3h_rep1.bam'
# b1 = 'aaa.rmdup.bam'
# b2 = 'aaa.proper_pair.bam'
# b1 = Bam(bam).rmdup('aaa.rmdup.bam')
# b2 = Bam(b1).proper_pair('aaa.proper_pair.bam')

# # call peaks
# Macs2(b2, "dm6", "results/peak", "aaa", atac=True).callpeak()

# peak = "results/peak/aaa_peaks.narrowPeak"
# f = cal_FRiP(peak, b2)
# print(f)

# len_dist
# x = frag_length(b2, "aaa.len_dist.txt")
# # print(x)

# bed_list = ['x.bed', 'y.bed']
# tiffout = 'aaa.tiff'
# x = peak_overlap2(bed_list, tiffout)
# print(x)

# bed_list = ['demo_rep1_peaks.narrowPeak', 'demo_rep2_peaks.narrowPeak']
# outdir = 'qc'
# # peak_idr(bed_list, outdir)
# peak_overlap2(bed_list, outdir)
