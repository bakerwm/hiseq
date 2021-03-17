from align import *

args = {
    'aligner': 'bowtie',
    'fq1': '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_1.fq.gz',
    'fq2': '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_2.fq.gz',
    # 'index_list': '/home/wangming/data/genome/dm6/bowtie_index/dm6',
    'genome': 'dm6',
    'to_rRNA': True,
#     'spikein': 'mm10',
    'outdir': 'test',
    'genome': 'dm6',
    'threads': 24,
    'unique_only': True,
    'extra_para': '--chunkmbs 128 -X 2000',
}

p = Align(**args).run()