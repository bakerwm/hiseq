from align import *

args = {
    'fq1': '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_1.fq.gz',
    # 'fq2': '/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/data/pe_control_rep1_2.fq.gz',
    'index_list': '/home/wangming/data/genome/dm6/STAR_index/dm6',
    'outdir': 'test',
    'genome': 'dm6',
    'threads': 24,
    'unique_only': True,
}

p = Align(**args)