import hiseq
import os, sys
import atac
import atac2

# rep_list = ['results/pe_rep1_1', 'results/pe_rep2_1']
d = '/data/yulab/wangming/work/yu_2019/projects/20191204_lxh_ATACseq/results/20191230_hiseq_pipeline/results'
rep_list = [os.path.join(d, i) for i in [
	'ATACseq_DaGal4Xsh3893_3h_rep1_1',
	'ATACseq_DaGal4Xsh3893_3h_rep2_1']]

# atac.Atac2(rep_list).run()

rep_list = hiseq.helper.listfiles(d, include_dir=True)
rep_list = [i for i in rep_list if os.path.isdir(i)]

# atac.AtacBatch2(rep_list).run()
config = 'config.txt'
design = 'design.json'
fq1 = "data/pe_rep1_1.fq.gz"
fq2 = "data/pe_rep1_2.fq.gz"
genome = "dm6"
outdir = "results/atac2"

rep_list = ['results/pe_rep1_1', 'results/pe_rep2_1']
input_dir = 'results'

# atac.AtacBatch2(rep_list).run()
# atac2.AtacConfig(design=design)
# atac2.AtacConfig(fq1=fq1, fq2=fq2, genome=genome, outdir=outdir)
# atac2.AtacConfig(rep_list=rep_list)
# atac2.AtacConfig(input_dir=input_dir)

# atac2.Atac(design=design).run()
# atac2.AtacSingle(fq1=fq1, fq2=fq2, genome=genome, outdir=outdir).run()
# atac2.AtacMerge(rep_list=rep_list).run()
# hiseq atac
# atac2.Atac(design=design).run()
atac2.Atac(config=config).run()