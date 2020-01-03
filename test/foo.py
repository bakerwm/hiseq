import hiseq
import os, sys
import atac

# rep_list = ['results/pe_rep1_1', 'results/pe_rep2_1']
d = '/data/yulab/wangming/work/yu_2019/projects/20191204_lxh_ATACseq/results/20191230_hiseq_pipeline/results'
# rep_list = [os.path.join(d, i) for i in [
# 	'ATACseq_DaGal4Xsh3893_3h_rep1_1',
# 	'ATACseq_DaGal4Xsh3893_3h_rep2_1']]

# atac.Atac2(rep_list).run()

rep_list = hiseq.helper.listfiles(d, include_dir=True)
rep_list = [i for i in rep_list if os.path.isdir(i)]

atac.AtacBatch2(rep_list).run()

