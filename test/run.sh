# hiseq  qc -i data/pe_rep1_1.fq.gz -o aaa -m 20
# hiseq qc -i data/pe_rep1_1.fq.gz -o ccc -m 20 --rmdup 
# hiseq qc -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o aaa --library-type 1
hiseq qc -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o aaa -m 20 -a AGATCGGAAGAGCACA --threads 4 --rmdup --cut-after-trim 9,-6
