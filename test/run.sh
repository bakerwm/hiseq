# hiseq  qc -i data/pe_rep1_1.fq.gz -o aaa -m 20
# hiseq qc -i data/pe_rep1_1.fq.gz -o ccc -m 20 --rmdup 
hiseq qc -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o aaa --library-type 2
# hiseq qc -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o aaa -m 20 -a AGATCGGAAGAGCACA --threads 4 --rmdup --cut-after-trim 9,-6

 
# hiseq align -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o aaa -g dm3 --aligner bowtie -n demo --unique-only --align-to-rRNA --threads 2

# fq1="data/pe_rep1_1.fq.gz"
# fq2="data/pe_rep1_2.fq.gz"
# outdir="results/align"
# hiseq align -i $fq1 --fq2 $fq2 -o $outdir -g dm3 -k mm10 --aligner bowtie2 -n demo --unique-only --align-to-rRNA --threads 2 --extra-para "-X 2000"

# hipipe-align.py -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o bbb -g dm3 --aligner bowtie2 -n demo --unique-only --align-to-rRNA 



# hiseq align -i data/pe_rep1_1.fq.gz --fq2 data/pe_rep1_2.fq.gz -o results/align -g dm6 --aligner bowtie2 -n demo --extra-para "-X 2000" --n-map 0 --unique-only --align-to-chrM 


