# hiseq



This is a package for High-throughput sequencing data analysis. including pipelines for `RNAseq`, `ChIPseq`, `ATACseq`, `smRNAseq` and some utilities for regular usage.


## submodules


+ qc    
+ align  
+ peak  
+ motif  
+ annotation  
+ report   
+ RNAseq   
+ ChIPseq   
+ ATACseq   
+ smRNAseq   


## to-do

1. Support `YAML` format  
2. add feature: `waiting`, check status every N minitues for each subfunction  
3. Fix RNAseq, kallisto, salmon, + DESeq2; kallisto + sleuth; ...   





1. `hiseq demx`

  - Check fq, MD5 from the same dir   
  - Record the working lines    
  - Re-run the program, skip-n lines,   
  - Exception, 

2. `hiseq align`   

  - Align with `--para-extra`





requirements 

pyyaml
toml
requests
bs4
pillow
crcmod 
trackhub 
python-Levenshtein 


pyfastx
bowtie
bowtie2
bwa
hisat2
STAR=2.5.2a 


