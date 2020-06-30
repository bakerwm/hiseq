
"""
From ENCODE ChIPSeq pipeline
source: jhttps://github.com/ENCODE-DCC/chip-seq-pipeline2  
Documentation: https://www.encodeproject.org/pages/pipelines/#DNA-binding


Overview:

0. FASTQ read quality filtering    
0a. Crop FASTQ, for 100PE/SE, crop 2 at 5 prim-end, save 98-100 bp read   
1a. Alignment - bowtie2  
  + parameter: --mm  / -X 2000 
  + SAMstats --sorted_sam_file - --outf flagstat.qc
1b. Post-alignment filtering
  + Remove unmapped, mate unmapped, not primary alignment    
  + reads failing platform    
  + Remove low MAPQ reads   
"""


## filtering
## samtools view -F 1804
## picard MarkDuplicates 