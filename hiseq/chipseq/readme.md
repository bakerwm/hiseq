
## Introduction of ChIPseq


Terms and Definitions https://www.encodeproject.org/data-standards/terms/


## Library complexity

Criteria: 

NRF: >0.8 for 10M  unique reads


```
PCR Bottlenecking Coefficient 1 (PBC1)

PBC1=M1/MDISTINCT where
M1: number of genomic locations where exactly one read maps uniquely
MDISTINCT: number of distinct genomic locations to which some read maps uniquely
PCR Bottlenecking Coefficient 2 (PBC2)

PBC2=M1/M2 where
M1: number of genomic locations where only one read maps uniquely
M2: number of genomic locations where two reads map uniquely

Non-Redundant Fraction (NRF) – Number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads.


Library Complexity
ChIP-seq Standards:

PBC1    PBC2    Bottlenecking level NRF Complexity  Flag colors
< 0.5   < 1 Severe  < 0.5   Concerning  Orange
0.5 ≤ PBC1 < 0.8    1 ≤ PBC2 < 3    Moderate    0.5 ≤ NRF < 0.8 Acceptable  Yellow
0.8 ≤ PBC1 < 0.9    3 ≤ PBC2 < 10   Mild    0.8 ≤ NRF < 0.9 Compliant   None
≥ 0.9   ≥ 10    None    > 0.9   Ideal   None
ATAC-Seq Standards:

PBC1    PBC2    Bottlenecking level NRF Complexity  Flag colors
< 0.7   < 1 Severe  < 0.7   Concerning  Orange
0.7 ≤ PBC1 ≤ 0.9    1 ≤ PBC2 ≤ 3    Moderate    0.7 ≤ NRF ≤ 0.9 Acceptable  Yellow
> 0.9   > 3 None    > 0.9   Ideal   None

```



**source**   

+ 1 



#### Annotation (ChIPSeeker)


#### ROSE for super-enhancer
http://younglab.wi.mit.edu/super_enhancer_code.html


#### integrate with RNA-seq
[Beta](http://cistrome.org/BETA/) from Shirley Liu's lab in Harvard. Tao Liu's previous lab.


#### Visualization

1. [deeptools](https://github.com/fidelram/deepTools), [ngs.plot](https://github.com/shenlab-sinai/ngsplot)    

2. You can also use bioconductor [Genomation](http://www.bioconductor.org/packages/release/bioc/vignettes/genomation/inst/doc/GenomationManual-knitr.html). It is very versatile.

3. [EnrichedHeatmaps](https://github.com/jokergoo/EnrichedHeatmap) from Zuguang Gu based on his own package ComplexHeatmaps. This is now my default go-to because of the flexiability of the package and the great user support. Thx!




#### remove duplicates

```
sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], outbam = output[0], log = log)
```



#### subtract input bigWig

```
bamCompare --bamfile1 {input[1]} --bamfile2 {input[0]} --normalizeUsing RPKM --ratio subtract --binSize 30 --smoothLength 300 -p 5  --extendReads 200 -o {output} 2> {log}
```




#### make bigWig

```
bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 5 --extendReads 200 -o {output} 2> {log}
```

#### call peaks

```
# macs1
macs14 -t {input.case}              -c {input.control} --keep-dup all -f BAM -g {config[macs_g]}              --outdir 08peak_macs1 -n {params.name1} -p {config[macs_pvalue]} &> {log.macs1}

# macs2
# nomodel for macs14, shift-size will be 100 bp (e.g. fragment length of 200bp)
        # can get fragment length from the phantompeakqual. Now set it to 200 bp for all.
        macs14 -t {input.case}              -c {input.control} --keep-dup all -f BAM -g {config[macs_g]}              --outdir 08peak_macs1 -n {params.name2} --nomodel -p {config[macs_pvalue]} &> {log.macs1_nomodel}

## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.case}              -c {input.control} --keep-dup all -f BAM -g {config[macs2_g]}              --outdir 09peak_macs2 -n {params.name} -p {config[macs2_pvalue]} --broad --broad-cutoff {config[macs2_pvalue_broad]} --nomodel &> {log}
```





## pipeline

**piPipes** 


1. Alignment

```
#to genome:

## SE
bowtie2 --very-sensitive-local -x genome -U fq -q $bowtie2PhredOption -p $CPU 
samtools view -uS -F0x4 

## PE
bowtie2 --very-sensitive-local -X 800 --no-mixed -x genome 
            -1 fq1 -2 fq2 -q $bowtie2PhredOption  -p $CPU 
samtools view -uS -f0x2
```


```
#to genome:

## SE
bowtie -S -v 3 -a -m 100 -best --strata 
samtools view -uS -F0x4 
csem 

## PE
bowtie2 -S -v 3 -a -m 100 --best --strata
samtools view -uS -f 0x2 
samtools view -f0x40 
csem b *bam ...

```


2. Call peaks

```
macs2 callpeak -f bam -t ip.bam -c input.bam -g genome {opt} -B --SPMR 
macs2 bdgcmp -t treat.bdg -c control.bdg -o ppois.bdg -m ppois 
bedGraphToBigWig ppois.bdg gsize ppois.bigWig 
macs2 bdgcmp -t treat.bdg -c control.bdg -o FE.bdg -m FE 
bedGraphToBigWig FE.bdg gsize fE.bigWig 
macs2 bdgcmp -t treat.bdg -c control.bdg -o logLR.bdg -m logLR -p 0.00001 
bedGraphToBigWig logLR.bdg gsize logLR.bigWig
```

Calculate the normalization scale: 1 M effective_depth. 

```
# total tags in treatment: 11963134
# tags after filtering in treatment: 8175251
# maximum duplicate tags at the same position in treatment = 1
# Redundant rate in treatment: 0.32
# total tags in control: 48266953
# tags after filtering in control: 41983618
# maximum duplicate tags at the same position in control = 1
# Redundant rate in control: 0.13
```

the larger one `after filtering in: `.

3. draw figures  

```
bash $DEBUG piPipes_aggregate_bw_on_beds.sh \
    $AGG_DIR \
    $EXT_LEN \
    $BW_OUTDIR/${PREFIX}.ppois.bigWig,$BW_OUTDIR/${PREFIX}.FE.bigWig,$BW_OUTDIR/${PREFIX}.logLR.bigWig && \
    touch .${JOBUID}.status.${STEP}.aggregate_beds
```


### map to TE consensus

```
bowtie2 -q -X 800 --no-mixed -q qual -1 fq1 -2 fq2 -x idx 
bamtobed -i bam | awk -v Q=10 '$5 > Q' > unique.bed

bowtie2 -q -q qual -U fq -X idx
bamtobed -i bam | awk -v Q=10 '$5> Q' > unique.bed

## summary graph
bedtools_piPipes genomecov -i unique.bed    -g $TRANSCRIPTOME_SIZES -bg -scale $NormScaleIP   

## create summarize plots

## normalized counts
```




### output

1. config/parameters      
2. input fq    
3. genome bam, bigWig, bdg, normalized (ppois, logLR, FE)   
4. TE bam, bigWig, bdg, counts, scatter, line, 
5. piRNA clusters: TSS, TES, gene_body
6. plots: 

foldEnrichment



### IDR example

URL: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html  


```
$ idr --samples Pou5f1_Rep1_sorted_peaks.narrowPeak Pou5f1_Rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file Pou5f1-idr \
--plot \
--log-output-file pou5f1.idr.log
```


output format: https://github.com/nboley/idr#output-file-format  
`5-th column` contains IDR value


Generate pseudo replicates, run IDR: `pseudo_rep_idr`

```
#!/bin/sh

# USAGE: sh pseudorep_idr.sh <input BAM rep1> <chip BAM rep1> <input BAM rep2> <chip BAM rep2> <NAME for IDR output>

# This script will take the BAM files and perform the following steps: 
    ## Merge BAMs for ChiP files,
    ## Shuffle reads and split into two new BAM files (pseudo-replicates), 
    ## Merge BAMs for Input files,
    ## Shuffle reads and split into two new BAM files (pseudo-replicates), 
    ## Call peaks on pseudo-replicates with MACS2 , 
    ## Sort peaks called on pseudo-replicates,
    ## IDR analysis using pseudo-replicate peak calls

# Please use the following SLURM directives:
    ## -t 0-12:00
    ## -p short
    ## --mem=40G

date 

inputFile1=`basename $1`
treatFile1=`basename $2`
inputFile2=`basename $3`
treatFile2=`basename $4`
EXPT=$5

NAME1=`basename $treatFile1 _full.bam`
NAME2=`basename $treatFile2 _full.bam`

# Make Directories
mkdir -p /n/scratch2/mm573/idr_chipseq/macs
mkdir -p /n/scratch2/mm573/idr_chipseq/pooled_pseudoreps
mkdir -p /n/scratch2/mm573/idr_chipseq/tmp

# Set paths
baseDir=/n/groups/hbctraining/ngs-data-analysis-longcourse/chipseq/bowtie2
macsDir=/n/scratch2/mm573/idr_chipseq/macs
outputDir=/n/scratch2/mm573/idr_chipseq/pooled_pseudoreps
tmpDir=/n/scratch2/mm573/idr_chipseq/tmp

#Merge treatment BAMS
echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME1}_${NAME2}_merged.bam $baseDir/${treatFile1} $baseDir/${treatFile2}
samtools view -H ${tmpDir}/${NAME1}_${NAME2}_merged.bam > ${tmpDir}/${EXPT}_header.sam

#Split merged treatments
nlines=$(samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}" # This will shuffle the lines in the file and split it
 into two SAM files
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}00 | samtools view -bS - > ${outputDir}/${EXPT}00.bam
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}01 | samtools view -bS - > ${outputDir}/${EXPT}01.bam

#Merge input BAMS
echo "Merging input BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam $baseDir/${inputFile1} $baseDir/${inputFile2}

#Split merged treatment BAM
nlines=$(samtools view ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
samtools view ${tmpDir}/${NAME1}input_${NAME2}input_merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}_input" # This will shuffle the lines in the file and split in two 
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_input00 | samtools view -bS - > ${outputDir}/${EXPT}_input00.bam
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_input01 | samtools view -bS - > ${outputDir}/${EXPT}_input01.bam


#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1 "
macs2 callpeak -t ${outputDir}/${EXPT}00.bam -c ${outputDir}/${EXPT}_input00.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${EXPT}01.bam -c ${outputDir}/${EXPT}_input01.bam -f BAM -g hs -n $macsDir/${NAME2}_pr -B -p 1e-3  2> $macsDir/${NAME2}_pr_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME1}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME1}_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_pr_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${NAME1}_pr_sorted.narrowPeak $macsDir/${NAME2}_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${EXPT}_pseudorep-idr --rank p.value --plot


# Remove the tmp directory
rm -r $tmpDir
```




> An example for our analysis is described below:    
> If starting with < 100K pre-IDR peaks for large genomes (human/mouse): For true replicates and self-consistency replicates an IDR threshold of 0.05 is more appropriate
Use a tighter threshold for pooled-consistency since pooling and subsampling equalizes the pseudo-replicates in terms of data quality. Err on the side of caution and use more stringent IDR threshold of 0.01



