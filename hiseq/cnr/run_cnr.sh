#!/bin/bash

################################################################################
# This script was designed to run ChIPseq analysis for IP and Input samples.   #
# support: ChIPseq, CUT&RUN, CUT&Tag reads
# libtype: SE or PE
# input:   required for ChIPseq analysis
# genome: dm6,mm10,hg38,hg19
# 
#
#
# genome data path were fixed !
################################################################################


## Arguments ##
[[ $# -lt 4 ]] && echo "run_cnr.sh <outdir> <fq1> <fq2>" && exit 1
genome=$1 # dm6, mm10, hg38, hg19
outdir=$2
fq1=$3
fq2=$4
prefix=$(basename ${fq1/_1.fq.gz})
CPU=4


## 0. create output dirs
outdir="${outdir}/${prefix}"
mkdir -p $outdir


## Global Variables ##
bowtie2_index_te="/data/biodata/genome/dm6/bowtie2_index/dm6_transposon"
bowtie2_index_genome="/data/biodata/genome/dm6/bowtie2_index/dm6"
CPU=4


exit_code=0
case $genome in
    dm6)
        bowtie2_index_te="/data/biodata/genome/dm6/bowtie2_index/dm6_transposon";
        bowtie2_index_genome="/data/biodata/genome/dm6/bowtie2_index/dm6"
        ;;
    mm10)
        bowtie2_index_te="/data/biodata/genome/dm6/bowtie2_index/dm6_transposon"; # fake TE
        bowtie2_index_genome="/data/biodata/genome/mm10/bowtie2_index/mm10"
        ;;
    hg19)
        bowtie2_index_te="/data/biodata/genome/dm6/bowtie2_index/dm6_transposon"; # fake TE
        bowtie2_index_genome="/data/biodata/genome/hg19/bowtie2_index/hg19"
        ;;
    hg38)
        bowtie2_index_te="/data/biodata/genome/dm6/bowtie2_index/dm6_transposon"; # fake TE
        bowtie2_index_genome="/data/biodata/genome/hg38/bowtie2_index/hg38"
        ;;
    *)
        echo "unknown genome: $genome";
        exit_code=1
        ;;
esac

[[ $exit_code -gt 0 ]] && echo "illegal arguments, exit" && exit 1


## 1. trimming reads
## cut the 3' adapter (TruSeq: AGATCGGAAGAGC, Nextera: CTGTCTCTTATACACATCT)
trim_log="${outdir}/${prefix}.trim.log"
clean_fq1="${outdir}/${prefix}.trim.1.fq.gz"
clean_fq2="${outdir}/${prefix}.trim.2.fq.gz"
echo "1. trimming reads"
[[ -f $clean_fq1 && -f $clean_fq2 ]] && \
    echo "file exists, trimming skipped" || \
    cutadapt -j $CPU -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20 -q 20 -O 3 -e 0.1 -n 3 --trim-n --max-n=0.1 -o $clean_fq1 -p $clean_fq2 $fq1 $fq2 > $trim_log


## 2. align to genome/te
## to TE, and transposon
## insert size: -I 10 -X 700 
## require proper paired reads: --no-mixed --no-discordant --dovetail
## suppress supplementary alignments by FLAG 2048
unmap_te="${outdir}/${prefix}.not_te.fq"
align_sam_te="${outdir}/${prefix}.te.sam"
align_bam_te_raw="${outdir}/${prefix}.te.raw.bam"
align_bam_te="${outdir}/${prefix}.te.bam"
align_log_te="${outdir}/${prefix}.te.log"
echo "2. align to Transposon consensus (TE)"
[[ -f $align_bam_te ]] && \
    echo "file exists, align to TE skipped" || \
    (bowtie2 --mm -p $CPU --local --very-sensitive --no-mixed --no-discordant --dovetail -I 10 -X 700 -x $bowtie2_index_te --no-unal --un-conc $unmap_te -1 $clean_fq1 -2 $clean_fq2 1>$align_sam_te 2> $align_log_te && \
    samtools view -Sub <(samtools view -H $align_sam_te; samtools view -F 2048 $align_sam_te | grep 'YT:Z:CP') | \
    samtools sort -o $align_bam_te_raw - && \
    samtools index $align_bam_te_raw && \
    sambamba markdup -r -t 2 --overflow-list-size 1000000 $align_bam_te_raw $align_bam_te)
[[ -f $align_sam_te ]] && rm $align_sam_te
[[ -f $align_bam_te_raw ]] && rm $align_bam_te_raw


## 3. align to genome: genome
not_te_fq1="${outdir}/${prefix}.not_te.1.fq"
not_te_fq2="${outdir}/${prefix}.not_te.2.fq"
align_sam_genome="${outdir}/${prefix}.genome.sam"
align_bam_genome_raw="${outdir}/${prefix}.genome.raw.bam"
align_bam_genome="${outdir}/${prefix}.genome.bam"
align_log_genome="${outdir}/${prefix}.genome.log"
echo "3. align to genome"
[[ -f $align_bam_genome ]] && \
    echo "file exists, align to genome skipped" || \
    (bowtie2 --mm -p $CPU --local --very-sensitive --no-mixed --no-discordant --dovetail -I 10 -X 700 -x $bowtie2_index_genome --no-unal -1 $not_te_fq1 -2 $not_te_fq2 1>$align_sam_genome 2> $align_log_genome && \
    samtools view -Sub <(samtools view -H $align_sam_genome; samtools view -F 2048 $align_sam_genome | grep 'YT:Z:CP') | \
    samtools sort -o $align_bam_genome_raw - && \
    samtools index $align_bam_genome_raw && \
    sambamba markdup -r -t 2 --overflow-list-size 1000000 $align_bam_genome_raw $align_bam_genome)
[[ -f $align_sam_genome ]] && rm $align_sam_genome
[[ -f $align_bam_genome_raw ]] && rm $align_bam_genome_raw
[[ -f $not_te_fq1 ]] && rm $not_te_fq1
[[ -f $not_te_fq2 ]] && rm $not_te_fq2


## 4. scale factor
## TE: scale factor: te / (te + genome)
echo "4. calculate scale factor"
chr_count_te="${outdir}/${prefix}.te.chr_count.txt"
chr_count_genome="${outdir}/${prefix}.genome.chr_count.txt"
[[ -f $chr_count_te ]] || samtools idxstats $align_bam_te | awk '{print $1,$3}' > $chr_count_te
[[ -f $chr_count_genome ]] || samtools idxstats $align_bam_genome | awk '{print $1,$3}' > $chr_count_genome
total_reads_te="${outdir}/${prefix}.te.total_reads.txt"
total_reads_genome="${outdir}/${prefix}.genome.total_reads.txt"
[[ -f $read_count_te ]] || samtools view -c $align_bam_te > $total_reads_te
[[ -f $read_count_genome ]] || samtools view -c $align_bam_genome > $total_reads_genome
scale_te=$(paste $total_reads_te $total_reads_genome | awk '{printf "%.10f",$1/($1+$2)}')
scale_genome=$(paste $total_reads_te $total_reads_genome | awk '{printf "%.10f",$2/($1+$2)}')


## 5. convert bam to bigwig
echo "5. convert bam to bigWig"
bw_te="${outdir}/${prefix}.te.norm.cpm.bigWig"
bw_log_te="${outdir}/${prefix}.te.bamcoverage.log"
bw_genome="${outdir}/${prefix}.genome.norm.cpm.bigWig"
bw_log_genome="${outdir}/${prefix}.genome.bamcoverage.log"
[[ -f $align_bam_te && ! -f $bw_te ]] &&
    bamCoverage -b $align_bam_te -bs 10 -p $CPU --scaleFactor $scale_te --normalizeUsing CPM -o $bw_te 2>$bw_log_te

[[ -f $align_bam_genome && ! -f $bw_genome ]] &&
    bamCoverage -b $align_bam_genome -bs 10 -p $CPU --scaleFactor $scale_genome --normalizeUsing CPM -o $bw_genome 2>$bw_log_genome


## 6. call peaks
## 7. IDR, cor, ...

