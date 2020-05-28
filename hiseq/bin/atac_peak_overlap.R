#!/use/bin/Rscript

# Date: 2019-12-29
# overlaps between replicates

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("peak_overlap.R <outdir> <peak1> <peak2> ...")
}

pdfout <- file.path(args[1], "peak_overlap.pdf")
peak_list <- args[-1]

## only 2, 3 samples 
if(length(peak_list) > 3) {
  peak_list <- peak_list[1:3]
} else if(length(peak_list) == 1) {
  return(NULL) # skip
} 

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hiseqr))

names <- paste0("rep", seq_len(length(peak_list)))
p     <- hiseqr::bed_venn(blist = as.list(peak_list), names = names) 

pdf(pdfout)
print(p)
dev.off()

