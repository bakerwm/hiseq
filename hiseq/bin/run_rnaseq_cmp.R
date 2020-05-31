#!/usr/bin/env Rscripts
# run GO analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

err_msg <- "Usage: Rscript run_GO.R <x>

Options:
	deseq_dirA    deseq_dir (a.vs.b) of groupA
	deseq_dirB    deseq_dir (a.vs.b) of groupB
	   feature    gene/te/ ...

Example:
1. Rscript run_GO.R deseq_dir gene 1
2. Rscript run_GO.R gene_list.txt dm6 outdir/
3. Rscript run_GO.R gene_list.txt dm6 outdir/ fc_exp.txt

"

if (length(args) < 4) {
	cat(err_msg)
  stop("arguments failed")
}

# args
dirA    <- args[1]
dirB    <- args[2]
feature <- args[3]
outdir  <- args[4]

library(hiseqr)
library(dplyr)
library(ggplot2)


# run pipe, Rmarkdown

hiseqr::compare_deseq_report(dirA, dirB, feature, outdir)


# END