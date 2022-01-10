#!/usr/bin/env Rscripts
# run GO analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

err_msg <- "Usage: Rscript run_deseq_pair.R <x> <y> <outdir>

Options:
	deseq_x   deseq_dir (a.vs.b) of groupA
	deseq_y    deseq_dir (a.vs.b) of groupB
	 outdir    outdir

Example:
"

if (length(args) < 3) {
	cat(err_msg)
  stop("arguments failed")
}

# args
x    <- args[1]
y    <- args[2]
outdir  <- args[3]

library(hiseqr)
library(dplyr)
library(ggplot2)


# run pipe, Rmarkdown
hiseqr::deseq_pair_report(x, y, outdir)


# END
