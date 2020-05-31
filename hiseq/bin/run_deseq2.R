#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print("Usage: Rscript run_deseq2.R <count.txt> <path_plot>")
  print("")
  print("Options:")
  print("    deseq_dir*   the output of featureCounts of Control samples, like: a.vs.b")
  print("      feature*   gene|te|...")
  stop("arguments failed")
}

# count.txt, path_out,
deseq_dir <- args[1]
feature   <- args[2]

library(hiseqr)
library(dplyr)
library(ggplot2)
library(ggrepel)

##----------------------------------------------------------------------------##
## custom: a.vs.b
hiseqr::rnaseq_pipe(deseq_dir, feature) # , ctl_vs_exp = TRUE)
#
# ##----------------------------------------------------------------------------##
# ## custom: b.vs.a
# hiseqr::rnaseq_pipe(deseq_dir, feature, ctl_vs_exp = FALSE)

## END
