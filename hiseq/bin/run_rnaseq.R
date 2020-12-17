#!/usr/bin/env Rscripts
# run RNAseq pipeline for control.vs.treatment
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript run_rnaseq.R <path>")
  print("")
  print("Options:")
  print("    path   The directory of a.vs.b")
  stop("arguments failed")
}

library(hiseqr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(fishualize)
library(ggthemes)

hiseqr::rnaseq_hub(args[1]) # ctl_vs_exp = TRUE)

## END
