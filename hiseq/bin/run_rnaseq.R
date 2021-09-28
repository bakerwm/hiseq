#!/usr/bin/env Rscripts
# run RNAseq pipeline for control.vs.treatment
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print("Usage: Rscript run_rnaseq.R <path>")
  print("")
  print("Options:")
  print("    path   The directory of a.vs.b")
  print("  run_go   Run enrich analysis or not, 1=yes, 0=no, default:[1]")
  stop("arguments failed")
}

rx_dir <- args[1]
run_go <- as.character(args[2])

library(hiseqr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(fishualize)
library(ggthemes)
library(clusterProfiler)

# hiseqr::rnaseq_hub(args[1]) # ctl_vs_exp = TRUE)
hiseqr::hiseq_deseq(args[1])
if(run_go == 1) {
  # hiseqr::rnaseq_enrich_hub(args[1])
  hiseqr::hiseq_enrich(args[1])
}

## END
