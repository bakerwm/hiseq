#!/usr/bin/env Rscripts


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript hiseq_kraken2_report.R <sample.dir> <out.dir>")
  print("")
  print("Option:")
  print("  sample.dir     The directory of HiSeq output")
  print("     out.dir     The directory to save html file")
  stop("arguments failed")
}

indir   <- args[1]
outdir  <- args[2]

library(hiseqr)
library(dplyr)
library(ggplot2)

hiseqr::hiseq_kraken2_report(indir, outdir)

## EOF
