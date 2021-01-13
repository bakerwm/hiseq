#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript atacseq_report.R <sample.dir> <out.dir>")
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
library(ggthemes)
library(ggrepel)
library(patchwork)

hiseqr::hiseq_report(indir, outdir)

## EOF
