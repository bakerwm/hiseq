#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  print("Usage: Rscript atacseq_report.R <prj.dir>")
  print("")
  print("Option:")
  print("  sample.dir     The directory of hiseq p7 output")
  stop("arguments failed")
}

indir   <- args[1]

library(hiseqr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(patchwork)

hiseqr::hiseq_p7_report(indir)

## EOF
