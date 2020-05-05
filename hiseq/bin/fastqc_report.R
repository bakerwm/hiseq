#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript qc_report.R <path>")
  print("")
  print("Option:")
  print("  path    The directory of fastqc output files, zip files")
  stop("arguments failed")
}

indir <- args[1]
outdir <- args[2]

if (! dir.exists(indir)) {
  stop("directory not exists")
}

if (! require("devtools")) {
  install.packages("devtools")
}

if (! require(hiseqr)) {
  devtools::install_github("bakerwm/hiseqr")
}

if (! require(fastqcr)) {
  #devtools::install_github("bakerwm/fastqcr")
  install.packages("fastqcr")
}

if (! require(dplyr)) {
  install.packages("dplyr")
}

library(hiseqr)


hiseqr::fastqc_report(indir = indir, outdir = outdir, preview = FALSE)

## save
print(paste0("Saving results in ", indir))

## EOF


