#!/usr/bin/env Rscripts

# run hiseq_enrich()
args <- commandArgs(trailingOnly = TRUE)

usage <- "Usage: Rscript run_enrich.R <x> [...]

Options:
  x  path to the rnaseq_rx dir,

Example:
$ Rscript run_enrich.R rnaseq_rx

"

if(length(args) < 1) {
  cat(usage)
  stop("arguments failed")
}

suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(clusterProfiler))

lapply(args, function(x) {
  if(hiseqr::is_hiseq_dir(x, "rnaseq_rx")) {
    hiseqr::hiseq_enrich(x)
  }
})

# EOF
