#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript run_deseq.R <config.json>")
  print("")
  print("Options:")
  print("  config.json, include count_ctl, count_exp, genome, outdir")
  stop("arguments failed")
}

# count.txt, path_out,
config_json <- args[1]
args        <- jsonlite::read_json(config_json)

# check required args
required_args <- c("count_ctl", "count_exp", "genome", "outdir")
if(! all(required_args %in% names(args))) {
  msg <- paste0("missing args: ", paste(required_args, collapse = ", "))
  stop(msg)
}

library(hiseqr)
library(dplyr)
library(ggplot2)
library(ggrepel)

hiseqr::deseq_hub2(count_ctl = unlist(args$count_ctl),
                   count_exp = unlist(args$count_exp),
                   genome    = args$genome,
                   outdir    = args$outdir)

## END
