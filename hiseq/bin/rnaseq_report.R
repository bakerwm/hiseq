#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript atac_report.R <sample.dir> <out.dir>")
  print("")
  print("Option:")
  print("  input.dir    The directory of RNAseq directory")
  print("    out.dir    The directory to save html file")
  stop("arguments failed")
}

indir   <- args[1]
outdir  <- args[2]

hiseqr::rnaseq_report(indir, outdir)

## EOF
