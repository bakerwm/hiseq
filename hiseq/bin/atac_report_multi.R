#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript atac_report_multi.R <project.dir> <out.dir>")
  print("")
  print("Option:")
  print("  qc.dir     The directory of ATAC project")
  print("  out.dir    The directory to save html file")
  stop("arguments failed")
}

indir   <- args[1]
outdir  <- args[2]

hiseqr::atac_report(indir, outdir)

## EOF
