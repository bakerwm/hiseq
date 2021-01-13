#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript atac_report_multi.R <project.dir> <out.dir>")
  print("")
  print("Option:")
  print("  prj_dir    The directory of ChIPseq project")
  print("  out_dir    The directory to save html file")
  stop("arguments failed")
}

indir   <- args[1]
outdir  <- args[2]

hiseqr::chipseq_report(indir, outdir)

## EOF
