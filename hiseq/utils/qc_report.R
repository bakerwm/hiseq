#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  print("Usage: Rscript qc_report.R <qc.dir>")
  print("")
  print("Option:")
  print("  qc.dir    The directory of fastqc output files")
  stop("arguments failed")
}

qc.dir <- args[1]
template <- NULL
if (file.exists(args[2])) {
  template <- args[2]
}

if (! dir.exists(qc.dir)) {
  stop("directory not exists")
}

if (! require("devtools")) {
  install.packages("devtools")
}

if (! require(goldclipReport)) {
  devtools::install_github("bakerwm/goldclipReport")
}

if (! require(fastqcr)) {
  devtools::install_github("bakerwm/fastqcr")
}

if (! require(dplyr)) {
  install.packages("dplyr")
}

library(goldclipReport)

#qc.report <- file.path(qc.dir, "report")

if (is.null(template)) {
  print("default template")
  # FastqcReport(qc.dir, qc.report, preview = FALSE)
  fastqc_report(qc.dir, qc.dir, preview = FALSE)
} else {
  print("custom template")
  # FastqcReport(qc.dir, qc.report, template = template, preview = FALSE)
  fastqc_report(qc.dir, qc.dir, template = template, preview = FALSE)
}
## save
print(paste0("Saving results in ", qc.dir))

## EOF


