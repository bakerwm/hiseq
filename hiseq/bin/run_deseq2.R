#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript run_deseq2.R <count.txt> <path_plot>")
  print("")
  print("Options:")
  print("    deseq_dir*   the output of featureCounts of Control samples, like: a.vs.b")
  print("      p_value    the pvalue cutoff, default: 0.1")
  stop("arguments failed")
}

# count.txt, path_out,
de_dir <- args[1]
pvalue <- args[2]

suppressPackageStartupMessages(library(goldclipReport))
library(dplyr)

## parse DE_idr for required files
## 1. config/arguments.txt: for the path to required files
config_pickle <- file.path(de_dir, "..", "config", "arguments.pickle")

if(! file.exists(config_pickle)){
  stop("arguments.pickle, not found in input_dir")
}

## read config files
# install.packages('reticulate')
library(reticulate)
pd <- import("pandas")
pd_data <- pd$read_pickle(config_pickle)

## required
required_names <- c("count_ctl", "count_exp", "prefix_ctl",
                    "prefix_exp", "deseqdir", "genome")

for(f in required_names){
  tag = ifelse(f %in% names(pd_data), pd_data[f], "failed")
  print(paste0(f, " : ", pd_data[f]))
}
stopifnot(all(required_names %in% names(pd_data)))

##----------------------------------------------------------------------------##
## custom: a.vs.b
outdir1 <- file.path(pd_data$deseqdir, "control.vs.treatment")
hiseqr::deseqHub2(count_ctl = pd_data$count_ctl,
                  count_exp = pd_data$count_exp,
                  outdir    = outdir1)

## publish quality figures
cnt_fix1 <- file.path(outdir1, "transcripts_deseq2.fix.xls")
stopifnot(file.exists(cnt_fix1))
print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", cnt_fix1))
tmp <- DESeq2_publish_plot(cnt_fix1, outdir1, save2pdf = TRUE)


##----------------------------------------------------------------------------##
## custom: b.vs.a
outdir2 <- file.path(pd_data$deseqdir, "treatment.vs.control")
hiseqr::deseqHub2(count_ctl = pd_data$count_exp,
                  count_exp = pd_data$count_ctl,
                  outdir    = outdir2)

## publish quality figures
cnt_fix2 <- file.path(outdir2, "transcripts_deseq2.fix.xls")
stopifnot(file.exists(cnt_fix2))
print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", cnt_fix2))
tmp <- DESeq2_publish_plot(cnt_fix2, outdir2, save2pdf = TRUE)


## END
