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
  print("       p_value   the pvalue cutoff, default: 0.1")
  stop("arguments failed")
}

# count.txt, path_out,
de_dir <- args[1]
pvalue <- args[2]

suppressPackageStartupMessages(library(goldclipReport))

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
# ## custom
deseqHub3(countA = pd_data$count_ctl,
          countB = pd_data$count_exp,
          organism = pd_data$genome,
          nameA  = pd_data$prefix_ctl,
          nameB  = pd_data$prefix_exp,
          outdir = pd_data$deseqdir,
          pvalue_cutoff = 0.1,
          readable=TRUE)


#-----------------------------------------------------------------------------##
## Generate publish quality figures
cnt_fix <- file.path(pd_data$deseqdir, "transcripts_deseq2.fix.xls")
stopifnot(file.exists(cnt_fix))
print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", cnt_fix))
tmp <- DESeq2_publish_plot(cnt_fix, pd_data$deseqdir, save2pdf = TRUE)


## END
