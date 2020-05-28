#!/usr/bin/env Rscripts
# run GO analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

err_msg <- "Usage: Rscript run_GO.R <x>

Options:
   deseq_dir* | gene_list*  deseq_dir (a.vs.b) or gene list, one per line, column name: Gene, log2FoldChange
     feature  |  organism*  gene/te or genome, eg: Drosophila melanogaster
  ctl_vs_exp  |    outdir*  1/0 or output dir

Caution:
  gene_list.txt, should be tab-separated, contain Gene or id column
  fc_exp.txt,    should be tab-separated, contain log2FoldChange column

Example:
1. Rscript run_GO.R deseq_dir gene 1
2. Rscript run_GO.R gene_list.txt dm6 outdir/
3. Rscript run_GO.R gene_list.txt dm6 outdir/ fc_exp.txt

"


if (length(args) < 1) {
	cat(err_msg)
  stop("arguments failed")
}

# args
x <- args[1]

library(hiseqr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)

# run type
x_type <- hiseqr::is_hiseq_dir(x)

if(is.null(x_type)) {
  # genes in table
  # require organism, outdir
  if(length(args) < 3) {
    message("require: x, organism, outdir, for go analysis")
    stop("exit...")
  }

  ## for input list
  df <- readr::read_delim(x, "\t", col_types = readr::cols())
  if("Gene" %in% colnames(df)) {
    gene_list <- df$Gene
  } else if("id" %in% colnames(df)) {
    gene_lsit <- df$id
  } else {
    stop(paste0("Gene|id not found in file: ", x))
  }

  ## args
  organism <- args[2]
  outdir   <- args[3]

  ## for foldChange
  foldChange <- NULL # default
  if(length(args) > 3) {
    df_fc <- readr::read_delim(args[4], "\t", col_types = readr::cols())
    if("log2FoldChange" %in% colnames(df_fc)) {
      foldChange <- setNames(df_fc$log2FoldChange, nm = df_fc$Gene)
      foldChange <- sort(foldChange, decreasing = TRUE)
    }
  }

  ## run GO
  hiseqr::go_pipe(gene_list, organism, outdir, foldChange)

} else if(x_type == "deseq_single") {
  feature = args[2] # gene, te
  ctl_vs_exp = args[3] # 1, 0

  if(is.na(feature)) feature <- "gene"
  if(is.na(ctl_vs_exp)) ctl_vs_exp <- TRUE
  # check
  if(is.character(ctl_vs_exp) & ctl_vs_exp == "0") {
    ctl_vs_exp <- FALSE
  } else {
    ctl_vs_exp <- TRUE
  }
  # run GO
  hiseqr::go_pipe(x, feature = feature, ctl_vs_exp = ctl_vs_exp)
} else {
  stop(paste0("unknown input: ", x))
}




