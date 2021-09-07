#!/usr/bin/env Rscripts
# run DESeq2 analysis for salmon directory
#
# Generate plots
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript run_deseq_salmon.R <project_dir>")
  print("")
  print("Options:")
  stop("arguments failed")
}

# count.txt, path_out,
project_dir <- args[1]
run_go <- ifelse(length(args) > 1, args[2], 0)
# project_dir <- "/data/yulab/wangming/work/devel_pipeline/hiseq/rnaseq/aaaaaa"

suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

# run DESeq2
# hiseqr::rnaseq_salmon_hub(project_dir)
output <- file.path(project_dir, "deseq")
hiseqr::hiseq_deseq(project_dir, outdir = output)
if(run_go == 1) {
  genome <- list_hiseq_file(project_dir, "genome", "rx")
  # hiseqr::rnaseq_enrich_hub(project_dir, organism = genome)
}




# pd <- hiseqr::read_hiseq(project_dir)
# # arguments
# wt_quant  <- list_hiseq_file(project_dir, "wt_quant", "rx") # px$wt_quant
# wt_name   <- list_hiseq_file(project_dir, "wt_name", "rx") # pd$args$wt_name
# mut_quant <- list_hiseq_file(project_dir, "mut_quant", "rx") # px$mut_quant
# mut_name  <- list_hiseq_file(project_dir, "mut_name", "rx") # pd$args$mut_name
# deseq_dir <- list_hiseq_file(project_dir, "deseq_dir", "rx")
# genome    <- list_hiseq_file(project_dir, "genome", "rx")
# fix_xls   <- list_hiseq_file(project_dir, "deseq_fix_xls", "rx")
#
# # check tx2gene
# salmon_index <- list_hiseq_file(project_dir, "salmon_index", "rx")
# if(is(salmon_index, "character")) {
#   tx2gene_csv <- file.path(salmon_index, "tx2gene.csv")
#   if(! file.exists(tx2gene_csv)) {
#     on.exit(glue::glue("tx2gene.csv not found: {tx2gene_csv}"))
#   }
# } else {
#   on.exit(glue::glue("salmon_index not found, {project_dir}"))
# }
#
# # prepare data
# sf_list <- c(wt_quant, mut_quant)
# names(sf_list) <- basename(dirname(dirname(sf_list)))
#
# # design
# samples <- data.frame(
#   name      = names(sf_list),
#   condition = c(rep(wt_name,  length = length(wt_quant)),
#                 rep(mut_name, length = length(mut_quant)))
# ) %>%
#   dplyr::mutate(condition = factor(
#     condition, levels = c(wt_name, mut_name)))
#
# # experiment design
# # x <- "~/data/genome/dm6/salmon_index/data/dm6.tx2gene.csv"
# tx2gene <- read.csv(tx2gene_csv)
#
# # load data
# txi <- tximport::tximport(sf_list, type = "salmon", tx2gene = tx2gene)
#
# # DEseq2 analysis
# dds_txi <- DESeq2::DESeqDataSetFromTximport(
#   txi, colData = samples, design = ~condition
# )
#
# # run
# hiseqr::deseq2_main(dds_txi, deseq_dir, organism = genome)
# hiseqr::make_publish_plots(project_dir)
# if(run_go == 1) {
#   hiseqr::rnaseq_enrich_hub(project_dir, organism = genome)
# }


## END
