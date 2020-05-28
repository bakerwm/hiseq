#!/use/bin/Rscript

# Date: 2019-12-29

# only support ATAC-seq: chrM + genome (dm6)
#
# to-do
# support general options

# create mapping stat plots

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("lendist.R <lendist.txt> <lendist.pdf>")
}

lendist <- args[1]
pdfout  <- args[2]

# # print(here::here())
#
# #func <- here::here("hiseq", "bin", "qc_report_function.R")
# basedir <- "/home/wangming/work/wmlib/hiseq"
# func <- file.path(basedir, "hiseq", "bin", "qc_report_function.R")
# source(func)
suppressPackageStartupMessages(library(hiseqr))


df <- hiseqr::read_frag(lendist)
p  <- hiseqr::frag_plot(df)

pdf(pdfout)
print(p)
dev.off()

##
