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

suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

df <- hiseqr::fragReader(lendist)
df$id <- basename(dirname(dirname(lendist)))
p  <- hiseqr::frag_plot(df)

pdf(pdfout)
print(p)
dev.off()

##
