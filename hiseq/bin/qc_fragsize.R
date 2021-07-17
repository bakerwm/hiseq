#!/use/bin/Rscript

# Date: 2020-09-17
# fragment length distribution

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("lendist.R <lendist.pdf> <lendist.txt> [txt2, ...]")
}

pdfout  <- args[1]
lendist <- args[-1] # multiple files

suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

df <- lapply(lendist, function(f){
  readr::read_csv(f, comment = '#', col_types = readr::cols())
}) %>%
  dplyr::bind_rows()
p <- hiseqr::fragsize_plot(df)

## fix plot height
# 2:4, 4:4.5, 6:5, 8:5.5, 10:6
n_files  <- length(lendist)
n_height <- floor(n_files/2)*.8 + 3.5

pdf(pdfout, width = 7, height = n_height)
print(p)
dev.off()
