#!/use/bin/Rscript

# Date: 2020-09-17
# fragment length distribution

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("lendist.R <lendist.pdf> <lendist.txt> [txt2, ...]")
}

pdfout  <- args[1]
lendist <- args[-1] # multiple files

library(readr)
library(dplyr)
library(ggplot2)
suppressPackageStartupMessages(library(hiseqr))

df_list <- lapply(lendist, function(f){
  readr::read_csv(f, comment = '#') %>%
    dplyr::filter(length < 1000)
})

df <- dplyr::bind_rows(df_list)

# df <- readr::read_csv(lendist) %>%
#   dplyr::filter(length < 1000)
p  <- hiseqr::frag_plot2(df)

pdf(pdfout, width = 10, heigh = 6)
print(p)
dev.off()
