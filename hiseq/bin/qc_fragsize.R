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
  # %>%
  #   dplyr::filter(length < 1000) %>%
  #   dplyr::mutate(count = ifelse(length < 50, 0, count))
}) %>%
  dplyr::bind_rows()

p <- hiseqr::fragsize_plot(df)
pdf(pdfout, width = 6, heigh = 4)
print(p)
dev.off()
