#!/usr/bin/Rscript

# Date: 2019-12-29

# calculate the correlation bewteen replicates, multiple BAM
# window = 500bp
# 
# deeptools computeMatrix => matrix
#
# R script => plot (ggcor -> plot)


args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("align.R <count.table> <out.pdf>")
}

countTable <- args[1]
corPdf <- args[2]

# print(here::here())
# func <- here::here("hiseq", "bin", "qc_report_function.R")
basedir <- "/home/wangming/work/wmlib/hiseq"
func <- file.path(basedir, "hiseq", "bin", "qc_report_function.R")
source(func)

p <- corPlot(countTable)

pdf(corPdf)
print(p)
dev.off()