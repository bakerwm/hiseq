#!/use/bin/Rscript

# Date: 2019-12-29

# only support ATAC-seq: chrM + genome (dm6)
#
# to-do
# support general options

# create mapping stat plots

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("lendist.R <lendist.txt> <lendist.pdf> {bar|line, x_max, x_min}")
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

## check type, x_min, x_max
plot_type <- ifelse(is.na(args[3]), "line", args[3])
x_min     <- ifelse(is.na(args[4]), 0, args[4])
x_max     <- ifelse(is.na(args[5]), 1000, args[5])

df <- hiseqr::read_frag(lendist)
p  <- hiseqr::frag_plot(df, type = plot_type,
                        x_min = x_min, x_max = x_max)

pdf(pdfout)
print(p)
dev.off()

##
