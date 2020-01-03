#!/use/bin/Rscript

# Date: 2019-12-29

# only support ATAC-seq: chrM + genome (dm6)
#
# to-do
# support general options

# create mapping stat plots

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("align.R <align.dir> <out.pdf>")
}

alignDir <- args[1]
alignPdf <- args[2]

# print(here::here())
# func <- here::here("hiseq", "bin", "qc_report_function.R")
basedir <- "/home/wangming/work/wmlib/hiseq"
func <- file.path(basedir, "hiseq", "bin", "qc_report_function.R")
source(func)

# stat files
flist   <- list.files(alignDir, "*.align.txt", all.files = TRUE, 
	full.names = TRUE, recursive = TRUE)
df <- alignStat(flist) %>% dplyr::select(-mito.pct)
p  <- alignPlot(df)

# save to file
rd <- paste0(alignPdf, '.Rdata')
save(df, p, file = rd)

pdf(alignPdf)
print(p)
dev.off()

##
