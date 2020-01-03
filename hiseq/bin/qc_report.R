#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript qc_report.R <qc.dir> <out.dir>")
  print("")
  print("Option:")
  print("  qc.dir     The directory of fastqc output files")
  print("  out.dir    The directory to save html file")
  stop("arguments failed")
}

qc_dir   <- args[1]
out_dir  <- args[2]
template <- "/data/yulab/wangming/work/wmlib/hiseq/hiseq/bin/qc_report_single.Rmd"

qc_report <- function(input, output, template = NULL, preview = TRUE) {
  ## input
  input <- normalizePath(input)

  if (is.null(template)) {
    report_template <- system.file("extdata",
                                   "fastqc_report_001.Rmd",
                                   package = "goldclipReport")
  } else {
    report_template <- template
  }
  
  ## output
  output <- normalizePath(output)
  outhtml <- file.path(output, "qc_report.html")
  if(! dir.exists(output)) dir.create(output, recursive = TRUE)

  rmarkdown::render(input = report_template, output_file = outhtml,
                    params = list(qc_dir = input))
  
  if (preview) {
    utils::browseURL(outhtml)
  }
}

qc_report(qc_dir, out_dir, template)

## save
print(paste0("Saving results in ", qc_dir))

## EOF


