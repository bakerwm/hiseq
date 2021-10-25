#!/usr/bin/env Rscripts

#Convert csv to html

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript csv_to_html.R <input.csv> <out.html>")
  print("")
  print("Option:")
  print("  input.csv    The directory of HiSeq output")
  print("output.html    The directory to save html file")
  stop("arguments failed")
}

input  <- args[1]
output <- args[2]

library(hiseqr)
print(glue::glue("input: {input}, output: {output}"))
hiseqr::csv_to_html_report(input, output)

## EOF
