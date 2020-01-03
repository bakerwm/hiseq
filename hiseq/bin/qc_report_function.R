
##----------------------------------------------------------------------------##
## VennDiagram
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

narrowPeakReader <- function(x) {
  ext <- c(signalValue = "numeric", pValue = "numeric",
           qValue = "numeric", peak = "integer")
  rtracklayer::import(x, format = "BED", extraCols = ext)
}

bed_intersect2 <- function(blist){
  stopifnot(length(blist) >= 2)
  bed1 <- blist[[1]]
  bed2 <- blist[[2]]
  # bed 2 gr
  gr1 <- narrowPeakReader(bed1)
  gr2 <- narrowPeakReader(bed2)
  gr12 <- findOverlaps(gr1, gr2)
  
  # intersect
  n1 <- paste("a", seq_len(length(gr1) - length(gr12)), sep = "")
  n2 <- paste("b", seq_len(length(gr12)), sep = "")
  n3 <- paste("c", seq_len(length(gr2) - length(gr12)), sep = "")
  
  ## list
  x <- list(rep1 = c(n1, n2), rep2 = c(n2, n3))
  # out
  return(x)
}
 
bed_intersect3 <- function(blist){
  stopifnot(length(blist) >= 3)
  bed1  <- blist[[1]]
  bed2  <- blist[[2]]
  bed3  <- blist[[3]]
  # bed to gr
  gr1   <- narrowPeakReader(bed1)
  gr2   <- narrowPeakReader(bed2)
  gr3   <- narrowPeakReader(bed3)
  # intersect
  gr12  <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  gr13  <- findOverlaps(gr1, gr3, ignore.strand = TRUE)
  gr23  <- findOverlaps(gr2, gr3, ignore.strand = TRUE)
  gr123 <- findOverlaps(gr12, gr3, ignore.strand = TRUE)
  # numbers
  n12   <- length(gr12)
  n13   <- length(gr13)
  n23   <- length(gr23)
  #
  n123  <- length(gr123)
  n12   <- n12 - n123
  n13   <- n13 - n123
  n23   <- n23 - n123
  n1    <- length(gr1) - n12 - n13 - n123
  n2    <- length(gr2) - n12 - n23 - n123
  n3    <- length(gr3) - n13 - n23 - n123
  # overlap
  out <- c(n1, n2, n3, n12, n13, n23, n123)
  names(out) <- c("n1", "n2", "n3", "n12", "n13", "n23", "n123")
  # out list
  p <- lapply(seq_len(7), function(i){
    paste(letters[i], seq_len(out[i]), sep = "")
  })
  names(p) <- names(out)
  # combine
  x <- list(
    rep1 = c(p$n1, p$n12, p$n13, p$n123),
    rep2 = c(p$n2, p$n12, p$n23, p$n123),
    rep3 = c(p$n3, p$n13, p$n23, p$n123))
  return(x)
}

bed_intersect4 <- function(blist){
  stopifnot(length(blist) >= 4)
  bed1 = blist[[1]]
  bed2 = blist[[2]]
  bed3 = blist[[3]]
  bed4 = blist[[4]]
  # bed to gr
  
  
}

vennplot <- function(x, names = NULL){
  # check category names
  if(is.null(names)) {
    names <- names(x)
  }
  # further
  if(is.null(names)) {
    names <- letters[seq_len(length(x))]
  }
  # main
  p <- ggVennDiagram::ggVennDiagram(x, label = "count", label_alpha = 0, category.names = names) +
    guides(fill = FALSE) +
    scale_fill_gradient(low = "white", high = "firebrick")
  return(p)
}

bed_venn <- function(blist, names = NULL){
  if(length(blist) == 2){
    x <- bed_intersect2(blist)
  } else if(length(blist) == 3) {
    x <- bed_intersect3(blist)
  } else if(length(blist) == 4) {
    x <- bed_intersect4(blist)
  } else {
    stop("only accept narrowpeaks: 2-4 files")
  }
  
  p <- vennplot(x, names)
  
  return(p)
}


##----------------------------------------------------------------------------##
## reads stat
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


readAlign1 <- function(x){
  hd <- c("total", "unmap", "unique", "multi", "map", "sample", "index")
  df <- readr::read_delim(x, "\t", col_names = hd, 
                          col_types = readr::cols(), comment = "#")
  df_chrM   <- df %>% 
    filter(grepl("chrM", index)) %>%
    select(sample, total, unique, multi)
  df_genome <- df %>% 
    filter(grepl("genome", index)) %>%
    select(sample, unique, multi, unmap)
  df_out <- merge(df_chrM, df_genome, by = "sample")
  names(df_out) <- c("sample", "total", "mito.u",
                     "mito.m", "dm6.u", "dm6.m",
                     "unmap")
  return(df_out)
}


readAlign2 <- function(x){
  df <- readr::read_delim(x, ",", col_names = TRUE, 
                          col_types = readr::cols()) %>% 
    dplyr::select(1, 2, 5, 6, 9:11)
  colnames(df) <- c("sample", "total", "mito.u", "mito.m",
                    "dm6.u", "dm6.m", "unmap")
  return(df)
}


alignStat <- function(flist){
  # check format
  if(all(grepl("*.align.txt", flist))) {
    tmp <- lapply(flist, readAlign1)
  } else if(all(grepl("*mapping_stat.csv", flist))) {
    tmp <- lapply(flist, readAlign2)
  } else {
    stop("File type unknown: *.align.txt, or *mapping_stat.csv")
  }
  df  <- dplyr::bind_rows(tmp) %>%
    mutate(mito.pct = (mito.u + mito.m) / total)
  return(df)
}


alignPlot <- function(df) {
  # plot
  groups <- c("mito.u", "mito.m", "dm6.u", "dm6.m", "unmap")
  group_colors <- c("darkgreen", "green2", "orange4", "orange", "grey50")
  
  df_plot <- df %>%
    mutate(sample = factor(sample, levels = rev(sample))) %>%
    select(-total) %>%
    tidyr::gather("group", "count", -1) %>%
    mutate(group = factor(group, levels = rev(groups))) 
  
  p <- df_plot %>%
    ggplot(aes(sample, count, fill = group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = rev(group_colors)) +
    scale_y_continuous(position = "right") +
    # geom_hline(yintercept = mito_mean) +
    xlab(NULL) + ylab("Percentage") +
    coord_flip() + 
    guides(fill = guide_legend(title = NULL,
                               label.position = "top",
                               reverse = TRUE)) +
    theme_bw() +
    theme(
      legend.position = "top"
    )
  
  return(p)
}


barplotCount <- function(df, label = NULL){
  stopifnot(all(c("sample", "count") %in% names(df)))
  
  p <- ggplot(df, aes(sample, count)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(position = "right") +
    coord_flip() +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "grey20"),
      plot.title = element_text(hjust = 0.5)
    )
  
  if(isTRUE(label)) {
    p <- p + 
      geom_text(aes(label = count), hjust = 1.2, color = "grey90")
  }
  
  return(p)
}


##----------------------------------------------------------------------------##
## reads stat
suppressPackageStartupMessages(library(ggcor))

corPlot <- function(x, sample = "all"){
  # df
  df1  <- readr::read_delim(x, "\t", col_names = TRUE, 
                            col_types = readr::cols()) %>%
    dplyr::select(-c(1:3))
  # all samples
  all_samples <- colnames(df1)
  
  # choose sample
  if(sample == "all") {
    df2 <- df1
  } else {
    df2 <- df1 %>%
      dplyr::select(contains(sample))
  }
  # check
  if(ncol(df2) == 0){
    stop("samples not found")
  }
  
  # names
  colnames(df2) <- gsub("ATACseq_DaGal4X|'|.not_MT_trRNA|.map_dm6|_1|\\.1", "", colnames(df2))
  
  # colors
  col <- ggcor:::.default_colors
  col <- rev(col)
  col[6] <- "#F2F2F2"
  
  # library(RColorBrewer)
  # cc = colorRampPalette(brewer.pal(n = 7, name = "RdYlGn"))
  p <- ggcor(df2, type = "lower", fill.colour = col, show.diag = TRUE) + 
    geom_color(data = get_data(type = "lower", show.diag = TRUE)) +
    geom_num(aes(num = r), data = get_data(type = "lower", show.diag = TRUE), 
             colour = "grey90", size = 3) +
    ggtitle(paste0(sample, " correlation (R)"))
  
  return(p)
}

## convert tiff to png
tiff2jpeg <- function(x, outdir = NULL) {
  suppressPackageStartupMessages(library("jpeg"))
  suppressPackageStartupMessages(library("tiff"))
  if(! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  jpeg_name <- gsub(".tiff$", ".jpeg", basename(x))
  jpeg_file <- file.path(outdir, jpeg_name)
  if(! file.exists(jpeg_file)) {
    img <- readTIFF(x, native=TRUE)
    writeJPEG(img, target = jpeg_file, quality = 8)
  }
  return(jpeg_file)
}

##----------------------------------------------------------------------------##
## Fragment length distribution
fragReader <- function(x){
  tmp <- lapply(x, function(f){
    fname <- gsub("frag_length.txt", "", basename(f))
    df <- read.delim(f, header = F, sep = "\t", 
                     col.names = c("length", "count")) %>%
      mutate(sample = fname)
    return(df)
  })
  # merge data.frame
  df <- bind_rows(tmp)
  return(df)
}

fragPlot <- function(df) {
  p <- df %>%
    group_by(sample) %>%
    mutate(frac = count / sum(count)) %>%
    ggplot(aes(length, frac, color = sample)) +
    geom_line(color = "red3", size = .5) + 
    xlab("Fragment length (bp)") + 
    ylab("Normalized read density") +
    scale_x_continuous(breaks = seq(0, 1000, by = 200),
                       labels = seq(0, 1000, by = 200),
                       limits = c(0, 1000), 
                       expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
  return(p)
}


