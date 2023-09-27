## -----------------------------------------------------------------------------
## Purpose: illustrate ASR output
##
## Author: Katerina Pratsinis
##
## Date Created: 2023-05-17
##
## Email: pratsink@ethz.ch
## -----------------------------------------------------------------------------

rm(list=ls())
# set path as needed
setwd("~/Documents/CBB/2023FS/LS/code/")

################## illustrate ASR trees for selected mutations ##################
# prep file and load libraries
suppressMessages(suppressWarnings({
  library(ggtree)
  library(ggplot2)
  library(ape)
  library(dplyr)
  library(lubridate)
  library(gridExtra)
}))

size= "comp_hf"

gene_dict <- c(
  'Rv0005' = 'gyrB',
  'Rv0006' = 'gyrA',
  'Rv0667' = 'rpoB',
  'Rv0668' = 'rpoC',
  'Rv0682' = 'rpsL',
  'Rv1483' = 'fabG1',
  'Rv1484' = 'inhA',
  'Rv1908c' = 'katG',
  'Rv2043c' = 'pncA',
  'Rv2416c' = 'eis',
  'Rv3794' = 'embA',
  'Rv3795' = 'embB',
  'Rv3919c' = 'gidB'
)

reformat_label <- function(label, lineage) {
  corr <- ifelse(lineage == '2', "pncA V139A", "pncA K96T")
  
  appendage <- ''
  if (endsWith(label, '_R_L')) {
    label <- substr(label, 1, nchar(label)-4)
    appendage <- paste0(", ", corr, " (L)")
  } else if (endsWith(label, '_L')) {
    label <- substr(label, 1, nchar(label)-2)
    appendage <- " (L)"
  } else if (endsWith(label, '_R')) {
    label <- substr(label, 1, nchar(label)-2)
    appendage <- paste0("; ", corr)
  }
  
  if (startsWith(label, 'rrs_')) {
    return(paste0("rrs", " ", strsplit(label, "_")[[1]][2], appendage))
  }
  
  if (startsWith(label, 'rpo')) {
    if (endsWith(label, '104NA_c')) {
      label <- substr(label, 1, nchar(label)-4)
      label <- paste0(label, "0R")
    }
    if (endsWith(label, '125NA_c')) {
      label <- substr(label, 1, nchar(label)-4)
      label <- paste0(label, "2L")
    }
    return(paste0(strsplit(label, "_")[[1]][1], " ", strsplit(label, "_")[[1]][3], " (C)"))
  }
  
  gene_code <- strsplit(label, "_")[[1]][1]
  if (gene_code %in% names(gene_dict)) {
    return(paste0(gene_dict[gene_code], " ", strsplit(label, "_")[[1]][2], appendage))
  }
  
  return(label)  # if no transformation found
}

plot_tree <- function(ds, lineage, mutation){
  #dates imported for tree scale
  dates <- read.delim(paste0("data/1_preprocessing_out/isolation_dates/dates_", ds, "_L", lineage, "_tempest.txt"), header=FALSE)
  tree <- read.nexus(paste0("data/2_ASR_out/pypastml_", ds, "_L", lineage, "_", size, '/named.tree_dated_', ds, '_L', lineage, '.nexus'))
  tab_file <- paste0("data/2_ASR_out/pypastml_", ds, "_L", lineage, "_", size, "/ancestral_states_main.tab")
  tab_df <- read.delim(tab_file, check.names = FALSE, row.names = 1)
  tab_df$node <- rownames(tab_df)
  
  gt <- ggtree(tree, layout = "rectangular", color='steelblue',
               mrsd = max(format(date_decimal(dates[[2]]), "%Y-%m-%d")))
  
  num_tips <- length(tree$tip.label)
  if(num_tips < 500) {
    line_width <- 1
    font_size <- 3
  } else if(num_tips < 1000) {
    line_width <- 0.6
    font_size <- 2
  } else {
    line_width <- 0.3
    font_size <- 1
  }
  
  line_width <- 0.2
  gt$data <- left_join(gt$data, tab_df, by=c("label"="node"))

  # data tree
  gt <- gt + geom_tree(aes(color = ifelse(get(mutation) == 1, "mutant", "wild-type")),
                       size = line_width, show.legend=T) +
    theme_tree2(text=element_text(size=20))
  if(ds=="SA"){
    country = "South Africa lineage "
  }else if(ds=="GA"){
    country = "Georgia lineage "
  }
  mut <- reformat_label(mutation, lineage)
  italic <- strsplit(mut, " ")[[1]][1]
  regular <- strsplit(mut, " ")[[1]][2]
  
  # italic/regular, colors, font size, scale, legend, title
  gt <- gt + scale_color_identity() + 
    ggtitle(bquote(paste(.(country), .(lineage), ", ", italic(.(italic))," ", .(regular)))) +
    theme_tree2(plot.title = element_text(hjust = 0.5, size = 14),) +
    scale_x_continuous(name="Year")
  
  gt <- gt + scale_color_manual(values = c("wild-type" = "steelblue", "mutant" = "red")) + 
    ggtitle(bquote(paste(.(country), .(lineage), ", ", italic(.(italic)), " ", .(regular)))) +
    theme_tree2(
      plot.title = element_text(hjust = 0.5, size = 14)
    )+ 
    scale_x_continuous(name="Year") +
    theme(
      legend.position = c(0.1, 0.95),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size=12),
      legend.key.size = unit(1.5, "cm"),
      legend.key.height = unit(0.5, "cm")
      
    )
  
  return(gt)
}

plots <- list(
  plot_tree("SA", "2", "Rv0667_S450L"),
  plot_tree("SA", "4", "Rv3795_Q497R"),
  plot_tree("GA", "2", "Rv1483_C15T"),
  plot_tree("GA", "4", "Rv1908c_S315T")
)

combined_plot <- do.call(grid.arrange, c(plots, ncol=2, top=NULL, bottom=NULL))

ggsave(filename="data/2_ASR_out/ASR-trees.png", plot=combined_plot, width=10, height=13)