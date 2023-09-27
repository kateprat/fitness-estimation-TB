## -----------------------------------------------------------------------------
## Purpose: illustrate fitness trees & significant mutations
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

suppressMessages(suppressWarnings(library(tictoc)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(phylotools)))
suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(ggtree)))
suppressMessages(suppressWarnings(library(stringr)))

# functions
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


read_data_files <- function(ds, lineage, size) {
  tab_file <- paste0("data/2_ASR_out/pypastml_", ds, "_L", lineage, "_", size, "/ancestral_states_main.tab")
  fit_file <- paste0("data/3_phyloTF_out/phyloTF_", ds, "_L", lineage, "_", size, "_a.csv")
  tree_file <- paste0("data/2_ASR_out/pypastml_", ds, "_L", lineage, "_", size, '/named.tree_dated_', ds, '_L', lineage, '.nwk')
  
  if (!(file.exists(tab_file) && file.exists(fit_file) && file.exists(tree_file))) {
    return(NULL)
  }
  
  list(
    tab_df = read.delim(tab_file, check.names = FALSE, row.names = 1),
    fit_df = read.csv(fit_file, check.names = FALSE, row.names = 1),
    tree = read.tree(tree_file)
  )
}

log_transform_branches <- function(tree, base = exp(1)) {
  tree$edge.length <- log(tree$edge.length + 1, base)
  return(tree)
}

calculate_node_fitness <- function(tab_df, fit_df) {
  node_fitness <- rep(0, nrow(tab_df))
  names(node_fitness) <- rownames(tab_df)
  node_fitness[which(rowSums(tab_df) == 0)] <- 1
  
  for (node in names(which(node_fitness != 1))) {
    mutated_sites <- names(tab_df[node,])[which(tab_df[node,] == 1)]
    fitness_values <- fit_df['main', mutated_sites]
    node_fitness[[node]] <- prod(fitness_values, na.rm = TRUE)
  }
  
  data.frame(node = names(node_fitness), fitness = node_fitness, stringsAsFactors = FALSE)
}

create_column <- function(p, tip_data, col_name, x_start, width, lineage) {
  p <- p + geom_rect(aes(xmin=x_start, xmax=x_start + width, ymin=-Inf, ymax=Inf), fill="grey", inherit.aes=FALSE)
  p <- p + geom_segment(data=filter(tip_data, !!sym(col_name) == 1),
                        aes(x=x_start, xend=x_start + width, y=y, yend=y),
                        color="blue3", inherit.aes=FALSE)
  
  display_label <- reformat_label(col_name, lineage)
  
  # split label
  split_label <- strsplit(display_label, " ")[[1]]
  italic_part <- split_label[1]
  regular_part <- if(length(split_label) > 1) paste(split_label[-1], collapse = " ") else ""
  
  top_y <- max(p$data$y, na.rm = TRUE)
  label_x <- (x_start + x_start + width) / 2
  
  p <- p + geom_text(aes(x = label_x, y = top_y, label = display_label),
                     angle = 90, vjust = 0.5, hjust = 1,
                     inherit.aes = FALSE, color = "white")

  # p <- p + geom_text(aes(x = label_x, y = top_y, label = italic_part),
  #                    angle = 90, vjust = 0.5, hjust = 3, 
  #                    inherit.aes = FALSE, color = "white", fontface = "italic")
  # 
  # p <- p + geom_text(aes(x = label_x, y = top_y, label = regular_part), 
  #                    angle = 90, vjust = 0.5, hjust = 1, 
  #                    inherit.aes = FALSE, color = "white")
  
  return(p)
}

fitness_tree <- function(t, tab_df, fit_df, ds, lineage, size, mutations, node ="root", xlim_val=c(0,50)) {
  t <- extract.clade(t, node)
  gt <- ggtree(t)
  gt$data <- left_join(gt$data, fit_df, by = c("label" = "node"))
  
  p <- ggtree(gt$data, aes(color=fitness)) + scale_color_viridis_c(name="Fitness")
  p <- p + theme(legend.position=c(0.05,0.9))
  # p <- p + ggtitle(paste0(ds, " Lineage ", lineage, " (", size, ")"))
  # p <- ggtree(gt$data, aes(color=fitness)) + scale_color_viridis_c(name="Fitness") # print node label
  # p <- p + geom_text(aes(label=label), vjust=1.5, hjust=0.5, size=3, color="black") #print node label
  
  tip_data <- gt$data[gt$data$isTip == TRUE, ]
  tab_df$label <- rownames(tab_df)
  tip_data <- left_join(tip_data, tab_df)
  
  x_start <- max(p$data$x) + 1
  width <- 10
  for (mut in mutations) {
    p <- create_column(p, tip_data, mut, x_start, width, lineage)
    x_start <- x_start + width + 1
  }
  return(p)
}

# main
tic("Runtime")

trees <- list()
plots <- list()

combinations <- expand.grid(ds=c('GA'), lineage=c('4'), size=c('comp_hf')) #choose SA or GA; 2 or 4; comp_hf or ten

for (i in 1:nrow(combinations)) {
  comb <- combinations[i,]
  ds <- comb$ds
  lineage <- comb$lineage
  size <- comb$size
  
  data <- read_data_files(ds, lineage, size)
  if (is.null(data)) next
  data$tree <- log_transform_branches(data$tree)

  if(ds=='SA' && lineage == '2'){
    mutations <- c("Rv0667_S450L", "Rv0667_L430P", "Rv0667_H445Y", "Rv0667_H445R", 
                   "Rv0667_L452P", "Rv0667_H445D", "Rv0667_H445L",
                   "Rv0667_D435Y", "Rv0667_D435V", "rpoC_Rv0668_D485Y",
                   "rpoC_Rv0668_V483G", "Rv0682_K88R", "Rv1483_C15T",
                   "Rv1484_I194T_L", "Rv2043c_D8N_L", "Rv2043c_Y103_L",
                   "Rv3795_M306V_L", "Rv3795_M306I_L", "Rv3919c_L79S_L")
  }else if(ds=='SA' && lineage == '4'){
    mutations <- c("Rv0006_A90V_L", "Rv0667_S450W", "Rv0667_S450L", "Rv0667_S441L",
                   "Rv0667_L430P", "Rv0667_L452P", "Rv0667_H445N", "Rv0667_H445D",
                   "Rv0667_D435Y", "Rv0667_D435V", "Rv0667_H445Y",
                   "rpoC_Rv0668_V125NA_c", "Rv0682_K43R_L", "Rv0682_K88Q_R_L",
                   "Rv1483_C15T", "Rv1908c_S315T", "Rv1908c_S315R_L", 
                   "Rv2043c_G97C_L", "Rv2043c_Q10P_L", "Rv2043c_H71Y_L",
                   "Rv3794_C16G_L", "Rv3795_G406A_L", "Rv3795_G406D_L",
                   "Rv3795_M306V_L", "Rv3795_Q497R")
  }else if(ds=='GA' && lineage == '2'){
    mutations <- c("Rv0006_D94A", "Rv0006_D94G_L", "Rv0006_D94Y_L", 
                   "rpoB_Rv0667_L731P_c", "Rv0667_L452P", "Rv0667_L430P", "rpoC_Rv0668_G332C_c",
                   "rpoC_Rv0668_D485N_c", "rpoC_Rv0668_V483G_c",
                   "Rv0682_K88R_L", "Rv1483_C15T", "Rv1908c_S315T",
                   "Rv2043c_T142M_L", "Rv2416c_G10A_L", "Rv3794_C12T_L",
                   "Rv3795_M306I_L", "Rv3795_G406S_L", "Rv3795_G406D",
                   "Rv3795_G406A", "Rv3795_D354A_L")
  }else if(ds=='GA' && lineage == '4'){
    mutations <-  c("Rv0667_D435V_L", "Rv0667_D435Y", "Rv0667_H445D", "Rv0667_H445Y",
                    "Rv0667_S450L", "Rv0682_K43R", "Rv0682_K88R_L", "Rv1483_C15T",
                    "Rv1908c_S315T", "Rv2416c_C12T_L", "Rv2416c_C14T_L", "Rv3794_C12T_L",
                    "Rv3795_D1024N", "Rv3795_G406D", "Rv3795_M306I", "Rv3795_M306V")
  }
  
  fit_df <- calculate_node_fitness(data$tab_df, data$fit_df)
  plot <- fitness_tree(data$tree, data$tab_df, fit_df, ds, lineage, size, mutations, xlim_val = c(0, 50)) #node="n91",
  print(plot)
  plots[[paste0(ds, "_L", lineage, "_", size)]] <- plot
  trees[[paste0(ds, "_L", lineage, "_", size)]] <- data$tree
}

saveRDS(plots, "data/6_fitness-trees_out/fitness_trees.rds")
toc()

