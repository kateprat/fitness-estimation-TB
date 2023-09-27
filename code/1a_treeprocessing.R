## -----------------------------------------------------------------------------
## Purpose: clean and date ML tree, prep for ASR
##
## Author: Katerina Pratsinis
##
## Date Created: 2023-01-27
##
## Email: pratsink@ethz.ch
## -----------------------------------------------------------------------------

rm(list=ls())
# set path as needed
setwd("~/Documents/CBB/2023FS/LS/code/")


################################## PREP SCRIPT ##################################
# clean environment and load libraries -----------------------------------------
list.of.packages <- c("tictoc", "phylotools", "data.table", "janitor",
                      "tidyverse", "ape", "Rlsd2", "ggtree", "tidytree", "lubridate", 
                      "TreeTools", "reshape2", "colorspace", "corrplot")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load necessary packages
suppressPackageStartupMessages({
  library(tictoc)
  library(phylotools)
  library(data.table)
  library(janitor)
  library(tidyverse)  # includes dplyr, tidyr, tibble, readr, ggplot2
  library(ape)
  library(Rlsd2)
  library(ggtree)
  library(tidytree)
  library(lubridate)
  library(TreeTools)
  library(reshape2)
  library(colorspace)
  library(corrplot)
})

# vars -------------------------------------------------------------------------
tic("Runtime")
start <- Sys.time()

ds = "SA" # "SA" or "GA"
lineage = 2
write_files <- T
save_plots <- F
outgroup <- "Mycobacterium_canettii"

# prep import
if(ds == "SA"){
  import=c("SouthAfrica","variable", "SA", "variable_positions")
}else if(ds=="GA"){
  import=c("Georgia","_polymorphic", "Georgia", "_var_aln")
}

# import files -----------------------------------------------------------------
# import to extract sequence length for lsd
lin_woDR <-
  suppressWarnings(phylotools::read.fasta(paste0("data/", import[1], "_L", lineage,
                                                 "_wo_DR_outgr_", import[4], ".fasta")))

# import metadata to extract isolation dates
full_metadata <- read.delim(paste0("data/", import[3], "_metadata.txt"))
names(full_metadata) <- tolower(names(full_metadata))
full_metadata <- full_metadata[full_metadata$g_number %in% lin_woDR$seq.name,
                               c("g_number", "isolation_date")]
rownames(full_metadata) <- 1:nrow(full_metadata)

# to create plain tree for lsd
tree <- read.tree(paste0("data/0_tree-inference_out/iqtree_", ds, "_L", lineage,
                         "/", import[1], "_L", lineage, "_wo_DR_outgr_",
                         import[4], ".fasta.treefile"))
(ggtree(tree) +  geom_text(aes(label=label),
                           data=td_filter(isTip & node %in% nodeid(tree, "Mycobacterium_canettii")),
                           hjust=0) + ggtitle(paste0("IQ-TREE for Lineage ", lineage, " (", ds, ")")))

boottrees <- read.tree(paste0("data/0_tree-inference_out/iqtree_", ds, "_L",
                              lineage, "/", import[1], "_L", lineage,
                              "_wo_DR_outgr_", import[4], ".fasta.boottrees"))

#################################### SCRIPT ####################################
# extract isolation date from metadata -----------------------------------------
if(ds=="SA"){
  full_dates <- full_metadata %>% mutate(isolation_date = 2008+ isolation_date) # in decimal format
}else if(ds=="GA"){
  # Convert isolation_date yy.mm.dd to yyyy-mm-dd
  full_dates <- full_metadata %>% mutate(isolation_date = decimal_date(dmy(full_metadata$isolation_date)))
}


if(write_files){
  write.table(full_dates,
              file=paste("data/1_preprocessing_out/isolation_dates/dates_", ds,
                         "_L", lineage, "_tempest.txt", sep=""),
              quote = F, sep = "\t", row.names = F, col.names = F)
  cat(length(full_dates[,2]), "\n",
      file = paste("data/1_preprocessing_out/isolation_dates/dates_", ds, "_L",
                   lineage, "_LSD.txt", sep=""))
  write.table(full_dates,
              file=paste("data/1_preprocessing_out/isolation_dates/dates_", ds,
                         "_L", lineage, "_LSD.txt", sep=""),
              append = T,
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# tree handling: drop tip, date, resolve zero branch lengths -------------------
clockrate <- c(NA,5*10^-8,NA,5*10^-8)

rm_multif <- function(tree){
  tree.rm.multif <- multi2di(tree)
  tree.rm.multif$edge.length <- pmax(tree.rm.multif$edge.length,1/365)
  return(tree.rm.multif)
}

plot_tree <- function(tree, title){
  t <- ggtree(tree) + ggtitle(title)
  plot(t)
}

write_tree_file <- function(tree, filename){
  if(write_files){
    write.tree(tree, file=filename)
  }
}

tree_handling <- function(tree, main){
  # drop og tip used for phylotree inference in iqtree
  plain <- drop.tip(tree, which(tree$tip.label==outgroup))
  if(main==T){
    plot_tree(plain, paste0("IQ-TREE: Lineage ", lineage, " (", ds, ")"))
    write_tree_file(plain, paste0("data/1_preprocessing_out/plain_", ds, "_L",
                                  lineage, ".tre"))
  }
  # LSD
  dated_tree <- lsd2(inputTree = plain,
                     inputDate = paste0("data/1_preprocessing_out/isolation_dates/dates_",
                                        ds, "_L", lineage, "_LSD.txt"),
                     nullblen = 0, seqLen = nchar(lin_woDR[1,2]),
                     outFile = paste0("data/1_preprocessing_out/LSD_out/LSD_", ds, "_L", lineage),
                     givenRate = clockrate[lineage])
  dated_tree2 <- dated_tree$dateNexusTreeFile@phylo
  if(length(any(dated_tree2$edge.length == 0))>0){
    dated_tree2 <- rm_multif(dated_tree2)
  }
  if(main==T){
    t <- ggtree(dated_tree2) +
      plot_tree(dated_tree2, paste("Dated tree: Lineage ", lineage, " (", ds, ")", sep=""))
    write_tree_file(dated_tree2, paste0("data/1_preprocessing_out/dated_trees/dated_", ds, "_L", lineage, ".nwk"))
  }
  return(dated_tree2)
}

dated_tree <- tree_handling(tree = tree, main = T)

for(i in 1:length(boottrees)){
  dated_bt <- tree_handling(tree = boottrees[[i]], main = F)
  write_tree_file(dated_bt, paste0("data/1_preprocessing_out/dated_trees/dated_", ds, "_L", lineage, "_bt", i, ".tre"))
}
