## -----------------------------------------------------------------------------
## Purpose: prepare binary table for ASR
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
list.of.packages <- c("tictoc", "data.table", "tidyverse", "tidytree",
                      "corrplot")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(suppressWarnings({
  library(tictoc)
  library(data.table)
  library(tidyverse)
  library(tidytree)
  library(corrplot)
}))

# vars -------------------------------------------------------------------------
tic("Runtime")
start <- Sys.time()

ds = "SA" # "SA" or "GA"
lineage = 2 # 2 or 4
write_files <- T

threshold <- 0.005
# prep import
if(ds == "SA"){
  if(lineage == 4){
    threshold <- 0.01
  }
  import=c("SouthAfrica","variable", "SA", "variable_positions")
}else if(ds=="GA"){
  import=c("Georgia","_polymorphic", "Georgia", "_var_aln")
}

# helper functions -------------------------------------------------------------
amino_acid_mapping <- c(
  "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D", 
  "Cys" = "C", "Gln" = "Q", "Glu" = "E", "Gly" = "G", 
  "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K", 
  "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S", 
  "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V"
)

process_mutations <- function(mutation_df, lin_copy, lin_all) {
  site_count_before <- ncol(lin_copy)
  count <- 0
  for(row in rownames(mutation_df)){
    count <- count+1
    gen_pos <- as.character(as.integer(mutation_df[row,"GENOMIC_POSITION"]))
    ref_base <- mutation_df[row, "REF_BASE"]
    alt_base <- mutation_df[row, "ALT_BASE"]
    col_name <- paste0(mutation_df[row, "LOCUS"], "_", mutation_df[row, "SUBSTITUTION"])
    
    # Check if gen_pos exists in the column names of lin_copy
    if(nchar(ref_base)==1){
      cat(count, " ", gen_pos, " mutates alone from ", ref_base, " to ", alt_base, "\n")
      check_mut <- lin_copy[[gen_pos]]
      print(length(check_mut))
    }else if(nchar(ref_base)==2){
      next_col <- as.character(as.integer(gen_pos)+1)
      if(next_col %in%colnames(lin_all)){
        cat(count," ", gen_pos, " and ", next_col, " are both in mutation from ",
            ref_base, " to ", alt_base, "\n")
        check_mut <- paste(lin_copy[,gen_pos], lin_copy[,next_col], sep="")
      }else{
        cat(count, " ", gen_pos, " does not have its mutation neighbour ", next_col, "\n")
        next
      }
    }else{#manually checked all 4 data sets, none of them have 3 or more consecutive mutated sites from DR list
      cat(count, " ", gen_pos, " mutates with more than 1 neighbor\n", sep="")
      next
    }
    binary_check <- as.integer(check_mut == alt_base)
    binary_sum <- sum(binary_check)
    cat("This mutation happens ", binary_sum, " times\n")
    if(binary_sum>0){
      if(col_name %in% colnames(lin_copy)){
        lin_copy[,col_name] <- as.integer((lin_copy[,col_name] | binary_check) > 0)
      }else{
        lin_copy[,col_name] <- binary_check
      }
    }
  }
  binary_mutations <- lin_copy[,(site_count_before+1):ncol(lin_copy)]
  lin_copy <- lin_copy[,1:site_count_before]
  binary_mutations_hf <- binary_mutations[,-which(colSums(binary_mutations)/nrow(binary_mutations)<threshold)]
  
  return(list(binary_mutations = binary_mutations, 
              binary_mutations_hf = binary_mutations_hf))
}

linked_sites <- function(binary, text, min_linkage = 0.95) {
  total_mutations <- colSums(binary)
  for(site in colnames(binary)){
    mutated_rows <- which(binary[, site] == 1)
    sub_df <- binary[mutated_rows, ]
    sub_mutations <- colSums(sub_df)
    linked_sites <- names(sub_mutations)[sub_mutations / total_mutations[site] >= min_linkage]
    if(length(linked_sites)>1){
      binary_sub <- binary[, linked_sites]
      binary[,site] <- apply(binary_sub, 1, prod)
      names(binary)[which(names(binary)==site)] <- paste0(site, "_L")
      linked_sites <- setdiff(linked_sites, site)
      text <- paste0(text, site, ": ", paste(linked_sites, collapse = ", "), "\n")
    }
  }
  text <- paste0(text, "\n")
  return(list("Data" = binary, "Text" = text))
}

comp_sites <- function(comp, binary, text, min_linkage=1){
  muts <- colSums(comp)
  for(c_site in colnames(comp)){
    if(muts[c_site]<(nrow(binary))/100) next
    sub_df <- binary[which(comp[,c_site]==1),]
    sub_muts <- colSums(sub_df)
    linked_sites <- names(sub_muts)[sub_muts / muts[c_site] >= min_linkage]
    linked_sites <- setdiff(linked_sites, c_site)
    if(length(linked_sites)>0){
      text <- paste0(text, c_site, ": ", paste(linked_sites, collapse=", "), "\n")
      binary$new <- comp[,c_site]
      names(binary)[ncol(binary)] <- c_site
    }
  }
  text <- paste0(text, "\n")
  return(list("Data" = binary, "Text" = text))
}

# import files -----------------------------------------------------------------
df <- read.delim(paste0("data/", import[1], "_L", lineage, "_outgr_",
                        import[2], "_positions.txt"), row.names = 1)
mut_DR <- read.delim("data/mutations/Known_DR_mutations_TBRU_and_WHO_combined.202206.tsv",
                     row.names = 1)
comp_muts <- read.delim("data/mutations/Compensatory_mutations.tsv",
                        row.names=NULL)

# transform files --------------------------------------------------------------
df <- df[,-1]
lin_all <- as.data.frame(t(df))
rm(df)

# genpos associated with DR throughout TB
mut_DR <- mut_DR[mut_DR$GENOMIC_POSITION %in% colnames(lin_all),]
mut_DR <- mut_DR[order(mut_DR$GENOMIC_POSITION),]
mut_DR <- unique(mut_DR[mut_DR$ASSOCIATED_WITH_RESISTANCE=="YES",
                        c("GENOMIC_POSITION", "REF_BASE", "ALT_BASE", "LOCUS",
                          "SUBSTITUTION", "DRUG")])
rownames(mut_DR) <- 1:nrow(mut_DR)
# make sure there are no white spaces in REF_BASE and ALT_BASE
mut_DR$GENOMIC_POSITION <- trimws(mut_DR$GENOMIC_POSITION)
mut_DR$REF_BASE <- trimws(mut_DR$REF_BASE)
mut_DR$ALT_BASE <- trimws(mut_DR$ALT_BASE)

# make sure there are no white spaces in REF_BASE and ALT_BASE
comp_muts$GENOMIC_POSITION <- as.integer(trimws(comp_muts$Genomic.position.s.))
comp_muts$REF_BASE <- trimws(comp_muts$Ref..Allele.s.)
comp_muts$ALT_BASE <- trimws(comp_muts$Alt..Allele.s.)
comp_muts <- rename(comp_muts, "LOCUS" = "Gene")
comp_muts <- rename(comp_muts, "DRUG"="Drug")
comp_muts$SUBSTITUTION <- sapply(comp_muts$`AA.Substitution`, function(x) {
  parts <- regmatches(x, gregexpr("[[:alpha:]]{3}|\\d{1,4}", x))[[1]]
  parts[c(1, length(parts))] <- amino_acid_mapping[parts[c(1, length(parts))]]
  paste0(parts[c(1, length(parts)-1, length(parts))], collapse="")
})
comp_muts <- comp_muts[,colnames(mut_DR)]
comp_muts <- comp_muts[comp_muts$GENOMIC_POSITION %in% colnames(lin_all),]
comp_muts <- comp_muts[!comp_muts$GENOMIC_POSITION %in% mut_DR$GENOMIC_POSITION, ]
comp_muts <- comp_muts[order(comp_muts$GENOMIC_POSITION),]
comp_muts$GENOMIC_POSITION <- as.character(comp_muts$GENOMIC_POSITION)
rownames(comp_muts) <- 1:nrow(comp_muts)

# unique should be redundant since duplicates already removed
keep <- unique(append(unique(comp_muts$GENOMIC_POSITION),
                      unique(mut_DR$GENOMIC_POSITION)))
lin_all <- lin_all[,colnames(lin_all) %in% keep]

#################################### SCRIPT ####################################
# prep binary file for ASR in PastML -------------------------------------------
# includes Pearson correlation computation -------------------------------------
lin_copy <- lin_all[2:nrow(lin_all),]

# Process mut_DR
output <- process_mutations(mut_DR, lin_copy, lin_all)
binary_DR <- output$binary_mutations
binary_hf <- output$binary_mutations_hf

# check for correlation among DR sites------------------------------------------
cor_matrix <- cor(binary_hf)
corrplot(cor_matrix, method = "color", type = "lower", order = "hclust",
         #tl.pos = "n", # remove variable labels
         # title = paste("Correlation heatmap of lineage", lineage,
         #               "genomic sites"), mar=c(0,0,1,0)
)

high_corr <- which(cor_matrix >= 0.95 & cor_matrix <= 1 & lower.tri(cor_matrix),
                   arr.ind = TRUE)
high_corr_names <- paste0(colnames(binary_hf)[high_corr[,2]], " & ",
                          colnames(binary_hf)[high_corr[,1]])

# initiate text to summarize all combined sites --------------------------------
txt <- paste0("The following contains all correlations, linkages and compensatory",
              " mutations for the ", ds, " dataset, lineage ", lineage, ":\n\n")

# print genomic sites with high correlation-------------------------------------
cat("Genomic Sites with Correlation >= 0.95 in lineage", lineage, ":\n")
if(length(high_corr)>0){
  txt <- paste0(txt, "Correlations:\n", paste(high_corr_names, collapse = "\n"))
  cat(paste(high_corr_names, collapse = "\n"))
  names(binary_hf)[names(binary_hf) %in% colnames(binary_hf)[high_corr[,2]]] <-
    paste0(colnames(binary_hf)[high_corr[,2]],"_R")
  binary_hf <- binary_hf[,!colnames(binary_hf)%in%colnames(binary_hf)[high_corr[,1]]]
  txt <- paste0(txt, "\n\n")
}

################################## find links ##################################
txt <- paste0(txt, "Linked sites:\n")
DataText <- linked_sites(binary_hf, txt) # comment out if you do not want linkages bin ordinary binary_hf instead
binary_hf <- DataText$Data
txt <- DataText$Text

# save binaries & freqs --------------------------------------------------------
freq_hf <- colSums(binary_hf)

sorted_hf <- freq_hf[order(unlist(freq_hf), decreasing = TRUE)]
high_mutation_sites <- c()  # Initialize vector to keep track of selected sites
num_sites <- 0  # Initialize count of selected sites
sorted_hf <- freq_hf[order(unlist(freq_hf), decreasing = TRUE)]
high_mutation_sites <- c()  # Initialize vector to keep track of selected sites
num_sites <- 0  # Initialize count of selected sites

for (site in names(sorted_hf)) {
  if (num_sites == 0) {
    high_mutation_sites <- c(high_mutation_sites, site)
    num_sites <- num_sites + 1
  } else {
    linkage_vals <- sapply(high_mutation_sites, function(prev_site) {
      sum(binary_hf[, prev_site] & binary_hf[, site]) / sum(binary_hf[, prev_site])
    })
    if (all(linkage_vals < 0.8)) {
      high_mutation_sites <- c(high_mutation_sites, site)
      num_sites <- num_sites + 1
    }
  }
  if (num_sites == 10) {
    break
  }
}

binary_ten <- binary_hf[,colnames(binary_hf) %in% high_mutation_sites]

mut_freq <- data.frame("mut_count" = colSums(binary_hf),
                       "mut_freq" = round(colSums(binary_hf)/nrow(binary_hf),4))
mut_freq_all <- data.frame("mut_count" = colSums(binary_DR),
                           "mut_freq" = round(colSums(binary_DR)/nrow(binary_DR),4))

# expand binary files with interactions and linkages incl. comp muts -----------
output <- process_mutations(comp_muts, lin_copy, lin_all)
binary_comp <- output$binary_mutations
#binary_comp_hf <- output$binary_mutations_hf

txt <- paste0(txt, "Compensatory mutations occuring at high frequency and linked to DR mutations: \n")
DataText <- comp_sites(binary_comp, binary_hf, txt, min_linkage= 0.95)
binary_comp_hf <- DataText$Data
txt <- DataText$Text

comp_freq <- data.frame("mut_count" = colSums(binary_comp_hf),
                        "mut_freq" = round(colSums(binary_comp_hf)/nrow(binary_comp_hf),4))
comp_freq_all <- data.frame("mut_count" = colSums(binary_comp),
                            "mut_freq" = round(colSums(binary_comp)/nrow(binary_comp),4))

txt <- paste0(txt, "# genomes: ", nrow(binary_DR), "\n", "# DR sites: ",
              ncol(binary_DR), "\n", "# hf DR sites: ", ncol(binary_hf), "\n",
              "# comp sites linked to a DR site: ",
              ncol(binary_comp), "\n",
              "# comp sites linked to a hf DR site with hf: ",
              ncol(binary_comp_hf)-ncol(binary_hf), "\n")


if(write_files){
  write.csv(mut_freq, file=paste0("data/1_preprocessing_out/", ds, "_L",
                                  lineage, "_mut-freq.csv"))
  write.csv(mut_freq_all, file=paste0("data/1_preprocessing_out/", ds, "_L",
                                      lineage, "_mut-freq_all.csv"))
  write.csv(comp_freq, file=paste0("data/1_preprocessing_out/", ds, "_L",
                                   lineage, "_comp-freq.csv"))
  write.csv(comp_freq_all, file=paste0("data/1_preprocessing_out/", ds, "_L",
                                       lineage, "_comp-freq_all.csv"))
  write.csv(binary_ten, file=paste0("data/1_preprocessing_out/binary/binary_",
                                    ds, "_L", lineage, "_ten.csv"))
  write.csv(binary_hf, file=paste0("data/1_preprocessing_out/binary/binary_",
                                   ds, "_L", lineage, "_hf.csv"))
  write.csv(binary_comp_hf, file=paste0("data/1_preprocessing_out/binary/binary_",
                                        ds, "_L", lineage, "_comp_hf.csv"))
  write.table(txt, paste0("data/1_preprocessing_out/RLc_", ds, "_L", lineage,".txt"),
              row.names = F, col.names = F, quote = F)
}

toc()