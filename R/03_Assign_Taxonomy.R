# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")


# Load seq table
seqtab <- readRDS("./output/seqtab.nochim.clean.RDS")

# Database paths
its_db <- "./taxonomy/sh_general_release_dynamic_all_04.04.2024.fasta.gz"

# ASSIGN TAXONOMY ####
# using UNITE "All Eukaryotes" version 10.0
# https://dx.doi.org/10.15156/BIO/2959334
taxa <- assignTaxonomy(seqtab, its_db, multithread=(parallel::detectCores()-1))

# Save taxonomy file
saveRDS(taxa, file = "./output/ITS_taxonomy.RDS")
