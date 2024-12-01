# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
set.seed(666)

# Load seq tables
ITS_seqtab <- readRDS("./output/ITS_seqtab.nochim.clean.RDS")
SSU_seqtab <- readRDS("./output/SSU_seqtab.nochim.clean.RDS")
LSU_seqtab <- readRDS("./output/LSU_seqtab.nochim.clean.RDS")

# Database paths
its_db <- "./taxonomy/Eukaryome_General_ITS_v1.8_reformatted.fasta.gz"
ssu_db <- "./taxonomy/Eukaryome_General_SSU_v1.8_reformatted.fasta.gz"
lsu_db <- "./taxonomy/Eukaryome_General_LSU_v1.8_reformatted.fasta.gz"

# ASSIGN TAXONOMY ####

its_taxa <- assignTaxonomy(ITS_seqtab, its_db, multithread=(parallel::detectCores()-1))
ssu_taxa <- assignTaxonomy(SSU_seqtab, ssu_db, multithread=(parallel::detectCores()-1))
lsu_taxa <- assignTaxonomy(LSU_seqtab, lsu_db, multithread=(parallel::detectCores()-1))

# Save taxonomy file
saveRDS(its_taxa, file = "./output/ITS_taxonomy.RDS")
saveRDS(ssu_taxa, file = "./output/SSU_taxonomy.RDS")
saveRDS(lsu_taxa, file = "./output/LSU_taxonomy.RDS")
