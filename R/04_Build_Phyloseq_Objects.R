# SETUP ###

# packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")

# shortcuts and functions
ra <- function(x){x/sum(x)}

insertaesedis_function <- function(taxa){
  paste0("_",taxa[grep("Incertae_sedis",taxa)] %>% 
           unique %>% 
           str_remove(".__") %>% 
           str_remove(taxa[grep("Incertae_sedis",taxa)] %>% 
                        unique %>% 
                        str_remove(".__") %>% 
                        str_split("_") %>% 
                        map_chr(1)) %>% 
           str_remove("_") %>% 
           unique())
}

# IMPORT DATA ####

meta <- read_xlsx("./data/Metadata.xlsx") %>% clean_names()
its_seqtab <- readRDS("./output/ITS_seqtab.nochim.clean.RDS")
its_taxa <- readRDS("./output/ITS_taxonomy.RDS")

ssu_seqtab <- readRDS("./output/SSU_seqtab.nochim.clean.RDS")
ssu_taxa <- readRDS("./output/SSU_taxonomy.RDS")

lsu_seqtab <- readRDS("./output/LSU_seqtab.nochim.clean.RDS")
lsu_taxa <- readRDS("./output/LSU_taxonomy.RDS")

# DATA CLEANING ####

# clean taxonomy names
its_insertaesedis_strings <- insertaesedis_function(its_taxa)
ssu_insertaesedis_strings <- insertaesedis_function(ssu_taxa)
lsu_insertaesedis_strings <- insertaesedis_function(lsu_taxa)

its_taxa[,"Phylum"] <- its_taxa[,"Phylum"] %>% str_remove(its_insertaesedis_strings[1]) %>% str_remove(".__")
its_taxa[,"Class"] <- its_taxa[,"Class"] %>% str_remove(its_insertaesedis_strings[2]) %>% str_remove(".__")
its_taxa[,"Order"] <- its_taxa[,"Order"] %>% str_remove(its_insertaesedis_strings[3]) %>% str_remove(".__")
its_taxa[,"Family"] <- its_taxa[,"Family"] %>% str_remove(its_insertaesedis_strings[4]) %>% str_remove(".__")
its_taxa[,"Genus"] <- its_taxa[,"Genus"] %>% str_remove(its_insertaesedis_strings[5]) %>% str_remove(".__")

# clean metadata

# subset to sample ids found in ASV table only
row.names(its_seqtab) <- its_seqtab %>% row.names %>% str_remove("-ITS") %>% str_replace("-","_")
row.names(ssu_seqtab) <- ssu_seqtab %>% row.names %>% str_remove("-18S_FWD_filt.fastq.gz") %>% str_replace("-","_")
row.names(lsu_seqtab) <- lsu_seqtab %>% row.names %>% str_remove("-28S") %>% str_replace("-","_")
meta$sample_id
meta <- 
meta %>% 
  dplyr::filter(sample_id %in% c(row.names(its_seqtab),row.names(ssu_seqtab),row.names(lsu_seqtab))) %>% 
  unique.data.frame()
its_meta <- meta %>% filter(sample_id %in% row.names(its_seqtab))
ssu_meta <- meta %>% filter(sample_id %in% row.names(ssu_seqtab))
lsu_meta <- meta %>% filter(sample_id %in% row.names(lsu_seqtab))



# quick glimpse at possible issues
skimr::skim(meta)

# remove useless column
meta$x16s_analysis <- NULL

# make size fraction a factor (not ordered...yet)
meta$size_fraction_um <- 
  meta$size_fraction_um %>% 
  factor(levels = c("0.3-1","1-53",">53"),ordered = FALSE)

# BUILD PHYLOSEQ ####
tax <- tax_table(its_taxa)
met <- sample_data(its_meta)
sample_names(met) <- its_meta$sample_id
asv <- otu_table(its_seqtab,taxa_are_rows = FALSE)
sample_names(asv) <- its_meta$sample_id
ps <- phyloseq(asv,met,tax)

saveRDS(ps,"./output/ITS_physeq_object.RDS")

tax <- tax_table(ssu_taxa)
met <- sample_data(ssu_meta)
sample_names(met) <- ssu_meta$sample_id
asv <- otu_table(ssu_seqtab,taxa_are_rows = FALSE)
sample_names(asv) <- ssu_meta$sample_id
ps <- phyloseq(asv,met,tax)

saveRDS(ps,"./output/SSU_physeq_object.RDS")

tax <- tax_table(lsu_taxa)
met <- sample_data(lsu_meta)
sample_names(met) <- lsu_meta$sample_id
asv <- otu_table(lsu_seqtab,taxa_are_rows = FALSE)
sample_names(asv) <- lsu_meta$sample_id
ps <- phyloseq(asv,met,tax)

saveRDS(ps,"./output/LSU_physeq_object.RDS")

