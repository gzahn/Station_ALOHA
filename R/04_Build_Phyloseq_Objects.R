# SETUP ###

# packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")

# shortcuts
ra <- function(x){x/sum(x)}

# IMPORT DATA ####
meta <- read_xlsx("./data/Metadata.xlsx") %>% clean_names()
seqtab <- readRDS("./output/seqtab.nochim.clean.RDS")
taxa <- readRDS("./output/ITS_taxonomy.RDS")

# DATA CLEANING ####

# clean taxonomy names
insertaesedis_strings <- 
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

taxa[,"Phylum"] <- taxa[,"Phylum"] %>% str_remove(insertaesedis_strings[1]) %>% str_remove(".__")
taxa[,"Class"] <- taxa[,"Class"] %>% str_remove(insertaesedis_strings[2]) %>% str_remove(".__")
taxa[,"Order"] <- taxa[,"Order"] %>% str_remove(insertaesedis_strings[3]) %>% str_remove(".__")
taxa[,"Family"] <- taxa[,"Family"] %>% str_remove(insertaesedis_strings[4]) %>% str_remove(".__")
taxa[,"Genus"] <- taxa[,"Genus"] %>% str_remove(insertaesedis_strings[5]) %>% str_remove(".__")

# clean metadata

# subset to sample ids found in ASV table only
row.names(seqtab) <- seqtab %>% row.names %>% str_remove("-ITS") %>% str_replace("-","_")

meta <- 
meta %>% 
  dplyr::filter(sample_id %in% row.names(seqtab)) %>% 
  unique.data.frame()

# quick glimpse at possible issues
skimr::skim(meta)

# remove useless column
meta$x16s_analysis <- NULL

# make size fraction an ordered factor
meta$size_fraction_um <- 
  meta$size_fraction_um %>% 
  factor(levels = c("0.3-1","1-53",">53"),ordered = TRUE)

# BUILD PHYLOSEQ ####
tax <- tax_table(taxa)
met <- sample_data(meta)
sample_names(met) <- meta$sample_id
asv <- otu_table(seqtab,taxa_are_rows = FALSE)

ps <- phyloseq(asv,met,tax)

# quick plot of kingdoms
ps %>% 
  transform_sample_counts(ra) %>% 
  plot_bar(fill="Kingdom")

# CLEAN PHYLOSEQ ####

# do we want to dump non-fungal stuff?


# EXPORT FINAL PHYLOSEQ ####
saveRDS(ps,"./output/clean_phyloseq_object.RDS")


