# SETUP ####

# You must have cutadapt (v3.5+) and itsxpress (v2.1.0+) installed 
# and in your PATH to run this script

# Load packages
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel)

# Load functions
source("./R/functions.R")

# file paths
path_its <- "./data/raw/2023-0018_ITS"
path_18S <- "./data/raw/StALOHA_18S"
path_28S <- "./data/raw/StALOHA_28S"

# primer sequences
its1f <- "CTTGGTCATTTAGAGGAAGTAA"
its2 <- "GCTGCGTTCTTCATCGATGC"

# TRIM PRIMERS ####

# this function does a prefilter step to remove any reads with N
# those read files can be deleted afterward since the "cutadapt" directory will have the 
# reads used downstream for DADA2

remove_primers(directory = path_its,
               fwd_pattern = "_R1_",
               rev_pattern = "_R2_",
               fwd_primer = its1f,
               rev_primer = its2)
# remove intermediate files if you want
list.files(file.path(path_its,"filtN"),full.names = TRUE) %>% 
  file.remove()

# EXTARCT ITS1 ####
path_its_cutadapt <- file.path(path_its,"cutadapt")

# function creates new fastq files in an ITSx subdirectory for DADA2
run_itsxpress(directory = path_its_cutadapt,
              itsregion = "ITS1",
              taxa_group = "All",
              fwd_pattern = "_R1_")

