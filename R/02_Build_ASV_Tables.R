# SETUP ####

# Packages
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(parallel)

# Options
set.seed(666)


# File paths

path <- "./data/raw/2023-0018_ITS/ITSx" 
filtpath <- "./data/clean/ITS" 
if(!file_test("-d", filtpath)) dir.create(filtpath)

fns <- sort(list.files(path, full.names = TRUE, pattern = "ITS.fastq.gz$")) 

filts_f <- file.path(filtpath, paste0(sample.names, "_FWD_filt.fastq.gz"))

# sample names
sample.names <- basename(fns) %>% strsplit("_") %>% map_chr(1)

# inspect read quality
plotQualityProfile(fns[1]) + ggtitle("Example forward reads")


# FILTER AND TRIM ####
out <- filterAndTrim(fns, filts_f, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=2, # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=(parallel::detectCores()-1)) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./output/trackreads.RDS")


# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered")
p3 / p4


# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2) 
saveRDS(errF,"./output/errF.RDS")

plotErrors(errF, nominalQ=FALSE)


# dereplication
derepF <- derepFastq(filts_f, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepF) <- 
  names(derepF) %>% 
  str_remove("_FWD_filt.fastq.gz")

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
saveRDS(dadaFs,"output/dadaFs.RDS")


# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF","nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./output/ITS_read_counts_at_each_step.csv", row.names = TRUE)


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./output/seqtab.nochim.clean.RDS")
