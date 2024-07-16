
# remove_primers() ####
# Function to remove primers from raw amplicon files

remove_primers <- function(directory, # where raw sequences live
                           fwd_pattern="_R1_",
                           rev_pattern="_R2_",
                           fwd_primer="GTGCCAGCMGCCGCGGTAA",
                           rev_primer="GGACTACHVGGGTWTCTAAT"){
  
  library(tidyverse); packageVersion("tidyverse")
  library(dada2); packageVersion("dada2")
  library(purrr); packageVersion("purrr")
  library(Biostrings); packageVersion("Biostrings")
  library(ShortRead); packageVersion("ShortRead")
  library(parallels)
  
  
  # File parsing
  path <- directory # CHANGE to the subdirectory containing your demultiplexed fastq files
  
  # your filenames might have a different pattern for determining FWD and REV reads
  # Change "pattern" to accomodate your file names
  fnFs <- sort(list.files(path, pattern = fwd_pattern, full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = rev_pattern, full.names = TRUE))

  FWD <- fwd_primer # Sequence of FWD primer
  REV <- rev_primer  # Sequence of REV primer
  
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  # Prefilter to remove reads with ambiguous (N) bases ####
  fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
  fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
  filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # on Windows, set multithread = FALSE
  
  
  path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  fnFs.cut <- file.path(path.cut, basename(fnFs))
  fnRs.cut <- file.path(path.cut, basename(fnRs))
  
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", "--cores 0", # -n 2 required to remove FWD and REV from reads
                                 "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                                 fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
  
  # Discover primer matches, regardless of orientation ####
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])) %>% 
    print()
  
}


# run_itsxpress() ####
# Isolate ITS region 

# Valid taxa_group arguments: 
# Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All

run_itsxpress <- function(directory, # where cutadapted reads live
                          itsregion="ITS1", # must be "ITS1" or "ITS2"
                          taxa_group="All",
                          nthreads=(parallel::detectCores()-1),
                          fwd_pattern="_R1_"){
  
  # find the "cutadapted" files
  fwds <- list.files(directory,pattern = fwd_pattern,full.names = TRUE)
  # build names for outfiles
  outs <- paste0(tools::file_path_sans_ext(fwds) %>% 
                   tools::file_path_sans_ext(), 
                 "_ITS.fastq.gz")
  its_dir <- directory %>% str_replace('cutadapt','ITSx')

  if(!dir.exists(its_dir)){dir.create(its_dir)}
  
  outs <- file.path(its_dir,basename(outs))
  
  # build the ITSxpress command and run it on each file in turn

  for(i in 1:length(fwds)){
      itsxpress <- paste0("itsxpress --fastq ",fwds[i],
                        " --outfile ",outs[i],
                        " --region ",itsregion,
                        " --taxa ",taxa_group,
                        " --threads ",nthreads,
                        " --log ",outs[i],".log",
                        " --single_end")
    
    system(command = itsxpress)
  }
  
}