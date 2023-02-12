# DADA2_FASTQ-to-Phyloseq
This R code will parse raw FASTQ files for an amplicon sequencing dataset and generate a phyloseq object.

```{r}

# DADA2 workflow for processing 16S raw sequences
# Follows https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(ShortRead)
library(Biostrings)


#-#-#-#-#-#-#-

path <- "PR+Hum+CI/16S"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE)) # If set to be FALSE, then working directory must contain the files
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))

# Remove any forward files that don't have reverse counterparts, and vise versa
# (filterAndTrim will throw an error if fnFs and fnRs have any mismatches)
basefilenames_Fs <- sub("_L001_R1_001.fastq.gz","",basename(fnFs))
basefilenames_Rs <- sub("_L001_R2_001.fastq.gz","",basename(fnRs))
rm_from_fnFs <- basefilenames_Fs[which(!(basefilenames_Fs %in% basefilenames_Rs))]
rm_from_fnRs <- basefilenames_Rs[which(!(basefilenames_Rs %in% basefilenames_Fs))]

for(name in rm_from_fnFs) {
  print(paste(name, "does not have a reverse-reads counterpart. Omitting from this analysis."))
  fnFs <- fnFs[-which(fnFs == paste0(path, "/", name, "_L001_R1_001.fastq.gz"))]
}
for(name in rm_from_fnRs) {
  print(paste(name, "does not have a forward-reads counterpart. Omitting from this analysis."))
  fnRs <- fnRs[-which(fnRs == paste0(path, "/", name, "_L001_R2_001.fastq.gz"))]
}


# Identify primers - used 515F & 806R from Hiro's spreadsheet
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer seq
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...


# Get all orientations of primers, just to be safe
# Note - changed this for the dimensions project due to different sequencing primers
# Used

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
FWD.orients

# “pre-filter” the sequences just to remove those with Ns, but perform no other filtering
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 1, multithread = TRUE)

# (From tutorial) We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations. 
# Identifying and counting the primers on one set of paired end FASTQ files is
# sufficient, assuming all the files were created using the same library preparation,
# so we’ll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# If you see the reverse-complement of the forward primer in the reverse reads (cells [2,4] and [3,4]),
# it's because the ITS region is short and it is reading part of the forward primer.

# Remove primers using cutadapt

cutadapt <- "miniconda3/condabin/cutadapt.exe" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
# for(i in 1:10) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files; fnFs.filtN replaced by fnFs.filtN, etc.
                             "--minimum-length", "1")) # min length of cutadapted reads: >0 
}

# Count primers in first post-cutadapt sample (should all be 0):
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[50]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[50]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[50]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[50]]))

# Since they are zero, skip step to remove other orientations of primers

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), split="_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles of forward reads #1-2
plotQualityProfile(cutFs[1:12])

# Inspect read quality profiles of reverse reads #1-2
plotQualityProfile(cutRs[1:12])

# Filter and trim

# Assigning the filenames for the output of the filtered reads 
# to be stored as fastq.gz files.
filtFs <- file.path(path, "filtered", basename(fnFs.filtN))
filtRs <- file.path(path, "filtered", basename(fnRs.filtN))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

filtFs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R1_001.fastq.gz", full.names=TRUE)
filtRs.out <- list.files(paste(path, "filtered", sep="/"), pattern="_L001_R2_001.fastq.gz", full.names=TRUE)

# Learn the error rates
errF <- learnErrors(filtFs.out, multithread = TRUE)
errR <- learnErrors(filtRs.out, multithread = TRUE)

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs.out, verbose = TRUE)
derepRs <- derepFastq(filtRs.out, verbose = TRUE)
# Name the derep-class objects by the sample names
get.sample.name <- function(fname) {
  paste(strsplit(basename(fname), "_")[[1]][1:3], collapse="_")
}
sample.names <- unname(sapply(filtFs.out, get.sample.name))

# DADA2's core sample inference algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

          
# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,trimOverhang = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

saveRDS(seqtab,"PR+Hum+CI_seqtab.bac.rds")
saveRDS(seqtab.nochim,"PR+Hum+CI_seqtab.bac.nochim.rds")

# Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))

#format out to accommodate dropped samples
raw.sample.names <- unname(sapply(row.names(out), get.sample.name))

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

track2<-cbind(out,track[match(row.names(out),row.names(track)),])

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
write.csv(track2,"PR+Hum+CI_bac_summary.csv")

rownames(track) <- sample.names
head(track2)

```

```{r combine runs}

# Assign taxonomy using the UNITE database
silva.ref<-"silva_nr99_v138_train_set.fa.gz"
silva.species<-"silva_species_assignment_v138.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = TRUE, tryRC = TRUE)

taxa <- addSpecies(taxa, silva.species)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Save OTU table and taxonomic table as RDS files
# to hand off to dada2_to_phyloseq.R
saveRDS(seqtab.nochim, "PR+Hum+CI_16S_seqtab_nochim_taxa_3-2-22")
saveRDS(taxa, "PR+Hum+CI_16S_taxa_3-2-22.rds")


```

```{construct phyloseq object with sample data variables}

#read in environmental sample table

sample <- readRDS("PR+Hum+CI_16S_seqtab_nochim_taxa_3-2-22.rds")
SD <- read.csv("PR+H+CI_16S_SampleData.csv", row.names = 1)
SAM <- sample_data(SD, errorIfNULL = T)
taxa <- readRDS("PR+Hum+CI_16S_taxa_3-2-22.rds")

#create a phyloseq object
bac.ps <- phyloseq(otu_table(sample, taxa_are_rows=FALSE), sample_data(SAM), tax_table(taxa))

#filter out unwanted taxa (e.g., mitochondira and chloroplast sequences)
bac.ps.filt<-subset_taxa(bac.ps,Family!="Mitochondria")
bac.ps.filt<-subset_taxa(bac.ps.filt,Genus!="Chloroplast")

# Removing sequence rownames for display only
taxa.print <- tax_table(bac.ps.filt)
rownames(taxa.print) <- NULL
head(taxa.print)

#save the filtered dataset 
saveRDS(bac.ps.filt,"PHC.16S.filtered.3-2-22.rds")

#filter out low abundant sequences
bac.ps.filt2 = prune_taxa(taxa_sums(bac.ps.filt) > 10, bac.ps.filt) 
bac.ps.filt2 = prune_samples(sample_sums(bac.ps.filt2)>1000, bac.ps.filt2)

#save the filtered+pruned dataset
saveRDS(bac.ps.filt2, "PHC.16S.filt-prune.rds")

#rarefy the dataset
bac.ps.rare <- rarefy_even_depth(bac.ps.filt2, sample.size = min(sample_sums(bac.ps.filt2)),
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#save rarefied dataset
saveRDS(bac.ps.rare, "Rarefied-PHC-Bacteria.rds")
