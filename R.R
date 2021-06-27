# Importing GenBank data--------------------------------------------------------

# R package for phylogenentics and comparative methods
#install.packages("ape")
library(ape)

# Package for nucleotide sequence management
#install.packages("seqinr")
library(seqinr)

# Creating a vector of GenBank accession numbers we want
Fusarium_accession <- c("GQ154454", "GQ154455", "GQ154456", "GQ154457",
                        "GQ154458", "GQ154459", "GQ154460", "GQ154461",
                        "GQ154462", "GQ154463", "GQ154464", "GQ154465",
                        "GQ154466", "GQ154467", "GQ154468", "GQ154469")

# Reads squences and place them in a DNAbin object
Fusarium_16 <- read.GenBank(Fusarium_accession)

# List of DNAbin objects
str(Fusarium_16)

# Species list
attr(Fusarium_16, "species")

# Fatsa file format imported to R project
write.dna(Fusarium_16, file = "Fusarium_16.fasta", format = "fasta",
          append = FALSE, nbcol = 1, colsep = " ", colw = 80)

#To make accession numbers into species name, but all named Fusarium_oxysporum??
#Fusarium_sequences_GenBank_IDs <- paste(attr(Fusarium_sequences, "species"), 
# names (Fusarium_sequences), sep ="_RAG1_")

#Fusarium_16 <- read.fasta(file = "Fusarium_fasta_1.fasta", seqtype = "DNA",
#                          as.string = TRUE, forceDNAtolower = FALSE)

#write.fasta(sequences = Fusarium_16, names = Fusarium_sequences_GenBank_IDs,
#            nbchar = 10, file.out = "Fusarium_seq_seqinr_format.fasta")

# Formatting Consensus----------------------------------------------------------

# Make sure to download and move consensus file to R project folder!!!!
# path for downloaded file to R
path <- file.path(getwd(), "Consensus-22Jun2021.fasta")
#path2 <- system.file(path)
# to make a new file of consensus lowercase
Fusarium_1 <- readLines(path)
Fusarium_1_lower <- tolower(Fusarium_1)

writeLines(Fusarium_1_lower, "consensus_lower.fasta")

path2 <- file.path(getwd(), "consensus_lower.fasta")


# WAIT FUSARIUM-1 IS A CHARACTER SO ITS CHARACTER TO DNABIN???
install.packages("adegenet")
library("adegenet")
Fusarium_1 <- fasta2DNAbin(path, quiet = FALSE, chunkSize = 10, ###fasta2DNAbin taking too long????
                           snpOnly = FALSE)
#or
Fusarium_1 <- fasta2DNAbin(path, chunkSize = 10)

# Combine either DNAbin or Fasta
#Fusarium_16_fasta <- read.dna("Fusarium_16.fasta", format = "fasta") ################# combine as fasta
#Fusarium_1_fasta <- read.dna("Consensus-22Jun2021.fasta", format = "fasta") ############################ combine as fasta



### I copy pasted consensus_lower.fasta to the bottom of Fusarium_16.fasta and 
### saved it!!!!

# Multiple sequencing using msa-------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
library(msa)
#Downloading fasta file already downloaded in desktop with Biostrings
BiocManager::install("Biostrings")
library(Biostrings)

# to locate texshade file
system.file("tex", "texshade.sty", package = "msa")

Fusarium_seq_file <- system.file("Fusarium_16", "Fusarium_16.fasta",#### what is "example"????
                                 package = "msa")
View(Fusarium_seq_file)
Fusarium_sequences <- readAAStringSet(Fusarium_seq_file)

Fusarium_alignment <- msa(Fusarium_sequences)

# NRRL strain instead of accession numbers
