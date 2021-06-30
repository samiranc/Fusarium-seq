#Importing GenBank data---------------------------------------------------------
  
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

# Formatting Consensus----------------------------------------------------------

### Make sure to download and move consensus file to R project folder!!!! ###

# path for downloaded file to R
path <- file.path(getwd(), "Consensus-22Jun2021.fasta")
#path2 <- system.file(path)

# to make a new file of consensus lowercase
Fusarium_1 <- readLines(path)
Fusarium_1_lower <- tolower(Fusarium_1)
# to save lower case consensus
writeLines(Fusarium_1_lower, "consensus_lower.fasta")

### I copy pasted consensus_lower.fasta to the bottom of Fusarium_16.fasta ###
###                     and saved it!!!                                    ###

# Multiple sequencing using msa-------------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
library(msa)

#Downloading fasta file already downloaded in desktop with Biostrings
#BiocManager::install("Biostrings")
library(Biostrings)

# to locate texshade file
#system.file("tex", "texshade.sty", package = "msa")

path2 <- file.path(getwd(), "Fusarium_16.fasta")

Fusarium_sequences <- readAAStringSet(path2)

Fusarium_sequences

Fusarium_alignment <- msa(Fusarium_sequences)

#Or

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalW") # Same as last

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalOmega") #To protein???

Fusarium_alignment <- msa(Fusarium_sequences, "Muscle") # WORKED???????????????

Fusarium_alignment

print(Fusarium_alignment, show = "complete")
