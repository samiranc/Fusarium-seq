# Importing data----------------------------------------------------------------

# R package for phylogenentics and comparative methods
install.packages("ape")
library(ape)

# package for nucleotide sequence management
install.packages("seqinr")
library(seqinr)

sessionInfo()

# Creating a vector of GenBank accession numbers we want
Fusarium_accession <- c("GQ154454", "GQ154455", "GQ154456", "GQ154457",
                        "GQ154458", "GQ154459", "GQ154460", "GQ154461",
                        "GQ154462", "GQ154463", "GQ154464", "GQ154465",
                        "GQ154466", "GQ154467", "GQ154468", "GQ154469")

??read.GenBank

# read squences and place them in a DNAbin object
Fusarium_sequences <- read.GenBank(Fusarium_accession)

# list of DNAbin objects
str(Fusarium_sequences)

# species list
attr(Fusarium_sequences, "species")

# fatsa file format
?write.dna

write.dna(Fusarium_sequences, file = "Fusarium_fasta_1.fasta", format = "fasta",
          append = FALSE, nbcol = 6, colsep = " ", colw = 10)

#
path <- file.path(getwd(), )
