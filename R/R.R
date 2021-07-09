#Importing GenBank data---------------------------------------------------------
  
# R package for phylogenentics and comparative methods
#install.packages("ape")
library(ape)

# Package for nucleotide sequence management
#install.packages("seqinr")
library(seqinr)

# Creating a vector of GenBank accession numbers we want
Fusarium_accession <- c("AY320087", "GU170559", "FJ985328", "AY527534",
                        "AF008506", "FJ985326", "KX434896", "FJ985388",
                        "FJ985365", "AF008485", "AF008491", "FJ985412",
                        "FJ985331", "FJ985312", "AF008490", "FJ985270",
                        "FJ985395", "FJ985330", "FJ985300",
                        "AF008509", "FJ985407", "MF684761", "GQ154455",
                        "GQ154462", "GQ154458", "GQ154456", "GQ154460",
                        "GQ154461", "GQ154454", "FJ985337")

# Reads squences and place them in a DNAbin object
Fusarium_fasta <- read.GenBank(Fusarium_accession)

# List of DNAbin objects
str(Fusarium_fasta)

# Attributes and descriptions
attributes(Fusarium_fasta)
attr(Fusarium_fasta, "description")

# Fasta file format imported to R project
write.dna(Fusarium_fasta, file = "Fusarium_fasta.fasta", format = "fasta",
          append = FALSE, nbcol = 6, colsep = " ", colw = 10)


# To change fasta sequence names into NRRL numbers
#install.packages("phylotools")
library("phylotools")

old_name <- get.fasta.name("Fusarium_fasta.fasta")
new_name <- c("31852", "38302", "36114", "25603",
              "22550", "36107", "26029", "38338",
              "38279", "26035", "26022", "38585",
              "36120", "32882", "25609", "22543",
              "38445", "36118", "26448",
              "26380", "38542", "26622", "53544",
              "53543", "53541", "46589", "53542",
              "53540", "46585", "36251")
ref2 <- data.frame(old_name, new_name)
rename.fasta(infile = "Fusarium_fasta.fasta", ref_table = ref2,
             outfile = "Fusarium_renamed.fasta")

# Formatting NRRL 26992, wasn't found in the GenBank----------------------------
# Given by Dr. Elliot: http://isolate.fusariumdb.org/sequence.php?a=dv&id=3926

NRRL26992 <- "GACAAGACTCACCTTAACGTCGTCGTCATCGGCCACGTCGACTCTGGCAAGTCGACCACTGTGAGTACTCTCCTCGACAA
TGAGTTTATCTGCCATCGTCAATCCCGACCAAGACCTGGTGGGGTATTTCTCAAAGTCAACATACTGACATCGTTTCACA
GACCGGTCACTTGATCTACCAGTGCGGTGGTATCGATAAGCGAACCATCGAGAAGTTCGAGAAGGTTAGTCACTTTCCCT
TCGATCGCGCGTCCTTTGCCCATCGATTTCCCCTACGACTCGAAGCGTGCCCGCTACCCCGCTCGAGACCAAGAATCTTG
CAATATGACCGTAATTTTTTTGGTGGGGCACTTACCCCGCCACTTGAGCGACGGGAGCGTTTGCCCTCTTAACCATTCTC
ACAACCTCAATGAGTGCGTCGTCACGTGTCAAGCAGTCACTAACCATTCAACAATAGGAAGCCGCTGAGCTCGGTAAGGG
TTCCTTCAAGTACGCCTGGGTTCTTGACAAGCTCAAGGCCGAGCGTGAGCGTGGTATCACCATCGATATTGCTCTCTGGA
AGTTCGAGACTCCTCGCTACTATGTCACCGTCATTGGTATGTTGTCGCTCATGCTTCATTCTACTTCTCTTCGTACTAAC
ATATCACTCAGACGCTCCCGGTCACCGTGATTTCATCAAGAACATGA"

## Copy paste NRRL26992 under Fusarium_renamed.fasta and delete the /'s

# Formatting Queen Palm sequence------------------------------------------------
# Download and move Queen Palm file to R project folder

# Path for downloaded file to R
path <- file.path(getwd(), "Queen_Palm_DNA", "EC_006_S3-2_ITS4_B06.ab1.gb")

# To make a new file of queen palm lowercase and have no whitespace or numbers
Fusarium_1 <- readLines(path)
Fusarium_1_lower <- tolower(Fusarium_1)
Fusarium_1_lower <- gsub(" ", "", Fusarium_1_lower)
Fusarium_1_lower <- gsub("[[:digit:]]+", "", Fusarium_1_lower)

# to save lower case consensus
writeLines(Fusarium_1_lower, "EC_lower")

# Copy paste EC_lower to the bottom of Fusarium_renamed.fasta and save it

# Multiple sequencing using msa-------------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
library(msa)

#Downloading fasta file already downloaded in desktop with Biostrings
#BiocManager::install("Biostrings")
library(Biostrings)

path2 <- file.path(getwd(), "Fusarium_renamed.fasta")
Fusarium_sequences <- readAAStringSet(path2)

Fusarium_sequences

Fusarium_alignment <- msa(Fusarium_sequences, "Muscle")

Fusarium_alignment

print(Fusarium_alignment, show = "complete")

# Phylogenetic Tree from msa----------------------------------------------------

Fusarium_alignment2 <- msaConvert(Fusarium_alignment, type ="seqinr::alignment")

#Compute distance matrix from seqinr package
d <- dist.alignment(Fusarium_alignment2)
as.matrix(d)[2:5, drop = FALSE]

# Phlyogenetic tree from nj() function from ape package
FusTree <- nj(d)
plot(FusTree, main = "Phylogenetic Tree of Fusarium oxysporum")

# Phylogenetic Tree using phangorn----------------------------------------------
#install.packages("phangorn")
library("phangorn")

data <- as.phyDat(Fusarium_alignment2, type = "DNA")
Fus_acctran <- acctran(FusTree, data)
plot(midpoint(Fus_acctran))
