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

# Species list
attr(Fusarium_fasta, "species")

# Fatsa file format imported to R project
write.dna(Fusarium_fasta, file = "Fusarium_fasta.fasta", format = "fasta",
          append = FALSE, nbcol = 6, colsep = " ", colw = 10)

# Formatting Consensus----------------------------------------------------------

### Make sure to download and move consensus file to R project folder     ###
### or documents!!!!                                                      ###

# Path for downloaded file to R
path <- file.path(getwd(), "EC_006_S3-2_ITS4_B06.ab1.gb")

# To make a new file of consensus lowercase
Fusarium_1 <- readLines(path)
Fusarium_1_lower <- tolower(Fusarium_1)
# to save lower case consensus
writeLines(Fusarium_1_lower, "EC_lower")

### I copy pasted consensus_lower.fasta to the bottom of Fusarium_fasta    ###
###                      and saved it!!!                                   ###

# Multiple sequencing using msa-------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
library(msa)

#Downloading fasta file already downloaded in desktop with Biostrings
#BiocManager::install("Biostrings")
library(Biostrings)

path2 <- file.path(getwd(), "Fusarium_fasta.fasta")

Fusarium_sequences <- readAAStringSet(path2)

Fusarium_sequences

Fusarium_alignment<- msa(Fusarium_sequences)

#Or

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalW") # Same as last

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalOmega") #To protein???

Fusarium_alignment <- msa(Fusarium_sequences, "Muscle") # WORKED???????????????

Fusarium_alignment

print(Fusarium_alignment, show = "complete")

# To install tinytex to use msaPrettyPrint function
install.packages("tinytex")
library(tinytex)
install.packages("knitr")
library(knitr)
# to locate texshade file
p <-system.file("tex", "texshade.sty", package = "msa")

msaPrettyPrint(Fusarium_alignment, output = "pdf", showNames = "none",
               showLogo = "none", askForOverwrite = FALSE)
texify.exe("Fusarium_alignment.tex", clean=TRUE)

texi2pdf("Fusarium_alignment.tex", clean = FALSE, quiet = TRUE,
         texi2dvi = getOption("texi2dvi"),
         texinputs = NULL, index = TRUE)

path5

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Program Files/MiKTeX 2.9/miktex/bin/x64", sep=.Platform$path.sep))
