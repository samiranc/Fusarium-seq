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


####
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

#############################################
# To isolate descriptions for NRRL #'s
Fusarium_IDs <- paste(attr(Fusarium_fasta, "description"), names(Fusarium_fasta),
                      sep = "_RAG1_")
Fusarium_IDs

# Read our fasta into seqinr
Fusarium_seqinr <- read.fasta(file = "Fusarium_fasta.fasta", seqtype = "DNA",
                              as.string = TRUE, forceDNAtolower = FALSE)

Fusarium_seqinr

# Rewrtie fasta with description as names
write.fasta(sequences = Fusarium_seqinr, names = Fusarium_IDs,
            nbchar = 10, file.out = "Fusarium_seqinr.fasta")

# Formatting Consensus----------------------------------------------------------

### Make sure to download and move queen of palm file to R project folder ###
###                     or documents!!!!                                  ###

# Path for downloaded file to R
path <- file.path(getwd(), "EC_006_S3-2_ITS4_B06.ab1.gb")

# To make a new file of consensus lowercase
Fusarium_1 <- readLines(path)
Fusarium_1_lower <- tolower(Fusarium_1)
# to save lower case consensus
writeLines(Fusarium_1_lower, "EC_lower")

### I copy pasted EC_lower to the bottom of Fusarium_fasta                 ###
###                      and saved it!!!                                   ###

# Multiple sequencing using msa-------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
library(msa)

#Downloading fasta file already downloaded in desktop with Biostrings
#BiocManager::install("Biostrings")
library(Biostrings)

path2 <- file.path(getwd(), "Fusarium_renamed.fasta")
Fusarium_sequences <- readAAStringSet(path2)

Fusarium_sequences

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalW") # Nothing

Fusarium_alignment <- msa(Fusarium_sequences, "ClustalOmega") #To protein?

Fusarium_alignment <- msa(Fusarium_sequences, "Muscle") # WORKED????

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

Sys.setenv(PATH = paste(Sys.getenv("PATH"),
                        "C:/Program Files/MiKTeX 2.9/miktex/bin/x64",
                        sep=.Platform$path.sep))

# Phylogenetic Tree-------------------------------------------------------------

Fusarium_alignment2 <- msaConvert(Fusarium_alignment, type = "seqinr::alignment")

#Compute distance matrix from seqinr package
library(seqinr)
d <- dist.alignment(Fusarium_alignment2)
as.matrix(d)[2:5, drop = FALSE]

# Phlyogenetic tree from nj() function from ape package
library(ape)
FusTree <- nj(d)
plot(FusTree, main = "Phylogenetic Tree of Fusarium oxysporum")
