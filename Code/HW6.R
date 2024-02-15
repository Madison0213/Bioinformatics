
#Madison Stevens
#HW 6
#2/15/24

#I used genes from the Palm Warbler on HW5 and used the amino acid seq.

BiocManager::install("GenomicAlignments")
install.packages(UniprotR)
install.packages(protti)
install.packages("r3dmol")
library("GenomicAlignments")
library("UniprotR")
#library("protti") won't work!
library("r3dmol")
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
library(seqinr)
library(ape)
library(phangorn)
library(devtools)
devtools::install_github("jpquast/protti", dependencies = TRUE)
library(protti)

setwd("Data/")
setwd("Bioinformatics/")
setwd("Homework06/")
#Set Working directory for Data and Homework folder ####

getwd()
#Checked that directory was in place ####

#read DNA seq
mySequences01 <- readDNAStringSet("sequence1.fasta")

# 3. transform into AA seq
amino_acid_sequences <- Biostrings::translate(mySequences01)
amino_acid_sequences
as.character(amino_acid_sequences)

output_file <- "amino_acid_sequence.fasta"
writeXStringSet(amino_acid_sequences, file = output_file,
                format = "fasta", width = 60)

#4. Read this file into R using the appropriate function ####
#accession_numbers<- read.table("Accession_numbers.txt")
#cannot open file, won't work

#5. Sample list of accession numbers ####
accession_numbers <- c("S0ASK7", "O21399", "Q94WR7", "A0A068L9A8", "P24984")

# Convert the list to a character string
accession_string <- paste(accession_numbers, collapse = ",")

# Print the formatted string
print(accession_string)

#6. Reading accession numbers into GetProteinGOInfo ####
AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
str(AccessionNumbersGO)
#write.csv(AccessionNumbersGO, "AccessionNumbersGO.csv", row.names = FALSE)

#7. Plot your results
PlotGoInfo(AccessionNumbersGO)
#Stop here, but if plotgo does not work, do the rest.

#8. Save plot of GO terms
#saved to the output folder within this repository

#9. interesting GO terms of your gene
#Biological Processes:oxidative phosphorylation
#Molecular Functions:metal ion binding, heme binding, cytochrome-c oxidase activity
#Cellular Components:REspiratory chain complex IV, mitochonrail inner membrane,
#mitochondrial respiratory chain complex III & IV

#10. Use GetPathology_Biotech() and Get.diseases() to find information on any diseases or pathologies associated with your gene ####
GP<-GetPathology_Biotech(accession_numbers)
GP
#need as variable to use for get.diseases
#NA on all counts
Get.diseases(GP, directorypath = getwd())
#NA


#11. We are going to access structural information using the protti package ####
viewtibble <- fetch_uniprot(accession_numbers)
View(viewtibble)


#12. Pull any available structural information from the Protein DataBase
#NA in table xref_pdb
FP1<-fetch_pdb("1ZMR")
view(FP1)

FP2<-fetch_pdb("2HWG")
view(FP2)


#13. Get information on any available 3D structures for your gene
FAlpha<-fetch_alphafold_prediction(accession_numbers)
View(FAlpha)
#image of S0ASK7 is in the output folder