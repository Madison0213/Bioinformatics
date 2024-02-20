
#Midterm 1: Take Home Practical
#Madison Stevens
#2024-02-19

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
library(protti)
library("GenomicAlignments")
library("UniprotR")
library("r3dmol")

# 1. Import and align your DNA sequences
combinedseq <- readDNAStringSet("Data/sequences.fasta")
combinedseq
print(combinedseq, show="complete")

#Number 6 has 641bp while all the others have 642bp

#### Creating an msa alignment ####
msaalignment <- msa(combinedseq)
msaalignment

print(msaalignment, show="complete")

---------------------------------------------------------
# 2. Are any of the samples different from the rest? 

# Identify variations or mutations
mutations <- consensus != combinedseq
print("Mutations:")
print(mutations)
# (1)FALSE (2)FALSE (3)FALSE  (4)TRUE (5)FALSE  (6)TRUE (7)FALSE (8)FALSE 
#(9)FALSE  (10)TRUE (11)FALSE (12)FALSE (13)FALSE (14)FALSE (15)FALSE (16)FALSE 
#(17)FALSE (18)FALSE (19)FALSE (20)FALSE
# 4, 6, and 10 have mutations

#If so, what kinds of mutations do you observe in this individual (or individuals)?

#most likely silent mutations

--------------------------------------------------------------------
# 3. What is the gene? hbb gene for beta globin
# What is the accession number of the best match to your search? LC121775
# I used individual 6 to find this accession number
# I used BLAST to find the gene and then GenBank for the accession number
--------------------------------------------------------------------
# 4. Who is the most different?

#Create a distance matrix
MSACom <- msaConvert(msaalignment, type="seqinr::alignment")
d <- dist.alignment(MSACom)
print(d)

#Visualize in Phylogenetic tree
HSTree <- nj(d)
plot(HSTree, main="Phylogenetic Tree of Homo sapains Gene Sequences")
#phylogenetic tree is in the output folder
#Individual 6 is the most different

# Choose one DNA sequence to translate to protein
individual6 <- combinedseq[[6]]  
# Replace with the index of the sequence you want to translate

# Translate the selected DNA sequence to a protein sequence
amino_acid_sequences <- Biostrings::translate(individual6)
amino_acid_sequences
AA_seq <- as.character(amino_acid_sequences)

output_file <- "individual6.fasta"

writeXStringSet(AAStringSet(AA_seq), file = output_file,
                format = "fasta", width = 60)
#file is named "individual6" and is in the data folder
---------------------------------------------------------------
# 5. Uniprot for protein match and accession number
  #Hemoglobin subunit beta
  #A0A0J9YWK4
---------------------------------------------------------------
# 6. Diseases associated with this gene? Does individual 6 have it?
#name accession number to variable
accession_number <- ("A0A0J9YWK4")

# Use GetPathology_Biotech() and Get.diseases() to find information on any diseases or pathologies associated with your gene ####
GP<-GetPathology_Biotech(accession_number)
GP

Get.diseases(GP, directorypath = getwd())
#NULL
#There are no diseases associated with this gene.

#I went on OMIM and searched the "Hemoglobin subunit beta"
#I found the diseases: thalassemia, erythrocytosis, heinz body anemia, 
#methemoglobinemia, and sickle cell disease.
#I don't think we can know what disease they have, but they could have any one 
#of these diseases.
---------------------------------------------------------------
# 7. 3-D structure of protein 
FAlpha<-fetch_alphafold_prediction(accession_number)
View(FAlpha)  

#I had to use https://alphafold.ebi.ac.uk/
# added screenshot to output folder
  
  




---------------------------------------------------------------  
  # CLEAN UP #####
# Clear environment
rm(list = ls()) 
# Clear packages
# requires the package pacman to work# CLEAN UP #####
# Clear environment
p_unload(all)  # Remove all add-ons
# Clear console
cat("\014")  # ctrl+L
# Clear mind :)
