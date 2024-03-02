#load in relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

#Read the csv into a dataframe 
MethylPotential <- read.csv("metabolism_summary_trimmed_genes.csv")
#Read the taxonmy for the MAGs into a dataframe
MAG_taxonomy <- read.csv("MAG_taxonomy.csv")

#Group by Methyl-Potential and sum the values for all columns of MAGs
#Using where(is.numeric) will only sum the columns that are numeric. More robust than manually selecting the columns
MethylPotentialSummed <- MethylPotential %>% group_by(module) %>% 
  summarise(across(where(is.numeric), sum, na.rm =  TRUE))

#Transpose MethylPotentialSummed to switch the rows and columns
TransMethylPotentialSummed <- as.data.frame(t(MethylPotentialSummed))

#Remove first row and move categories to column headers
colnames(TransMethylPotentialSummed) <- c("Methyl-N", "Methyl-O", "Methyl-S")
TransMethylPotentialSummed <- TransMethylPotentialSummed[-1, ]

#Save a vector as the row names of Trans Methyl and add the vector as a column. Move the column to the first position
mag_names <- row.names(TransMethylPotentialSummed)
TransMethylPotentialSummed$MAGs <- mag_names
TransMethylPotentialSummed <- TransMethylPotentialSummed %>% select(MAGs, everything())

#Melt the dataframe into a long format for use in a heatmap in ggplot using the new column "MAGs" as the identifier to sort by
MeltMethylPotential <- melt(TransMethylPotentialSummed, id.vars = "MAGs")

#Convert values to a numeric vector
MeltMethylPotential$value <- as.numeric(as.character(MeltMethylPotential$value))

#Adds the taxonomy from MAG_taxonomy to MeltMethylPotential based on when MAGs are the same
Taxonomy_Methyl_Potential <- MeltMethylPotential %>% left_join(MAG_taxonomy, by = "MAGs")

#Add a column with combined MAG and taxonomy labels
Taxonomy_Methyl_Potential$combined_labels <- paste(Taxonomy_Methyl_Potential$MAGs, Taxonomy_Methyl_Potential$Taxonomy, sep = "\n")

#Plot heatmap
ggplot(Taxonomy_Methyl_Potential, aes(x=variable, y=combined_labels , fill=value))+
  geom_tile(color='black', lwd=.15, linetype=1)+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt")
  ) +
  scale_fill_continuous(low="white", high="purple", na.value="white")+
  xlab(NULL) + labs(title = "Methylotrophic Substrate Potential", fill = "# of Genes")+
  scale_y_discrete(limits=rev)+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))

#Create vector specifying the MtxB genes 
MtxB_Genes <- c("MttB", "MtbB", "MtmB", "MtgB", "MttB10", "MV8460", 
                "MtpB", "MtcB", "MtyB", "MtsA", "MtpA", "MtsD", "MtsF", 
                "MtsH", "MtaB", "MtoB", "MtvB", "OdmB", "VmdB", "Dhaf_4610", "ELI_2003")
#I could not tell you how the fuck this works, but by god does it work. Filters 
#out any row that doesn't include at least one character from Mtx_B genes in the gene_description column
MtxB_Filtered_Methyl_Potential <- MethylPotential[grepl(paste(MtxB_Genes, collapse = "|"),MethylPotential$gene_description), ]

#Group by Methyl-Potential and sum the values for all columns of MAGs
#Using where(is.numeric) will only sum the columns that are numeric. More robust than manually selecting the columns
MtxB_Filtered_Methyl_Potential <- MtxB_Filtered_Methyl_Potential %>% group_by(module) %>% 
  summarise(across(where(is.numeric), sum, na.rm =  TRUE))

#Transpose MethylPotentialSummed to switch the rows and columns
MtxB_Filtered_Methyl_Potential <- as.data.frame(t(MtxB_Filtered_Methyl_Potential))

#Remove first row and move categories to column headers
colnames(MtxB_Filtered_Methyl_Potential) <- c("Methyl-N", "Methyl-O", "Methyl-S")
MtxB_Filtered_Methyl_Potential <- MtxB_Filtered_Methyl_Potential[-1, ]

#Save a vector as the row names of Trans Methyl and add the vector as a column. Move the column to the first position
mag_names <- row.names(TransMethylPotentialSummed)
MtxB_Filtered_Methyl_Potential$MAGs <- mag_names
MtxB_Filtered_Methyl_Potential <- MtxB_Filtered_Methyl_Potential %>% select(MAGs, everything())

#Melt the dataframe into a long format for use in a heatmap in ggplot using the new column "MAGs" as the identifier to sort by
MtxB_Filtered_Methyl_Potential <- melt(MtxB_Filtered_Methyl_Potential, id.vars = "MAGs")

#Convert values to a numeric vector
MtxB_Filtered_Methyl_Potential$value <- as.numeric(as.character(MtxB_Filtered_Methyl_Potential$value))

#Adds the taxonomy from MAG_taxonomy to MeltMethylPotential based on when MAGs are the same
MtxB_Filtered_Methyl_Potential <- MtxB_Filtered_Methyl_Potential %>% left_join(MAG_taxonomy, by = "MAGs")

#Add a column with combined MAG and taxonomy labels
MtxB_Filtered_Methyl_Potential$combined_labels <- paste(MtxB_Filtered_Methyl_Potential$MAGs, MtxB_Filtered_Methyl_Potential$Taxonomy, sep = "\n")

#Plot heatmap
ggplot(MtxB_Filtered_Methyl_Potential, aes(x=variable, y=combined_labels , fill=value))+
  geom_tile(color='black', lwd=.15, linetype=1)+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt")
  ) +
  scale_fill_continuous(low="white", high="purple", na.value="white")+
  xlab(NULL) + labs(title = "MtxB Genes Present", fill = "# of Genes")+
  scale_y_discrete(limits=rev)+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
