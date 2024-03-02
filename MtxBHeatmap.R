#load in relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

#Read the csv into a dataframe and change one of the column headers to be easier to read
MethylPotential <- read.csv("TR_Metabolism_Summary_METHYL_MtxBTrim.csv")
colnames(MethylPotential)[colnames(MethylPotential) == "X.Rumen.1."] <- "Rumen.1"

#Group by Methyl-Potential and sum the values 
MethylPotentialSummed <- MethylPotential %>% group_by(module) %>% 
  summarise(Rumen_1 = sum(Rumen.1), Rumen_2 = sum(Rumen.2), 
            Rumen_3 = sum(Rumen.3), Rumen_4 = sum(Rumen.4), Wetland_1 = sum(Wetland.1), Wetland_2 = sum(Wetland.2))

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
#Change MAGs so they look nicer on the plot 
MeltMethylPotential$MAGs <- c("Rumen 1", "Rumen 2", "Rumen 3", "Rumen 4", "Wetland 1",
                              "Wetland 2", "Rumen 1", "Rumen 2", "Rumen 3", "Rumen 4", 
                              "Wetland 1", "Wetland 2","Rumen 1", "Rumen 2", "Rumen 3", 
                              "Rumen 4", "Wetland 1", "Wetland 2")

#Plot heatmap

ggplot(MeltMethylPotential, aes(x=variable, y=MAGs , fill=value))+
  geom_tile(color='black', lwd=.15, linetype=1)+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt")
  ) +
  scale_fill_continuous(low="white", high="purple", na.value="white")+
  xlab(NULL) + labs(title = "MtxB Genes Present", fill = "# of Genes")+
  scale_y_discrete(limits=rev)
