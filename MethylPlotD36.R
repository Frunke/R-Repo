#load in relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

#Set working directory to location of csv and load csv into dataframe
setwd("C:/Users/teaga/Downloads/Grad School/WrightonLab/")
MethaneDF <- read.csv("TML_Enrichment_GC_2024.csv")

#Move the top row to the column names and remove the redundant row
colnames(MethaneDF) <- MethaneDF[1, ]
MethaneDF <- MethaneDF[-1, ]

#Create a vector that stores "TRUE" if the column header is a CH4 Area (Also the first column)
#Use the vector to make a subset that only has CH4 area
CH4_Columns <- grepl("^CH4 Area", names(MethaneDF))
CH4_Columns[1] <- TRUE
Methane_subset <- MethaneDF[, CH4_Columns]

#Add meaningful headers to the dataframe
colnames(Methane_subset) <- c("Sample_Names", 0, 6, 15, 22, 29, 36, 43, 50)

#Save a vector with every column but the first
#Use that vector to convert every column to numeric
columns_to_numeric <- names(Methane_subset)[-1]
Methane_subset[, columns_to_numeric] <- lapply(Methane_subset[, columns_to_numeric], as.numeric)

Fed_subset <- Methane_subset[,-8]

#Sum all rows to a new column. Don't include the first column of sample names
Fed_subset$Total <- rowSums(Fed_subset[ ,-1], na.rm = TRUE)

#Add column for soil grams
Soil_grams <- c(3.06, 3.06, 3.02, 3.27, 3.08, 3.36, 3.25, 3.11, 3.06, 3.38, 2.95, 2.99, 2.95, 3.11, 2.97, 3.05)
Fed_subset$Grams_Soil <- Soil_grams

#Convert to CH4 Percent normalized to grams of soil added
Fed_subset$Soil_Normalized <- (Fed_subset$Total/Fed_subset$Grams_Soil)
Fed_subset$Percent_SoilNormalized <- ((Fed_subset$Soil_Normalized + 65.482) / 41035 ) * 100

#Add categories for the samples so I can group them and average.
Fed_subset$ExpGroup <- c("OWC.Control", "OWC.Control", "OWC.Control", "OWC.Control", 
               "OWC.TML", "OWC.TML", "OWC.TML", "OWC.TML", "STM.Control", 
               "STM.Control", "STM.Control", "STM.Control", "STM.TML", "STM.TML", "STM.TML", "STM.TML")


#Make a new DF (Stripped) from the columns 10,11 in MethaneDF
#Pipe the columns into group_by which organizes them by the ExpGroup
#Summarize the new groups by creating a new column "ave" that is the mean of the groups by Percent_SoilNormalized
#Also make a solumn that is sd which is the standard deviation of the Percent_SoilNormalized
Methane_Cumulative <- Fed_subset[c(11,12)] %>% group_by(ExpGroup) %>% summarise(ave=mean(Percent_SoilNormalized), se=(sd(Percent_SoilNormalized)/sqrt(length(Percent_SoilNormalized))))

#Make the dataplot
#load in ggplot and use the DF 
#geom_col makes the plot with the x axis being the Experimental group and the y axis being the averages
#fill sets the color of the barplot
#theme_minimal changes the background to a clean one
#scale_x_discrete sets the titles of the bars
#geom_errorbar adds error bars to the x-axis with the minimum being the average - standard deviation and
#the maximum being average + standard deviation
#Ylab and xlab > label the axes
ggplot(Methane_Cumulative) + geom_col(aes(x=ExpGroup, y=ave, fill = ExpGroup)) +
  geom_errorbar(aes(x=ExpGroup, ymin=ave-se, ymax=ave+se), width=0.25) + 
  ylab("Cumulative Percent Methane per Gram of Soil") + 
  xlab("Experimental Groups") +
  scale_x_discrete(labels = NULL)+
  scale_fill_manual(values = c(OWC.Control = "#a1d99b",
                               OWC.TML = "#31a354",
                               STM.Control = "#9ecae1",
                               STM.TML = "#3182bd"),
                    name = "",
                    labels = c("OWC Control", "OWC TML", "STM Control", "STM TML")) +
  theme_minimal() + 
  scale_y_continuous(breaks = seq(0,3, by = 0.25))

#LINE PLOT
#Copy the dataframe to preserve the initial data
Methane_lineplot <- Methane_subset
#Divide the values by the grams of soil to normalize
Methane_lineplot$Grams_soil <- Fed_subset$Grams_Soil

Methane_lineplot[,2:8] <- Methane_lineplot[,2:8] / Methane_lineplot$Grams_soil
#Convert to a percentage using Laura's standard curve
Methane_lineplot[,2:8] <- ((Methane_lineplot[,2:8] + 65.842) / 41035) * 100


Methane_lineplot$ExpGroup <- Fed_subset$ExpGroup

#Melt the dataframe to a long format by the sample.names
Methane_lineplot <- melt(Methane_lineplot[c(1:8)],id.vars = "Sample_Names")

#Adds the ExpGroup as a column based on column 1 (Sample_names) and column 12 (ExpGroup) from Methane_subset
#by=c() function means match sample names to sample names to add the expgroup
Methane_lineplot <-  inner_join(Fed_subset[c(1,12)], Methane_lineplot, by=c("Sample_Names"="Sample_Names"))

#Create new df grouped by exp_group and variable and create means and sd
#na.rm = TRUE is vital, otherwise it only sums/se's the values where all 4 replicates have values
Methane_lineplot <- Methane_lineplot %>% group_by(ExpGroup, variable) %>% summarise(mean=mean(value, na.rm = TRUE),se=sd((value)/ sqrt(n()), na.rm = TRUE))

#Change NAs and NaNs to zeros
Methane_lineplot <- Methane_lineplot %>% mutate_at(vars(value), ~ifelse(is.na(.x) | is.nan(.x), 0, .x))

#joins the group from "group.df" to average_time
group <- read.csv("Group.csv")
Methane_lineplot <- inner_join(Methane_lineplot, group, by=c("ExpGroup"="X"))


#Dataplot for linegraph by time
#geom_line calls the line graph and sets the x axis as Days, y axis as the mean, and the group by ExpGroup
#geom_point puts points in the graph, looks nicer
#geom_errorbar adds the error bars using the sd made earlier
#facet_wrap() seperates the line graphs by STM and OWC (A variable)
plotly::ggplotly(
ggplot(Methane_lineplot) +
  geom_line( aes(x=variable, y=value, group=Sample_Names, color=Sample_Names),linetype=1)+
  geom_point(aes(x=variable, y=value, group=Sample_Names, color=Sample_Names), shape=15)+
  theme_classic()+ 
  ylab("Percent Methane per Gram of Soil")+
  xlab("Day of Enrichment")+
  scale_color_manual(values = c(OWC_E08 = "black"))+
  ggtitle("Trimethyllysine (TML) stimulation of CH4 in 2018 M2 soils") + theme(axis.text.x=element_text(angle=45, hjust=1))+
  facet_wrap(~A)
)
