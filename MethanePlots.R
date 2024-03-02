#load in relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

#Set working directory to location of csv and load csv into dataframe
setwd("C:/Users/teaga/Downloads/Grad School/WrightonLab/")
MethaneDF <- read.csv("TML_Enrichment_GC_2024.csv")

#Remove non CH4 Area columns and top row
MethaneDF <- MethaneDF[, -c(2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20)]
MethaneDF <- MethaneDF[-1, ]

#Add meaningful headers to the dataframe
colnames(MethaneDF) <- c("Sample Names", "CH4.D00", "CH4.D06", "CH4.D15", "CH4.D22", "CH4.D29")

#Change columns to numeric to allow math. Use class() to ensure that the columns are now correctly converted
MethaneDF$CH4.D00 <- as.numeric(as.character(MethaneDF$CH4.D00))
MethaneDF$CH4.D06 <- as.numeric(as.character(MethaneDF$CH4.D06))
MethaneDF$CH4.D15 <- as.numeric(as.character(MethaneDF$CH4.D15))
MethaneDF$CH4.D22 <- as.numeric(as.character(MethaneDF$CH4.D22))
MethaneDF$CH4.D29 <- as.numeric(as.character(MethaneDF$CH4.D29))
class(MethaneDF$CH4.D00)
class(MethaneDF$CH4.D06)
class(MethaneDF$CH4.D15)
class(MethaneDF$CH4.D22)
class(MethaneDF$CH4.D29)

#Sum all rows to a new column. Don't include the first column of sample names
MethaneDF$CH4.Total <- rowSums(MethaneDF[ ,-1], na.rm = TRUE)

#Add column for soil grams
Soil_grams <- c(3.06, 3.06, 3.02, 3.27, 3.08, 3.36, 3.25, 3.11, 3.06, 3.38, 2.95, 2.99, 2.95, 3.11, 2.97, 3.05)
MethaneDF$Grams.Soil <- Soil_grams

#Convert to CH4 Percent normalized to grams of soil added
MethaneDF$CH4_SoilNormalized <- (MethaneDF$CH4.Total/MethaneDF$Grams.Soil)
MethaneDF$CH4_Per_SoilNormalized <- ((MethaneDF$CH4_SoilNormalized + 65.482) / 41035 ) * 100

#Add categories for the samples so I can group them and average. Then add the vector as a column to the DF
Exp_group <- c("OWC.Control", "OWC.Control", "OWC.Control", "OWC.Control", 
               "OWC.TML", "OWC.TML", "OWC.TML", "OWC.TML", "STM.Control", 
               "STM.Control", "STM.Control", "STM.Control", "STM.TML", "STM.TML", "STM.TML", "STM.TML")
MethaneDF$ExpGroup <- Exp_group

#Make a new DF (Stripped) from the columns 10,11 in MethaneDF
#Pipe the columns into group_by which organizes them by the ExpGroup
#Summarize the new groups by creating a new column "ave" that is the mean of the groups by CH4_Per_SoilNormalized
#Also make a solumn that is sd which is the standard deviation of the CH4_Per_SoilNormalized
MethaneDF.Stripped <- MethaneDF[c(10,11)] %>% group_by(ExpGroup) %>% summarise(ave=mean(CH4_Per_SoilNormalized), sd=sd(CH4_Per_SoilNormalized))

#Make the dataplot
#load in ggplot and use the DF 
#geom_col makes the plot with the x axis being the Experimental group and the y axis being the averages
#fill sets the color of the barplot
#theme_minimal changes the background to a clean one
#scale_x_discrete sets the titles of the bars
#geom_errorbar adds error bars to the x-axis with the minimum being the average - standard deviation and
#the maximum being average + standard deviation
#Ylab and xlab > label the axes
ggplot(MethaneDF.Stripped) + geom_col(aes(x=ExpGroup, y=ave), fill = "forestgreen") +
  geom_errorbar(aes(x=ExpGroup, ymin=ave-sd, ymax=ave+sd), width=0.2) + ylab("Cumulative Average Percent CH4 per gram of Soil") + xlab("Experimental Groups") +
  theme_minimal() + scale_x_discrete(labels=c("OWC Control", "OWC TML", "STM Control", "STM TML"))

#LINE PLOT
#Divide the CH4 area by grams of soil to get the normalized value
#Convert to percentage using Laura's CH4 standard curve
MethaneDF$CH4.D06 <- (MethaneDF$CH4.D06 / MethaneDF$Grams.Soil)
MethaneDF$CH4.D15 <- (MethaneDF$CH4.D15 / MethaneDF$Grams.Soil)
MethaneDF$CH4.D22 <- (MethaneDF$CH4.D22 / MethaneDF$Grams.Soil)
MethaneDF$CH4.D29 <- (MethaneDF$CH4.D29 / MethaneDF$Grams.Soil)
MethaneDF$CH4.D29 <- ((MethaneDF$CH4.D29 + 65.842) / 41035) *100
MethaneDF$CH4.D22 <- ((MethaneDF$CH4.D22 + 65.842) / 41035) *100
MethaneDF$CH4.D15 <- ((MethaneDF$CH4.D15 + 65.842) / 41035) *100
MethaneDF$CH4.D06 <- ((MethaneDF$CH4.D06 + 65.842) / 41035) *100

#change column to make it easier to work with
colnames(MethaneDF)[1] <- "Sample.Names"
#Make a new trimmed DF
TimeMethaneDF <- select(MethaneDF, Sample.Names, CH4.D00, CH4.D06, CH4.D15, CH4.D22, CH4.D29, Grams.Soil, ExpGroup)
#Melt the dataframe to a long format by the sample.names
MeltTime <- melt(TimeMethaneDF[c(1:6)],id.vars = "Sample.Names")
#I don't remember why I needed this
MeltTime <-  inner_join(TimeMethaneDF[c(1,8)], MeltTime, by=c("Sample.Names"="Sample.Names"))

#Create new df grouped by exp_group and variable and create means and sd
average_time <- MeltTime %>% group_by(ExpGroup, variable) %>% summarise(mean=mean(value),sd=sd(value))

#Change NAs to zeros
average_time[is.na(average_time)] <- 0

#Change CH4.D to just D.
average_time$variable <- gsub("CH4.D00", "D00", average_time$variable)
average_time$variable <- gsub("CH4.D06", "D06", average_time$variable)
average_time$variable <- gsub("CH4.D15", "D15", average_time$variable)
average_time$variable <- gsub("CH4.D22", "D22", average_time$variable)
average_time$variable <- gsub("CH4.D29", "D29", average_time$variable)

#joins the group from "group.df" to average_time
average_time <- inner_join(average_time, group, by=c("ExpGroup"="X"))


#Dataplot for linegraph by time
#geom_line calls the line graph and sets the x axis as Days, y axis as the mean, and the group by ExpGroup
#geom_point puts points in the graph, looks nicer
#geom_errorbar adds the error bars using the sd made earlier
#facet_wrap() seperates the line graphs by STM and OWC (A variable)
ggplot(average_time) +
  geom_line( aes(x=variable, y=mean, group=ExpGroup, color=ExpGroup),linetype=1)+
  geom_point(aes(x=variable, y=mean, group=ExpGroup,color=ExpGroup), shape=15)+
  geom_errorbar(aes(x=variable,ymin=mean-sd, ymax=mean+sd), width=0.1)+
  theme_classic()+ 
  ylab("Cumulative CH4/gram soil")+
  xlab("Day of Enrichment")+
  ggtitle("Trimethyllysine (TML) stimulation of CH4 in 2018 M2 soils") + theme(axis.text.x=element_text(angle=45, hjust=1))+
  facet_wrap(~A)