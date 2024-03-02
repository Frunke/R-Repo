#load in relevant libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(plotly)

df <- read.csv("TML_Enrichment_GC_2024.csv")

#Move the top row to the column names and remove the redundant row
colnames(df) <- df[1, ]
df <- df[-1, ]

#Create a vector that stores "TRUE" if the column header is a CH4 Area (Also the first column)
#Use the vector to make a subset that only has CH4 area
CH4_Columns <- grepl("^CH4 Area", names(df))
CH4_Columns[1] <- TRUE
df <- df[, CH4_Columns]

#Add meaningful headers to the dataframe
colnames(df) <- c("Sample_Names", 0, 6, 15, 22, 29, 36, 43, 50)

#Save a vector with every column but the first
#Use that vector to convert every column to numeric
columns_to_numeric <- names(df)[-1]
df[, columns_to_numeric] <- lapply(df[, columns_to_numeric], as.numeric)

#Add column for soil grams
Soil_grams <- c(3.06, 3.06, 3.02, 3.27, 3.08, 3.36, 3.25, 3.11, 3.06, 3.38, 2.95, 2.99, 2.95, 3.11, 2.97, 3.05)
df$Grams_Soil <- Soil_grams

#Add categories for the samples so I can group them and average.
df$ExpGroup <- c("OWC.Control", "OWC.Control", "OWC.Control", "OWC.Control", 
                         "OWC.TML", "OWC.TML", "OWC.TML", "OWC.TML", "STM.Control", 
                         "STM.Control", "STM.Control", "STM.Control", "STM.TML", "STM.TML", "STM.TML", "STM.TML")
df2 <- df

rownames(df) <- NULL

#Divide the methane area values by the grams of soil
#Mutate function is needed since it otherwise causes NA errors in all values
df <- df %>% 
  mutate(across(2:9, ~ ifelse(Grams_Soil != 0, ./Grams_Soil, 0)))

#Divide the methane area values by Laura's standard curve to get the percentage values
df <- df %>%
  mutate(across(2:9, ~ ifelse(. != 0, ((. + 65.842) / 41035) * 100, 0)))

#Melt the dataframe to a long format by the sample.names
df <- melt(df[c(1:9)],id.vars = "Sample_Names")

#Adds the ExpGroup as a column based on column 1 (Sample_names) and column 12 (ExpGroup) from Methane_subset
#by=c() function means match sample names to sample names to add the expgroup
df <-  inner_join(df2[c(1,11)], df, by=c("Sample_Names"="Sample_Names"))

#joins the group from "group.df" to average_time
group <- read.csv("Group.csv")
df <- inner_join(df, group, by=c("ExpGroup"="X"))

#Convert variable to numeric
df$variable <- as.numeric(as.character(df$variable))

df <- na.omit(df)

df <- df %>% 
  mutate(ExpGroup = case_when(
    ExpGroup == "OWC.Control" ~ "OWC Control",
    TRUE ~ ExpGroup
  ))

df <- df %>% 
  mutate(ExpGroup = case_when(
    ExpGroup == "OWC.TML" ~ "OWC Fed TML Until D36",
    TRUE ~ ExpGroup
  ))

df <- df %>% 
  mutate(ExpGroup = case_when(
    ExpGroup == "STM.TML" ~ "STM TML",
    TRUE ~ ExpGroup
  ))

df <- df %>% 
  mutate(ExpGroup = case_when(
    ExpGroup == "STM.Control" ~ "STM Control",
    TRUE ~ ExpGroup
  ))

df <- df %>% 
      mutate(ExpGroup = ifelse(Sample_Names == "OWC_E08", "OWC Continually Fed TML", ExpGroup))
df <- df %>% 
  mutate(A = ifelse(A == "OWC", "Old Woman Creek", ifelse(A == "STM", "Stordalen Mire", A)))

df$ExpGroup <- factor(df$ExpGroup, levels = c(
  "OWC Control", 
  "OWC Fed TML Until D36", 
  "OWC Continually Fed TML", 
  "STM Control", 
  "STM TML"
))

#Dataplot for linegraph by time
#geom_line calls the line graph and sets the x axis as Days, y axis as the mean, and the group by ExpGroup
#geom_point puts points in the graph, looks nicer
#geom_errorbar adds the error bars using the sd made earlier
#facet_wrap() seperates the line graphs by STM and OWC (A variable)
#Plotly makes the linegraph much smoother
plotly::ggplotly(
  ggplot(df) +
    geom_line( aes(x=variable, y=value, group=Sample_Names, color=ExpGroup),linetype=1)+
    geom_point(aes(x=variable, y=value, group=Sample_Names, color=ExpGroup), shape=15)+
    theme_classic()+ 
    ylab("Percent Methane per Gram of Soil")+
    xlab("Days Since Start of Enrichment")+
    labs(color = "Experimental Group")+
    theme(
      axis.title.x = element_text(size = 20),  
      axis.title.y = element_text(size = 20),
      legend.title = element_text(size = 20), 
      legend.text = element_text(size = 20),
      strip.text = element_text(size = 16),
      plot.margin = margin(b = 0.5, l = 1.5, unit = "cm"))+
    scale_color_manual(values = c("OWC Control" = "#a1d99b", "OWC Fed TML Until D36" = "#31a354",
                                  "OWC Continually Fed TML" = "#1c5e31", "STM Control" = "#9ecae1",
                                  "STM TML" = "#3182bd"))+
    facet_wrap(~A, scales = "free_x")
)


