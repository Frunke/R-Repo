################
###### ggtreeEXTRA manual: https://yulab-smu.top/treedata-book/chapter10.html
################

library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(readxl)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(cowplot)
library(TDbook)
library(treeio)


# read in GTDB trees (from data_processing)
setwd("C:/Users/teaga/Downloads/Grad School/WrightonLab/Presentation/MAG_tree/teagan_MAG_tree" )
arc_tree = read.tree("gtdbtk.ar53.decorated.tree")
arc_tree

# read in taxonomy (from data_processing)
annotation = read.delim("classification_wOutgroup.txt",sep="\t",header = FALSE)
colnames(annotation)=c("MAG","tax")
annotation=annotation%>%separate(tax,into=c("d","p","c","o","f","g","s"),sep=";",remove=FALSE)

#read in data on linking MAG to source
site <- read.delim('site.txt')

#read in data on alling MAGs isolates
mag_for_isolation<- read.delim('isolate.txt')

# make archaeal tree
arc_dat=as.data.frame(arc_tree$tip.label)
colnames(arc_dat)=c("MAG")
arc_dat=left_join(arc_dat,annotation,by="MAG")

# undecorated tree - section 4.2.2 of https://yulab-smu.top/treedata-book/chapter4.html
mgen_tree <- ggtree(arc_tree, color="black", size=0.5, linetype="solid", layout="circular",  branch.length = "branch.length")
mgen_tree 

#color palette for sites - 
Sourcecol <- c("GTDB"="black", "OWC"="sienna1", "STM"="khaki", 'PPR'="springgreen2",
               "PPR_Enrichment"="darkolivegreen4", "OWC_Enrichment"="chocolate4",
               "Rumen"="cadetblue3", "Yojoa"="violetred2")
IsolateCol <- c("no"="white", "yes"="grey33")





#script with dummy data for pangenome analyses
arc=mgen_tree +
  geom_fruit(data=mag_for_isolation,geom=geom_tile,
             mapping=aes(y=MAG,x=0,fill=mag_for_isolation),width = 0.05,offset=0.03,show.legend=TRUE,color = "white",lwd = 0.25,linetype = 1)+
  scale_fill_manual(values=IsolateCol)+
  new_scale_fill()+
  geom_fruit(data=site,geom=geom_tile,
                 mapping=aes(y=MAG,x=0,fill=site),width = 0.05,offset=0.08,show.legend=TRUE,color = "black",lwd = 0.25,linetype = 1)+
  scale_fill_manual(values=Sourcecol)
  
 arc
 
 
 #add taxonomy to tree 
 
 arc <- arc %<+% annotation + new_scale_fill() +
   geom_tippoint(mapping=aes(color= o)) 
 arc
 









