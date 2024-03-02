library(readxl)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(dplyr)
library(tidyr)


setwd("./16S_Figures/")

###################################
##
## Figuring out sequences per sample
##
###################################
full_data = read.delim("feature_table_wGTDBtax.tsv",sep="\t",header = TRUE)
full_data%>%gather(-taxonomy,-X.OTU.ID,key="sample",value="count")%>%
  filter(sample!="NTC1",sample!="NTC2",sample!="Blank2",sample!="Blank3")%>%
  group_by(sample)%>%summarise(sum=sum(count))%>%
  ggplot(aes(x=sum))+
  geom_histogram()+
  theme_classic()

# plotting as a scatter
full_data%>%gather(-taxonomy,-X.OTU.ID,key="sample",value="count")%>%
  filter(sample!="NTC1",sample!="NTC2",sample!="Blank2",sample!="Blank3")%>%
  group_by(sample)%>%summarise(sum=sum(count))%>%
  ggplot()+
  geom_jitter(aes(x=1,y=sum))+
  scale_y_log10()

full_data%>%gather(-taxonomy,-X.OTU.ID,key="sample",value="count")%>%
  filter(sample!="NTC1",sample!="NTC2",sample!="Blank2",sample!="Blank3")%>%
  group_by(sample)%>%summarise(sum=sum(count))%>%
  arrange(sum)%>%
  ggplot(aes(x=1:nrow(.),y=sum))+
  geom_line()


full_data%>%gather(-taxonomy,-X.OTU.ID,key="sample",value="count")%>%
  filter(sample!="NTC1",sample!="NTC2",sample!="Blank2",sample!="Blank3")%>%
  group_by(sample)%>%summarise(sum=sum(count))%>%
  arrange(sum)
# 14k looks like a good threshold
# what do i lose?
# 1 sample= +CT, not important to Teagan rotation

###################################
##
## Using the feature table rarified to 14k sequences
##
###################################
##read in 16S feature table 
data = read.delim("feature_table_wGTDBtax_14k.tsv", sep="\t",header=TRUE)

## checking goods coverage
data%>%gather(-taxonomy,-X.OTU.ID,key="sample",value="count")%>%
  ungroup()%>%
  group_by(sample)%>%
  summarise(n_seqs=sum(count),
            n_sings=sum(count==1),
            goods=100*(1-n_sings/n_seqs))

data=as.data.frame(data)
rownames(data)=data[,1]
data = data[,-1]

# make data with features as columns and samples as rows
t_data=t(data[,1:15])
log_t_data=log(t_data+1,10)

##NMDS on ASV features
set.seed(13)
NMDS_Bray_data_16S <-metaMDS(t_data, distance = "bray", k=2,
                         autotransform = FALSE, noshare = 0.1, trace = 1,trymax = 100)
bray_dist_16S = metaMDSdist(t_data, distance = "bray", k=2,
                        autotransform = FALSE, noshare = 0.1, trace = 1)
NMDS_Bray_data_16S$stress
# stress = 0.038
ord.data = as.data.frame(vegan::scores(NMDS_Bray_data_16S, display="sites"))

# adding metadata
metadata=read.delim("enrichment_metadata.txt",sep="\t",header=TRUE)
ord.data2=ord.data
ord.data2$sample=row.names(ord.data2)
ord.data2=ord.data2%>%left_join(.,metadata,by=c("sample"="sample.id"))

# bare plot
ggplot(ord.data2)+
  geom_point(aes(x=NMDS1,y=NMDS2)) +
  geom_text(aes(x=NMDS1,y=NMDS2,label=sample))+
theme_classic()


# by Treatment
ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Treatment),size=2)+
  stat_ellipse(aes(color=Treatment))+
  theme_classic()
# significant by treatment? no
set.seed(13)
adonis2(bray_dist_16S ~ ord.data2$Treatment, perm=999) #0.944 no
mrpp(bray_dist_16S, ord.data2$Treatment) # 0.96 no


# by Host
ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=Host),size=2)+
  stat_ellipse(aes(color=Host))+
  theme_classic()
# significant by Host? yes
set.seed(13)
adonis2(bray_dist_16S ~ ord.data2$Host, perm=999) #0.001 yes
mrpp(bray_dist_16S, ord.data2$Host) # 0.001 yes


# by RME
ord.data2%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color=RME),size=2)+
  stat_ellipse(aes(color=RME))+
  theme_classic()
# significant by RME? yes
set.seed(13)
adonis2(bray_dist_16S ~ ord.data2$RME, perm=999) #0.001 yes
mrpp(bray_dist_16S, ord.data2$RME) # 0.001 yes


###################################
##
## Looking at overall community
##
###################################

##read in 16S feature table 
data = read.delim("feature_table_wGTDBtax_14k.tsv", sep="\t",header=TRUE)
data%>%
  gather(-X.OTU.ID,-taxonomy,key="sample",value="count")%>%
  separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%
  filter(p!="Unassigned")%>%
  filter(is.na(c)==0)%>%
  group_by(genus,sample)%>%
  summarise(genus_sum=sum(count)) %>%
  left_join(.,metadata,by=c("sample"="sample.id"))%>%
  ggplot()+
  geom_bar(aes(x=reorder(sample,-tube),y=genus_sum,fill=genus),position="stack",stat="identity",show.legend = FALSE)+
  coord_flip()+
  theme_classic()+
  facet_grid(~Treatment)


# looking at archaeal community
data%>%
  gather(-X.OTU.ID,-taxonomy,key="sample",value="count")%>%
  separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(f,g,sep=";"))%>%
  mutate(species=paste(d,p,c,o,f,g,s,sep=";"))%>%
  filter(d=="d__Archaea")%>%
  filter(p!="p__Thermoproteota")%>%
  group_by(genus,sample)%>%
  summarise(sum=sum(count)) %>%
  left_join(.,metadata,by=c("sample"="sample.id"))%>%
  ggplot()+
  geom_bar(aes(x=reorder(sample,-tube),y=sum/14000,fill=genus),position="stack",stat="identity",show.legend = TRUE)+
  coord_flip()+
  theme_classic()+
  facet_grid(Treatment~Host)



# calculating methanogen count totals
mgen_totals=data%>%
  gather(-X.OTU.ID,-taxonomy,key="sample",value="count")%>%
  separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(f,g,sep=";"))%>%
  mutate(species=paste(d,p,c,o,f,g,s,sep=";"))%>%
  filter(d=="d__Archaea")%>%
  filter(p!="p__Thermoproteota")%>%
  group_by(sample)%>%
  summarise(total=sum(count)) 

# plotting methanogen community abundance
mod_data <- data%>%
  gather(-X.OTU.ID,-taxonomy,key="sample",value="count")%>%
  separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(f,g,sep=";"))%>%
  mutate(species=paste(d,p,c,o,f,g,s,sep=";"))%>%
  filter(d=="d__Archaea")%>%
  filter(p!="p__Thermoproteota")%>%
  group_by(genus,sample)%>%
  summarise(sum=sum(count)) %>%
  left_join(.,mgen_totals,by="sample")%>%
  mutate(abund=sum/total)%>%
  left_join(.,metadata,by=c("sample"="sample.id"))


  ggplot(mod_data2)+
  geom_bar(aes(x=sample,y=abund,fill=genus),position="stack",stat="identity",show.legend = TRUE)+
  theme_minimal()+
  ylab("16S Relative Abundance")+
  xlab("Rumen Reactors")+
  labs(fill = "Taxonomy")+
  scale_fill_manual(values=c("#e0f3db","#e5f5f9","#ccece6","#99d8c9","#66c2a4",
                             "#41ae76","#238b45","#006d2c","#00441b","#67000d",
                             "#fb6a4a","#dd3497","#f768a1","#a6bddb","#74a9cf",
                             "#3690c0","#045a8d","purple","#023858","#fec44f","#8c96c6"))
  
  MAG_data <- mod_data2 %>% filter(sample %in% c("R9", "R10", "R11", "R12"))
  
  ggplot(MAG_data)+
    geom_bar(aes(x=sample,y=abund,fill=genus),position="stack",stat="identity",show.legend = TRUE)+
    theme_minimal()+
    ylab("16S Relative Abundance")+
    xlab("Rumen Reactors")+
    labs(fill = "Taxonomy")+
    scale_fill_manual(values=c("#e0f3db","#e5f5f9","#ccece6","#99d8c9","#66c2a4",
                               "#41ae76","#238b45","#006d2c","#00441b","#67000d",
                               "#fb6a4a","#dd3497","#f768a1","#a6bddb","#74a9cf",
                               "#3690c0","#045a8d","purple","#023858","#fec44f","#8c96c6"))
  
  
  # lookign at ASVs
l=data%>%
  gather(-X.OTU.ID,-taxonomy,key="sample",value="count")%>%
  separate(taxonomy,into=c("d","p","c","o","f","g","s"),sep=";")%>%
  mutate(genus=paste(f,g,sep=";"))%>%
  mutate(species=paste(d,p,c,o,f,g,s,sep=";"))%>%
  filter(g=="g__UBA71")%>%
  mutate(abund=count/14000)%>%
  arrange(-abund)
  
