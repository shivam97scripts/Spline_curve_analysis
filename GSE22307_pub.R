#set path to working directory
#setwd("/Users/shivam/Desktop/experiments/Practice_GSE22307/")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))
install.packages(c("dplyr", "tidyverse", "ggplot2", "ggpubr", "stringr", "gridExtra"))

library(GEOquery)
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(gridExtra)

# data retival 

# matrix retival
gse<-getGEO("GSE22307", GSEMatrix = T)
eSet <- gse[[1]]
raw_data<-exprs(eSet)

normalize_data<-normalizeBetweenArrays(raw_data, method = "quantile")
normalize_data1<-log2(normalize_data+1)

# probe annotation data
feature_data<-fData(eSet)

# change probes name to symbols
probe_ID<-feature_data$ID
geneSymbol<-feature_data$`Gene Symbol`

# Replace rownames of normalized_data with corresponding gene symbols
rownames(normalize_data1) <- geneSymbol[match(rownames(normalize_data1), probe_ID)]

# phenotype data
pheno_data<-pData(eSet)
col_data<-pheno_data[,c("title","time point (day):ch1","geo_accession")]
col_data<-col_data%>%mutate(group = ifelse(`time point (day):ch1`==0, "Untreated", paste0( `time point (day):ch1`, " Days of DSS treatment")))
col_data$day_of_DSS_treatment<-as.numeric(col_data$`time point (day):ch1`)

# OC for quntile normalization process
color_pal<- c("green", "pink","red3","red4")
boxplot(log2(raw_data+1), col= color_pal[as.factor(col_data$day_of_DSS_treatment)]) 
boxplot(log2(normalize_data+1), col= color_pal[as.factor(col_data$day_of_DSS_treatment)])

# calculate Z socre
norm_z_score<- t(scale(t(normalize_data1),center =T, scale =T ))
norm_z_score<-as.data.frame(na.omit(pmax(pmin(norm_z_score, 2), -2)))
norm_z_score$gene<-rownames(norm_z_score)

# long data
data_long<-norm_z_score %>%
  gather(key ='samples', value = 'Normalized_expressoin',-gene)%>%
  left_join(.,col_data, by=c("samples"="geo_accession"))


# Derive gene signature 
GIN_26<- read.csv("dge_update_25jun24.txt", header = F)
GIN26_vals<-str_to_title(GIN_26$V1)

Inflam_expression<-data_long[data_long$gene %in% GIN26_vals,]


# Create spline curve

Inf_dynamics<-function(g1){
  {test_expression<-data_long[data_long$gene == g1,]}
  ggplot() +
  geom_smooth(data = test_expression,
              aes(x = day_of_DSS_treatment, y = Normalized_expressoin),
              formula = y ~ s(x, bs = "cs", k = 4),
              method = "gam",
              color = "darkgreen",
              fill = "lightgreen",
              se = TRUE)+
  geom_smooth(data = Inflam_expression,
              aes(x = day_of_DSS_treatment, y = Normalized_expressoin),
              formula = y ~ s(x, bs = "cs", k = 4),
              method = "gam",
              color = "red",
              fill = "pink",
              se = TRUE)+
  theme_classic()+
  annotate("text",x =6, y = 2.7, label = "Inflammation gene expression", color = "red", hjust = 1)+
  annotate("text",x =6, y = 2.1, label = as.character(substitute(g1)), color = "darkgreen", hjust = 1)
}

#Excecute Inf_dynamics

s1<-Inf_dynamics("Sox9")
s2<-Inf_dynamics("Lgr5")
s3<-Inf_dynamics("S100a9")
s4<-Inf_dynamics("Aldob")
s5<-Inf_dynamics("Muc2")
s6<-Inf_dynamics("Vil1")
s7<-Inf_dynamics("Fabp1")
s8<-Inf_dynamics("Tgfb1")
s9<-Inf_dynamics("Smad3")
s10<-Inf_dynamics("Tgfbr1")
s11<-Inf_dynamics("Tgfbr2")
s12<-Inf_dynamics("Fgf1")
s13<-Inf_dynamics("Fgfr1")

grid.arrange(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
  
  

  


