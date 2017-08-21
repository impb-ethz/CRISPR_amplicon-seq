library(ggplot2)
library (ggjoy)
library(RColorBrewer)
library(reshape2)
library(ggthemes)
library(cowplot)

data<-read.csv("/Users/devangmehta/Dropbox/PhD/AVCas9_PACBIO/ProteinData/Distributions/distributionAC3_3wpi.csv",header = TRUE)
head(data)

datalong<-melt(data)

head(datalong)
setwd('/Users/devangmehta/Dropbox/PhD/AVCas9_PACBIO/')
tiff("protein_align_AC3_3wpi.tiff",height = 30, width = 20, units = 'cm', compression = "lzw", res = 300)
p<-ggplot(datalong,aes(x=value,y=variable,fill=variable,color=variable))+geom_joy2(scale=2)+scale_fill_manual(values = c("X_WT" = "gray", "X_155" = "red","X_92"="#c6dbef","X_111"="#9ecae1","X_114"="#6baed6","X_118"="#4292c6","X_120"="#2171b5","X_130"="#08519c","X_139"="#08306b"))+scale_color_manual(values=c("grey45","grey45","grey45","grey50","grey80","grey80","grey80","grey80","grey45"))
p+labs(x="Protein length (amino acids)", y="")+annotate("rect",xmin=0,xmax=7, ymin=-Inf,ymax=Inf,fill="firebrick1",alpha=0.4)+annotate("rect",xmin=4,xmax=6, ymin=-Inf,ymax=Inf,fill="firebrick3",alpha=0.4)+theme(axis.title=element_text(size=16, face="bold"),axis.text.x=element_text(size=10),axis.line = element_line(),panel.grid.major.x=element_line(color = "lightgrey"),legend.position="none",panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_blank())
dev.off()