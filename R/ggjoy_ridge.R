library(ggplot2)
library (ggjoy)
library(RColorBrewer)
library(reshape2)
library(ggthemes)
data<-read.csv("/Users/lbc65--lbc65.csv",header = TRUE)
datalong<-melt(data=data,
               id.vars=c("height","sequences"),
               variable.name="Nucleotides",
               value.name = "Reads"
               )

head(datalong)

setwd('/Users/PairwiseAlignImages/whole')
tiff("155_3.tiff",height = 30, width = 12, units = 'cm', compression = "lzw", res = 300)

p<-ggplot(data=datalong,mapping=aes(x=Nucleotides,y=Reads*3,height=rev(height),group=height,color=rev(height)))+scale_color_continuous_tableau("Red")
p+geom_ridgeline(scale=2,alpha="0")+scale_x_discrete(breaks = c("N10","N20","N30","N40","N50","N60","N70","N80"),expand=c(0, 0))+annotate("rect",xmin="N42",xmax="N62", ymin=-Inf,ymax=Inf,fill="firebrick1",alpha=0.4)+annotate("rect",xmin="N45",xmax="N46", ymin=-Inf,ymax=Inf,fill="firebrick3",alpha=0.4)+theme(axis.title=element_text(size=16, face="bold"),axis.text.x=element_text(size=10),axis.line = element_line(),panel.grid.major.x=element_line(color = "lightgrey"),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none",panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_blank())

dev.off()
