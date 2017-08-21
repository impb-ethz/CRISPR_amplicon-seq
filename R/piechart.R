library(ggplot2)
library(RColorBrewer)
library(readxl)

data<-read_xlsx("/Users/devangmehta/Dropbox/PhD/AVCas9_PACBIO/ProteinData/Distributions/Pie_3wpi.xlsx", col_names = TRUE)
head(data)

setwd('/Users/devangmehta/Dropbox/PhD/AVCas9_PACBIO/Pie/3wpi/')
tiff("AC3_130.tiff",height = 10, width = 10, units = 'cm', compression = "lzw", res = 300)

p<-ggplot(data,aes(x=1,y=Percent,fill=Type))+geom_bar(stat="identity", colour='black')+scale_fill_manual(values=c("#d8b365","#f5f5f5","#5ab4ac"))
p<-p+coord_polar(theta='y')
y.breaks <- cumsum(data$Percent)-data$Percent/2
y.breaks

p <- p +theme_minimal()+
  theme(axis.ticks=element_blank(),  # the axis ticks
        axis.title=element_blank(),  # the axis labels
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        legend.key = element_blank(),
        legend.position="none")

p

dev.off()
p