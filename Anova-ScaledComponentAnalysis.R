library(data.table)
library(ggplot2)
library(MetStaT)

data3=fread("Optical_with_correctedReservoirs.csv", drop="V1")
data3$replicate=factor(data3$replicate, levels = c( "Reservoir","A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
data3$col_no=factor(data3$col_no, levels = c("Reservoir", "C1", "C2", "C3"))
data3=subset(data3, subset=!(sample_date>"S10" & replicate=="O"))

pca_data3=unique(data3[sample_date%in% c("S08",  "S11", "S13", "S14",  "S19")], by="sample")

ASCA.Calculate(data=pca_data3, levels=(pca_data3[,c("sample_date", "replicate", "col_no")]))

pca_data3$replicate=as.numeric(pca_data3$replicate)
pca_data3$col_no=as.numeric(pca_data3$col_no)
pca_data3$sample_date=factor(pca_data3$sample_date)
pca_data3$sample_date=as.numeric(pca_data3$sample_date)
pca_data3[is.na(pca_data3)]<-0
pca_data3=data.frame(pca_data3)
df1=as.matrix(pca_data3[,-c(1,10,11,12,16)])
df2=as.matrix(pca_data3[,c("sample_date", "replicate", "col_no")])

pca=ASCA.Calculate(data=df1, levels=df2)
pca=ASCA.Calculate(data=as.matrix(pca_data3[,-c(1,10,11,12,16)]), 
                   levels=as.matrix(pca_data3[,c("sample_date", "replicate", "col_no")]),
                   equation.elements="1,2,3,12,13")

ASCA.Plot(pca)

