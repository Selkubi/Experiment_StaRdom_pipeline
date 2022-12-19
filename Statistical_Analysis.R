library(data.table)
library(ggplot2)

samples=fread("4comp-NonNormalized.txt", select=c(1:5))
data=fread("Optical_with_correctedReservoirs.csv", drop="V1")
data$replicate=factor(data$replicate, levels = c("reservoir", "Reservoir","A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
data$col_no=factor(data$col_no, levels = c("Reservoir", "C1", "C2", "C3"))

#Merge the tables for the PCA master table
pca_data=unique(data, by="sample")#Exclude the replicated reservoir values for the pca

pca_data[is.na(pca_data)]<-0
wine.pca <- prcomp(pca_data[,-c("sample", "replicate","sample_date","col_no", "E4_E6")], scale. = TRUE) 
summary(wine.pca)

# PCA plots
#### PCA with automated prcomp calcualted rotations ####
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)

# Plotting according to the sampling date
ggplot(pca_data, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=pca_data$sample_date),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")

# mean values for each samping date in order to plot spiders or hulls
pca_results=cbind(pca_data,wine.pca$x)
pca_results[, c("sample_date2", "replicate", "col_no2"):=tstrsplit(sample, "_")]
cols=c("PC1", "PC2")
pca_results[, paste0(cols, "_median") := lapply(.SD, median, na.rm=T), .SDcols = cols, by=.(sample_date)]

ggplot(pca_results, aes(x=pca_results$PC1_median, y=pca_results$PC2_median))+
  geom_point(aes(color=pca_data$sample_date),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")

# Plotting according to the sampling date with before reversal after reversal facets
ggplot(pca_data, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  facet_grid(~sample_date>"S10", scales="free_x", labeller=as_labeller(c('FALSE'="Before Reversal", 'TRUE'="After Reversal")))+
  geom_point(aes(color=pca_data$sample_date, shape=col_no),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), 
               inherit.aes=FALSE, arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", aes(x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4), inherit.aes=FALSE),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")
