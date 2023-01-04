library(data.table)
library(ggplot2)

samples=fread("4comp-NonNormalized.txt", select=c(1:5))
data=fread("Optical_with_correctedReservoirs.csv", drop="V1")
data$replicate=factor(data$replicate, levels = c("reservoir", "Reservoir","A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
data$col_no=factor(data$col_no, levels = c("Reservoir", "C1", "C2", "C3"))
data=subset(data, subset=!(sample_date>"S10" & replicate=="O"))

#### PCA with the reservoir samples included ####
#Merge the tables for the PCA master table
pca_data=unique(data, by="sample")#Exclude the replicated reservoir values for the pca

pca_data[is.na(pca_data)]<-0
pca_data=pca_data[!sample_date%in%c("S05", "S06")] # Exlude the sample that should not be in the PCA computation
wine.pca <- prcomp(pca_data[,-c("sample", "replicate","sample_date","col_no", "E4_E6")], scale. = TRUE) 
summary(wine.pca)

# PCA plots
# PCA with automated prcomp calcualted rotations #
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)

# Plotting according to the sampling date
ggplot(pca_data, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=pca_data$sample_date, shape=col_no),  size=4)+
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


# Plotting according to the sampling date with before reversal after reversal facets
ggplot(pca_data, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  facet_grid(~sample_date>"S10", scales="free_x", labeller=as_labeller(c('FALSE'="Before Reversal", 'TRUE'="After Reversal")))+
  geom_point(aes(color=pca_data$sample_date, shape=col_no),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), 
               inherit.aes=FALSE, arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  guides(fill="legend")

# Plotting the different columns colored accoording to the column
ggplot(pca_data, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  facet_grid(~col_no, scales="free_x")+
  geom_point(aes(color=sample_date, shape=col_no, fill=sample_date),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  scale_fill_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 18)))+
  scale_shape_manual(values=c(21,22,23,24))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), 
               inherit.aes=FALSE, arrow = arrow(length = unit(1/2, "picas")),color = "black")

# Faceting according to the sampling date

ggplot(pca_results[sample_date%in%c("S08", "S11", "S13", "S14", "S19")], aes(x=PC1, y=PC2))+
  facet_grid(~sample_date, scales="free_y")+
  geom_point(aes(fill=col_no, shape=col_no), alpha=0.7, size=4)+
  scale_fill_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 4)))+
  scale_shape_manual(values=c(21,22,23,24))
# Also try x=PC2 and y =PC1 to see the distance from the reservoir clearly
#####

#### PCA with fewer sampling days with Reservoir ####
pca_data3=unique(data[sample_date%in% c("S08", "S10", "S11", "S13", "S14", "S17", "S19")], by="sample")
pca_data3[is.na(pca_data3)]<-0

wine.pca3 <- prcomp(pca_data3[,!c("sample_date", "col_no", "replicate", "sample", "E4_E6")], scale. = TRUE) 
summary(wine.pca3)

pca_results3=cbind(pca_data3, wine.pca3$x)
# PCA plots
#### PCA with automated prcomp calcualted rotations ####
PCAloadings = data.frame(Variables = rownames(wine.pca3$rotation), wine.pca3$rotation)

# Plotting according to the sampling date
ggplot(data=pca_results3,aes(x=PC1, y=PC2))+
  geom_point(aes(color=sample_date, shape=col_no),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 16)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")

# Faceting according to the sampling date
ggplot(pca_results3, aes(x=PC1, y=PC2))+
  facet_grid(~sample_date, scales="free_y")+
  geom_point(aes(fill=col_no, shape=col_no),  size=4)+
  scale_fill_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 4)))+
  scale_shape_manual(values=c(21,22,23,24))
# Also try x=PC2 and y =PC1 to see the distance from the reservoir clearly

#### PCA with mean values and no Reservoir ####
pca_data2=data[col_no!="Reservoir"]
#pca_data=pca_data[sample_date%in%c("S10","S13", "S16", "S19")]
pca_data2[is.na(pca_data2)]<-0
pca_data2=pca_data2[!sample_date%in%c("S05", "S06")] # Exlude the sample that should not be in the PCA computation

cols=c("bix", "b", "t","a","m","c","fi","hix","a254", "a300","E2_E3","S275_295", "S350_400","S300_700","SR")
pca_data_means=data.table(aggregate(dplyr::select(pca_data2, cols), by=pca_data2[,c("sample_date", "col_no")], FUN=mean, na.rm=T))

wine.pca2 <- prcomp(pca_data_means[,!c("sample_date", "col_no")], scale. = TRUE) 
summary(wine.pca2)

# PCA plots
PCAloadings = data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)

# Plotting according to the sampling date
ggplot(pca_data_means, aes(x=wine.pca2$x[,1], y=wine.pca2$x[,2]))+
  geom_point(aes(color=sample_date, shape=col_no),  size=4)+
  scale_color_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 16)))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")

# mean values for each samping date in order to plot spiders or hulls
pca_results2=cbind(pca_data_means,wine.pca2$x)

# Faceting according to the sampling date
ggplot(pca_results2, aes(x=PC1, y=PC2))+
  facet_grid(~sample_date, scales="free_y")+
  geom_point(aes(fill=sample_date, shape=col_no),  size=4)+
  scale_fill_manual(values=ggthemes::tableau_div_gradient_pal()(seq(0, 1, length = 16)))+
  scale_shape_manual(values=c(21,22,23,24))
# Also try x=PC2 and y =PC1 to see the distance from the reservoir clearly
#### End ####

#### Oxygen preliminary ####
library(vegan)
library(dplyr)
subset_pca=pca_results2[sample_date%in%c("S10", "S13", "S16", "S19")]
plot.new()

centroid_S10=t(summary(ordihull(ord=subset_pca[sample_date=="S10"][,c("PC1","PC2")],  display="species",groups=subset_pca[sample_date=="S10"]$col_no)))[,1:2]
centroid_S14=t(summary(ordihull(ord=subset_pca[sample_date=="S13"][,c("PC1","PC2")],  display="species",groups=subset_pca[sample_date=="S13"]$col_no)))[,1:2]
centroid_S17=t(summary(ordihull(ord=subset_pca[sample_date=="S16"][,c("PC1","PC2")],  display="species",groups=subset_pca[sample_date=="S16"]$col_no)))[,1:2]
centroid_S19=t(summary(ordihull(ord=subset_pca[sample_date=="S19"][,c("PC1","PC2")],  display="species",groups=subset_pca[sample_date=="S19"]$col_no)))[,1:2]

# The change from column 1 to column 3 with direction arrows. This is useful to show S19 approaches S09. But the furthest is the S11
ggplot()+
  geom_segment(data=data.table(centroid_S10), aes(x=PC1, y=PC2,xend=lead(PC1), yend=lead(PC2)), arrow=arrow(), col="red", lwd=2)+
  geom_segment(data=data.table(centroid_S14), aes(x=PC1, y=PC2, xend=lead(PC1), yend=lead(PC2)), col="blue", arrow=arrow(),lwd=2)+
  geom_segment(data=data.table(centroid_S17), aes(x=PC1, y=PC2,xend=lead(PC1), yend=lead(PC2)), arrow=arrow(), col="green", lwd=2)+
  geom_segment(data=data.table(centroid_S19), aes(x=PC1, y=PC2,xend=lead(PC1), yend=lead(PC2)), arrow=arrow(), col="black", lwd=2)

# There might be an oxygen gradient (or any other environmental variable), so check it from here with the smoth contour curves
oxygen_sample_data = fread("Experiment_StaRdom_pipeline/oxygen_sample_data.csv", dec=",")
oxygen_sample_mean = aggregate(oxygen_sample_data[,c("Oxygen", "Temp", "Pressure")] , FUN=median, na.rm=T,
                               by=oxygen_sample_data[,c("Sample_Name", "Column_Number")])

pca_oxygen=pca_results2[oxygen_sample_mean, on=.(sample_date=Sample_Name, col_no=Column_Number)][]

hd_v <- envfit(pca_oxygen[,c("PC1", "PC2")] ~ Oxygen, pca_oxygen)
hd_s <- ordisurf(pca_oxygen[,c("PC1", "PC2")], pca_oxygen$Oxygen, xlim=c(-2,6))
summary(hd_s)
plot(hd_v, col="darkgreen", p.max = 0.1)
ordihull(ord=wine.pca2$x[,c(1:2)],  display="sites", label=T, groups=pca_data_means$sample_date,, show.groups = c("S02", "S08", "S11","S13","S19"))
### End of oxygen preliminary###

#### Statistical Tests ####
# Paired t-test for the optical data
pca_data3

subset_data=pca_data3

dt=subset_data[, .(mean = mean(fi), sd=sd(fi)), by = .(sample_date,col_no)]

ggpubr::ggboxplot(subset_data, x = "sample_date", y = "fi", 
          color = "col_no", palette = c("#00AFBB", "#E7B800", "red", "green"),
          ylab = "bix", xlab = "sample_date")

# Compute the t test after checking normality
shapiro.test(scale(pca_data[sample_date=="S11"]$SR,center=TRUE,scale=TRUE))

hist(log(pca_data$SR), freq=FALSE)
qqnorm(pca_data[,SR],main="Normal Q-Q Plot of male");qqline(pca_data[,SR])

zdata<-scale(pca_data[,-c("sample_date", "col_no", "replicate", "sample", "E4_E6")],center=TRUE,scale=TRUE)
hist(scale(pca_data[sample_date=="S11"]$SR,center=TRUE,scale=TRUE), freq=F)

t.test(subset_data[sample_date=="S08"]$bix, subset_data[sample_date=="S11"]$bix, paired = TRUE)

# a 2-way anova with 2 factors to separate replicate from sample_date
res.aov2=aov(bix~sample_date, data=subset_data[col_no=="C3"])
summary(res.aov2)
TukeyHSD(res.aov2)

# Manova for testing all the variables at once. 
hist(pca_data[m<0.015]$m)
lapply(subset_data[,2:6],  shapiro.test)

# a 2-way anova with 2 factors to separate replicate from sample_date



