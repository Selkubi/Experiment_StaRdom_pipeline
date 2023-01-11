library(data.table)
library(ggplot2)

data3=fread("Optical_with_correctedReservoirs.csv", drop="V1", stringsAsFactors = T)
cols=c("sample", "sample_date")
data3[,paste0(cols):=lapply(.SD, FUN=as.character), .SDcols=cols]
data3=subset(data3, subset=!(sample_date>"S10" & replicate=="O"))

cell_counts= read.csv("C:/Users/c7701233/Nextcloud/Column-Experiment/EEA/Cell_Counts_standardized.csv")
pca_data3=unique(data3[sample_date%in% c("S08", "S10", "S11", "S12", "S13", "S14", "S17",  "S16", "S19")], by="sample")

reps=as.character(unique(pca_data3[sample_date=="S11"]$replicate)) #takes the replicates of that sampling date
anova_first=pca_data3[sample_date%in% c("S08", "S11") & replicate %in% reps] # copies only the matching replicates to the data frame

shapiro.test(anova_first[col_no=="C1"]$SR)
aov=aov(fi~sample_date*col_no, data=anova_first)
summary(aov)
TukeyHSD(aov)

anova_first=anova_first[replicate!="Reservoir"]
pairwise.t.test(x=anova_first$fi, g=interaction(anova_first$sample_date, anova_first$col_no), 
                p.adjust.method = "bonferroni", paired = F, pool.sd = T)

# Anova uing only matching replicates
anova_second=pca_data3[sample_date%in% c("S10", "S13", "S16", "S19") & replicate!="Reservoir"] # copies only the matching replicates to the data frame
anova_second=anova_second[!(sample_date=="S13" & replicate%in%c("A","B", "C", "F", "G", "I", "P"))]

aov2=aov(E2_E3~sample_date*col_no, data=anova_second)
summary(aov2)
TukeyHSD(aov2)
