library(data.table)
library(ggplot2)

data3=fread("Optical_with_correctedReservoirs.csv", drop="V1", stringsAsFactors = T)
cols=c("sample", "sample_date")
data3[,paste0(cols):=lapply(.SD, FUN=as.character), .SDcols=cols]
data3=subset(data3, subset=!(sample_date>"S10" & replicate=="O"))

pca_data3=unique(data3[sample_date%in% c("S08",  "S11", "S13", "S14",  "S19")], by="sample")

# Testing variance homogeneity
car::leveneTest(fi ~ sample_date*col_no, pca_data3)

# Testing for normality
mvnormtest::mshapiro.test(t(cbind(pca_data3$bix,pca_data3$hix,pca_data3$fi)))

# Though this can be ignored due to large number of observations
mov_additive=manova(cbind(bix,hix,fi)~sample_date+col_no, data=pca_data3)
mov_interaction=manova(cbind(bix,hix,fi)~sample_date*col_no, data=pca_data3)
mov_interaction_error=manova(cbind(bix,hix,fi)~sample_date*col_no+replicate, data=pca_data3)

summary(mov_additive)
summary(mov_interaction)
summary(mov_interaction_error)

summary.aov(mov_additive)

library(AICcmodavg)
model.set=list(mov_additive, mov_interaction, mov_interaction_error)
model.names=c("mov_additive", "mov_interaction", "mov_interaction_error")

AIC(mov_additive)
