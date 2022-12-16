library(data.table)
library(ggplot2)
setwd('..')

#### StaRdom Optical Indice plots####
# Take the data in as data table
data=data.table(read.delim("absorbance_indices.csv", sep=",", row.names = "X"))
replicate_columnNo_Check = data[,.SD][replicate=="Reservoir" | replicate=="reservoir"]
data[grep(x=data$replicate, pattern="reservoir", ignore.case = T)]$col_no="Reservoir"

reservoirs=data[col_no=="Reservoir"]
reservoirs=reservoirs[!sample %in%c("S03_Reservoir_Day1","S15_Reservoir_Day1")] #Delete the same day replicate
date_replicates = unique(data[col_no!=c("Reservoir"), c("replicate", "sample_date")])

#To count how many times to replicate for each sampling date
date_replicates$rep_count=c(1:nrow(date_replicates))
for (i in unique(date_replicates$sample_date)){
  date_replicates[sample_date==i]$rep_count=nrow(date_replicates[sample_date==i])
}

# Create rows that are copies of reservoir of a sampling date with the replicates of all of the col_no of those days
copied_reservoirs=reservoirs[rep(seq_len(nrow(reservoirs)), unique(date_replicates[,.(sample_date,rep_count)])$rep_count)]
all(copied_reservoirs$sample_date %in% date_replicates$sample_date) #To check they are in the same group order
copied_reservoirs$replicate=date_replicates$replicate

#copy the reservoir as a col_no
# Convert the column no, replicate and other variables to factor for easier plotting
data=rbind(data, copied_reservoirs)
data$replicate=factor(data$replicate, levels = c("reservoir", "Reservoir","A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
data$col_no=factor(data$col_no, levels = c("Reservoir", "C1", "C2", "C3"))

cols = c("bix", "fi", "hix", "SR", "b", "t", "a", "m", "c", "a254", "a300", "E2_E3", "E4_E6", "S275_295", "S350_400")
data[, paste0(cols, "_median") := lapply(.SD, median, na.rm=T), .SDcols = cols, by=.(sample_date, col_no)]

# PLot individual incies to check
ggplot(data)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=fi, group=replicate))+
  geom_line(aes(x=col_no, y=fi_median, group=sample_date), color="red", lwd=2)

#facet plots for checking them all
cols = c("bix", "fi", "hix", "SR", "E2_E3") # add any variable you want from the cols list above
data_melted=melt(data, id.vars = c("sample", "sample_date", "replicate","col_no"), measure.vars = cols)
means_melted=melt(data, id.vars = c("sample", "sample_date", "col_no"), measure.vars = paste0(cols, "_median"))

# Check the median or the individual replicated by uncommenting the plot below
ggplot()+
  facet_grid(vars(variable), vars(sample_date), scale="free_y")+
  #geom_line(data=data_melted, aes(x=col_no, y=value, group=replicate))+
  geom_line(data=means_melted, aes(x=col_no, y=value, group=sample_date), color="red", lwd=1.5)
 
#### End of optical Indice plots

 #### PRAFAC component plots ####
samples=fread("4comp-NonNormalized.txt", select=c(1:5))
emission=fread("4comp-NonNormalized.txt", select=c(6:10), nrows=64)
excitation=fread("4comp-NonNormalized.txt", select=c(11:15), nrows=41)

# This was no need to copy the steps above, it just replicates the reservoir samples wince they are on the second table already
comp_data=samples[data[,.(sample, sample_date,replicate, col_no)], on = .(sample)]
cols=c("Comp.1", "Comp.2", "Comp.3", "Comp.4")
comp_data[, paste0(cols, "_median") := lapply(.SD, median, na.rm=T), .SDcols = cols, by=.(sample_date, col_no)]

ggplot(comp_data)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=Comp.4, group=replicate))+
  geom_line(aes(x=col_no, y=Comp.4_median, group=sample_date), color="red", lwd=2)




