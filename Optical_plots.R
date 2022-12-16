library(data.table)
library(ggplot2)
setwd('..')

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

ggplot(data)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=bix, group=replicate))
# hatayi duzeltmek icin replicate daki reservoir i 4 kre cogaltman lazim