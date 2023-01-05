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
write.csv(data, file="Optical_with_correctedReservoirs.csv")

cols = c("bix", "fi", "hix", "SR", "b", "t", "a", "m", "c", "a254", "a300", "E2_E3", "E4_E6", "S275_295", "S350_400")
data[, paste0(cols, "_median") := lapply(.SD, median, na.rm=T), .SDcols = cols, by=.(sample_date, col_no)]

# PLot individual incies to check
ggplot(data)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=fi, group=replicate))+
  geom_line(aes(x=col_no, y=fi_median, group=sample_date), color="red", lwd=2)

#facet plots for checking them all
data2=data[sample_date%in%c("S08", "S11", "S13","S14","S19")]
cols = c("bix", "fi", "hix", "SR", "a254",  "E2_E3")
data_median=data2[, lapply(.SD, median, na.rm=T), .SDcols = cols, by=c("sample_date", "col_no")]
data_sd=data2[,lapply(.SD, sd, na.rm=T), .SDcols = cols, by=c("sample_date", "col_no")]

data_median_melted=reshape2::melt(data_median, 
                  id.vars = c("sample_date", "col_no"), 
                  measure.vars =cols)
data_sd_melted=reshape2::melt(data_sd, 
                                  id.vars = c("sample_date", "col_no"), 
                                  measure.vars =cols)
data_median_melted$sdmin=data_median_melted$value+data_sd_melted$value
data_median_melted$sdmax=data_median_melted$value-data_sd_melted$value

# Check the median or the individual replicated by uncommenting the plot below
ggplot()+
  facet_wrap(~variable, scale="free")+
  geom_line(data=data_median_melted, aes(x=col_no, y=value, group=sample_date, color=sample_date), lwd=1.2)+
  geom_point(data=data_median_melted, aes(x=col_no, y=value, group=sample_date,fill=sample_date, shape=sample_date), size=3)+
  scale_color_manual(values=c("#5EADD1", "#DB2A2A", "#FFBF1C", "#E0E334", "#001D7A"))+
  scale_fill_manual(values=c("#5EADD1", "#DB2A2A", "#FFBF1C", "#E0E334", "#001D7A"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  scale_shape_manual(values=c(21,22,23,24,25))
  #geom_errorbar(data=data_median_melted,aes(x=col_no, y=value, ymin=sdmin, ymax=sdmax, color=sample_date), width=.2,position=position_dodge(0.05))

# Radar plot to observe changes per column before and after reversal

library(ggradar)
data3=data_median[col_no=="C1"][,c(1,3:8)]
plotting_data=data3[,paste0(colnames(data3[,c(2:7)]), "_scaled"):=lapply(.SD, FUN=scale, scale=T, center=T), .SDcols=cols][,c(1,8:13)]
plotting_data$sample_date=c("Control", "Day 1", "Day3", "Day7", "Day17")
plotting_data$sample_date=factor(plotting_data$sample_date, levels=c("Control", "Day 1", "Day3", "Day7", "Day17"))
col1=ggradar(plotting_data, 
        values.radar = c("-1.8", "0", "1.8"),
        grid.min=-2.4,
        grid.max=2.4, grid.mid = 0, 
        axis.labels = cols, fill=T, fill.alpha=0.1, plot.title="Column 1", legend.position = "none")

data3=data_median[col_no=="C2"][,c(1,3:8)]
plotting_data=data3[,paste0(colnames(data3[,c(2:7)]), "_scaled"):=lapply(.SD, FUN=scale, scale=T, center=T), .SDcols=cols][,c(1,8:13)]
plotting_data$sample_date=c("Control", "Day 1", "Day3", "Day7", "Day17")
plotting_data$sample_date=factor(plotting_data$sample_date, levels=c("Control", "Day 1", "Day3", "Day7", "Day17"))
Col2=ggradar(plotting_data, 
        values.radar = c("-1.8", "0", "1.8"),
        grid.min=-2.4,
        grid.max=2.4, grid.mid = 0, 
        axis.labels = cols, fill=T, fill.alpha=0.1, plot.title="Column 2", legend.position = "none")

data3=data_median[col_no=="C3"][,c(1,3:8)]
plotting_data=data3[,paste0(colnames(data3[,c(2:7)]), "_scaled"):=lapply(.SD, FUN=scale, scale=T, center=T), .SDcols=cols][,c(1,8:13)]
plotting_data$sample_date=c("Control", "Day 1", "Day3", "Day7", "Day17")
plotting_data$sample_date=factor(plotting_data$sample_date, levels=c("Control", "Day 1", "Day3", "Day7", "Day17"))
col3=ggradar(plotting_data, 
        values.radar = c("-1.8", "0", "1.8"),
        grid.min=-2.4,
        grid.max=2.4, grid.mid = 0, 
        axis.labels = cols, fill=T, fill.alpha=0.1, plot.title="Column 3", legend.position = "none")

data3=data_median[,lapply(.SD, FUN= median), .SDcols=cols, by=sample_date]
plotting_data=data3[,paste0(colnames(data3[,c(2:7)]), "_scaled"):=lapply(.SD, FUN=scale, scale=T, center=T), .SDcols=cols][,c(1,8:13)]
plotting_data$sample_date=c("Control", "Day 1", "Day3", "Day7", "Day17")
plotting_data$sample_date=factor(plotting_data$sample_date, levels=c("Control", "Day 1", "Day3", "Day7", "Day17"))
col_all=ggradar(plotting_data, 
        values.radar = c("-1.8", "0", "1.8"),
        grid.min=-2.4,
        grid.max=2.4, grid.mid = 0, 
        axis.labels = cols, fill=T, fill.alpha=0.1, plot.title="Whole Column Set", legend.position = "none")


legend <- cowplot::get_legend(col1)
cowplot::plot_grid(col1,Col2,col3,col_all, legend=legend,   ncol=2, nrow = 3, 
                   layout_matrix = rbind(c(1,2), c(3,4),c(5,5)),
                   widths = c(5, 5), heights = c(0.1, 2.5, 2.5)) 
#### End of optical Indice plots

 #### PRAFAC component plots ####
samples=fread("4comp-NonNormalized.txt", select=c(1:5))
emission=fread("4comp-NonNormalized.txt", select=c(6:10), nrows=64)
excitation=fread("4comp-NonNormalized.txt", select=c(11:15), nrows=41)

# This was no need to copy the steps above, it just replicates the reservoir samples wince they are on the second table already
comp_data=samples[data[,.(sample, sample_date,replicate, col_no)], on = .(sample)]
cols=c("Comp.1", "Comp.2", "Comp.3", "Comp.4")
comp_data[, paste0(cols, "_median") := lapply(.SD, median, na.rm=T), .SDcols = cols, by=.(sample_date, col_no)]
median_values=aggregate(comp_data[,c("Comp.1", "Comp.2", "Comp.3", "Comp.4")], by=comp_data[,c("sample_date", "col_no")], FUN=median, na.rm=T)

ggplot(comp_data)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=Comp.2, group=replicate))+
  geom_line(aes(x=col_no, y=Comp.2_median, group=sample_date), color="red", lwd=2)


# Plot of different components in each column
ggplot(median_values)+
  facet_grid(~col_no, scale="free")+
  geom_line(aes(x=sample_date, y=Comp.1, group=col_no), color="red")+
  geom_line(aes(x=sample_date, y=Comp.2, group=col_no), color="blue3")+
  geom_line(aes(x=sample_date, y=Comp.3, group=col_no), color="black")+
  geom_line(aes(x=sample_date, y=Comp.4, group=col_no), color="green3")+
  geom_point(aes(x=sample_date, y=Comp.1, group=col_no), color="red")+
  geom_point(aes(x=sample_date, y=Comp.2, group=col_no), color="blue3")+
  geom_point(aes(x=sample_date, y=Comp.3, group=col_no), color="black")+
  geom_point(aes(x=sample_date, y=Comp.4, group=col_no), color="green3")+
  geom_vline(xintercept = "S10", color="red", linetype="dashed")

DOC_consumption = fread("C:/Users/c7701233/Nextcloud/Column-Experiment/DOC_measurements/DOC_git/Expmeriment_DOC_dataprep/DOC_consumption.csv", drop="V1")

# Correction as 'C consumed per column'
#Data created from comps and the related corrections

C_correction=rbind(comp_data[DOC_consumption[,.(sample_date, replicate, C1consumed)], C_consumed:=C1consumed, 
                on=.(replicate=replicate, sample_date=sample_date)][col_no=="C1"],
      comp_data[DOC_consumption[,.(sample_date, replicate, C2consumed)], C_consumed:=C2consumed, 
                on=.(replicate=replicate, sample_date=sample_date)][col_no=="C2"],
      comp_data[DOC_consumption[,.(sample_date, replicate, C3consumed)], C_consumed:=C3consumed, 
                on=.(replicate=replicate, sample_date=sample_date)][col_no=="C3"])
C_correction[,c("Comp.1_median","Comp.2_median","Comp.3_median","Comp.4_median"):=NULL]

division=function(Comp,C_consumed){
  result=Comp/C_consumed
  return(as.matrix(result))
}

cols=c("Comp.1", "Comp.2", "Comp.3", "Comp.4")
C_correction[,paste0(cols, "_corrected"):=lapply(.SD, division, C_consumed=C_correction[,.(C_consumed)]), .SDcols=cols]

# Replot all using the corrected parafac components


# This was no need to copy the steps above, it just replicates the reservoir samples wince they are on the second table already
comp_data=C_correction
cols=colnames(comp_data)[endsWith(colnames(comp_data), "corrected")]

comp_data[, paste0(cols, "_median") := lapply(.SD, mean, na.rm=T), .SDcols = cols, by=.(sample_date, col_no)]
median_values=setDT(aggregate(comp_data[, ..cols], by=comp_data[,c("sample_date", "col_no")], FUN=median, na.rm=T))

ggplot(comp_data[!sample_date%in%c("S07", "S06")])+
  facet_grid(~sample_date)+
  #geom_line(aes(x=col_no, y=Comp.2_corrected, group=replicate, color=replicate))+
  geom_line(aes(x=col_no, y=Comp.4_corrected_median, group=sample_date), color="red", lwd=2)


# Plot of different components in each column
ggplot(median_values[!sample_date%in%c("S07", "S06")])+
  facet_wrap(~col_no, scale="free_y")+
  geom_line(aes(x=sample_date, y=Comp.1_corrected, group=col_no), color="red")+
  geom_line(aes(x=sample_date, y=Comp.2_corrected, group=col_no), color="blue3")+
  geom_line(aes(x=sample_date, y=Comp.3_corrected, group=col_no), color="black")+
  geom_line(aes(x=sample_date, y=Comp.4_corrected, group=col_no), color="green3")+
  geom_point(aes(x=sample_date, y=Comp.1_corrected, group=col_no), color="red")+
  geom_point(aes(x=sample_date, y=Comp.2_corrected, group=col_no), color="blue3")+
  geom_point(aes(x=sample_date, y=Comp.3_corrected, group=col_no), color="black")+
  geom_point(aes(x=sample_date, y=Comp.4_corrected, group=col_no), color="green3")+
  geom_vline(xintercept = "S10", color="red", linetype="dashed")

