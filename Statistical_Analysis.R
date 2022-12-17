library(data.table)
library(ggplot2)

samples=fread("4comp-NonNormalized.txt", select=c(1:5))
data=fread("Optical_with_correctedReservoirs.csv", drop="V1")
data$replicate=factor(data$replicate, levels = c("reservoir", "Reservoir","A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
data$col_no=factor(data$col_no, levels = c("Reservoir", "C1", "C2", "C3"))

