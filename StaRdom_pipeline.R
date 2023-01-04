library(staRdom) # folder
library(data.table)
cores <- detectCores(logical = FALSE)
source("functions_FLEELab.R") # alternate funcitons needed with the data from our current machines

setwd("..")


Samples <- eem_read("Horiba", recursive = TRUE, import_function = import_fluoromax4_reverse)

#eem_overview_plot(Samples, spp=9, contour = TRUE)
Samples[[1]][["x"]]

# THis takes the name of the folders in the Hitachi folder, then I pasted it as the path of the absorbance files
BiblioDir = list.dirs(path = "Hitachi", full.names = TRUE, recursive = FALSE)
abs=list.files(path=BiblioDir, pattern=".DX")

#abs=sub(paste0(meta[,2],"/", meta[,1], ".dx"), pattern="Horiba", replacement="Hitachi")
absorbance <- absorbance_Read("Hitachi_all", package="staRdom", dec=".", sep=" ", verbose=T, cores=cores) # load .dx tables in folder with the new absorbance_Read funciton

#Auto write a metatable with file locations
write.csv(eem_metatemplate_df(Samples, absorbance), file="metatable.csv", row.names=T)
meta <- read.table("metatable.csv", header = TRUE, sep = ",", dec = ".", row.names=1)
names(meta)=c("Sample", "eem_location", "abs_location")
### Check for problematic samples
problem <- eem_checkdata(eem_list=Samples, absorbance=absorbance, metadata=meta, error=FALSE)

#Absorbance baseline correction
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

eem_list <- eem_extend2largest(Samples, interpolation = 1, extend = FALSE)
eem_raman_area(eem_list, blanks_only = T, average=F)

# Remove the rest of the blank from the list
eem_list <- eem_remove_blank(eem_list)
eem_overview_plot(eem_list[201:206], spp=25, contour = TRUE)

eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 2, unit="absorbance")
eem_overview_plot(eem_list[201:206], spp=25, contour = TRUE)

eem_list<- eem_raman_normalisation2(eem_list, blank = "blank")
eem_overview_plot(eem_list[5], spp=9, contour = TRUE)

#remove blanks from the sample list
eem_list <- eem_extract(eem_list, c("blank"),ignore_case = TRUE)
absorbance <- dplyr::select(absorbance, -matches("nano|miliq|milliq|mq|blank", ignore.case = TRUE))

#remove and interpolate scatter
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(40,15,15,15)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)
eem_overview_plot(eem_list[3], spp=9, contour = TRUE)

eem_list <- eem_list %>%
  eem_setNA(sample = c(1:13,189,115), ex = 290, interpolate = FALSE) %>%
  eem_setNA(sample = c(240), ex = 280, interpolate = FALSE) %>%
  eem_setNA(sample = c(239), ex = 310, interpolate = FALSE) %>%
  eem_setNA(sample = c(44, 3), ex = 285, interpolate = FALSE) %>%
  eem_setNA(sample = c(132), ex = 275, interpolate = FALSE)%>%
eem_setNA(sample = c(248), ex = c(250,265), interpolate = FALSE) %>%
eem_setNA(sample = c(172), ex = c(260,280), interpolate = FALSE) %>%
eem_setNA(sample = c(162,121), ex = c(300), interpolate = FALSE) %>%
eem_setNA(sample = c(191), ex = c(300,340), interpolate = FALSE) %>%
eem_setNA(sample = c(94,88), ex = c(260), interpolate = FALSE) %>%
eem_setNA(sample = c(84,68), ex = c(270), interpolate = FALSE) %>%
eem_setNA(sample = c(58), ex = c(315), interpolate = FALSE) 

eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)

eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)
eem_overview_plot(eem_list[c(1:13,35,264,210,194,139,110)], spp=9, contour = TRUE)
eem_overview_plot(eem_list[c(3)], spp=9, contour = TRUE)

bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

# Needs to be converted to base or data.table
library(tidyverse)
indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

setDT(indices_peaks)[, c("sample_date", "replicate", "col_no") := tstrsplit(sample, "_")]

slope_parms <- abs_parms(absorbance, cuvl = 2, cores = cores, unit="absorbance")
abs_params=merge(indices_peaks, slope_parms, by="sample")

list(abs_params, abs_params$sample_date)

#For quick plotting of the wanted indices to see some trends
ggplot(abs_params)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=bix, group=replicate))

# Take the optical indice data out
write.csv(abs_params, file="absorbance_indices.csv", row.names=T)

##### PARAFAC analysis #####
#### Simple models ####

eem_list=eem_red2smallest (eem_list)
#eem_list %>% 
# eem_extract(sample = "", keep = TRUE) %>%
#ggeem(contour = TRUE)

dim_min <- 3
dim_max <- 7

nstart <- 50 # number of similar models from which best is chosen
maxit = 5000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-6 # tolerance in PARAFAC analysis

# calculating PARAFAC models, one for each number of components
pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("uncons", "uncons", "uncons"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, verbose=T)
pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf1, contour = TRUE)

# same model but using non-negative constraints
pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, verbose=T)
pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf1n, contour = TRUE)

# Checkin the correlatino between components (there shouldn't be any components with high correlation)
eempf_cortable(pf1n[[3]])
eempf_corplot(pf1n[[2]], progress = FALSE, normalisation = FALSE)

# Calculations with normalized sample data
pf2n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, verbose=T)
pf2n <- lapply(pf2n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf2n, contour = TRUE)

#check the correlation
#normalization can stay False when chekcing the correlation (see the manual for details)
eempf_corplot(pf2n[[4]], progress = FALSE, normalisation = F)

#### The best model seems to be model 3 (5 component model)

#### Finding and Excluding outliers ####
# calculate leverage
cpl <- eempf_leverage(pf2n[[2]])

# plot leverage (nice plot)
eempf_leverage_plot(cpl, qlabel=0.1)
#exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("S14_I_C1", "S13_I_C3", "S01_C_C3")
)

# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eem_list, exclude)

#New model withough the outlilers
pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("uncons", "uncons", "uncons"),maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, verbose=T)
pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf3, contour = TRUE)

pf3n <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores, verbose=T)
pf3n <- lapply(pf3n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf3n, contour = TRUE)

#Examine Resuduals
# When examining outliers, keep the excluded samples within the loop
eempf_corplot(pf3n[[4]], progress = FALSE, normalisation = F)
eempf_residuals_plot(pf3n[[4]], eem_list, residuals_only = TRUE,  spp = 9, cores = cores, contour = TRUE)

cpl <- eempf_leverage(pf5p[[1]])
eempf_leverage_plot(cpl, qlabel=0.1)

# These samples are excluded because they are known to have DOC anomalies
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("S14_I_C1", "S13_I_C3",
                             "S15_Reservoir_Day3", "S03_Reservoir_Day1", 
                             "S06_B_C1", "S06_B_C2", "S06_B_C3","S06_E_C1", "S06_E_C2", "S06_E_C3",
                             "S06_I_C1", "S06_I_C2", "S06_I_C3","S06_P_C1", "S06_P_C2", "S06_P_C3", "S06_Reservoir_Day2",
                             "S07_A_C1", "S07_A_C2", "S07_A_C3","S07_F_C1", "S07_F_C2", "S07_F_C3",
                             "S07_G_C1", "S07_G_C2", "S07_G_C3","S07_M_C1", "S07_M_C2", "S07_M_C3", "S07_Reservoir_Day2")
)
eem_list_ex <- eem_exclude(eem_list, exclude)

# recaluclating the model with higher accuracy. This takes a lot of time, so try to narrow the moedl comp range as much as possible
####
dim_min <- 4
dim_max <- 7

ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
nstart = 50 # number of random starts
maxit = 10000 # increase number of maximum interations

pf4 <- eem_parafac(eem_list_ex, comps = 4, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
pf4p <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf4p, contour = TRUE)

pf5 <- eem_parafac(eem_list_ex, comps = 5, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
pf5p <- lapply(pf5, eempf_rescaleBC, newscale = "Fmax")

pf6 <- eem_parafac(eem_list_ex, comps = 6, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
pf6p <- lapply(pf6, eempf_rescaleBC, newscale = "Fmax")

pf7 <- eem_parafac(eem_list_ex, comps = 7, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
pf7p <- lapply(pf6, eempf_rescaleBC, newscale = "Fmax")

eempf_leverage_plot(eempf_leverage(pf5p[[1]])) # [[1]] means the 4th model in the list, 6 component model in that case
eempf_corplot(pf5p[[1]], progress = FALSE)

eempf_comp_load_plot(pf4p[[1]], contour = TRUE)
eempf_residuals_plot(pf6p[[1]], eem_list_ex, residuals_only = TRUE,  spp = 9, cores = cores, contour = TRUE)

# Final residuals and sample plots
eempf_residuals_plot(pf6[[1]], eem_list, select = eem_names(eem_list)[35:39], cores = cores, contour = TRUE)

# Split half analysis
split_half <- splithalf(eem_list_ex, comp=5, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)
split_half2 <- splithalf(eem_list_ex, comp=4, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)
split_half3 <- splithalf(eem_list_ex, comp=6, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)
split_half4 <- splithalf(eem_list_ex, comp=7, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)


splithalf_plot(split_half3)

# 4 component model is better than 5 so better stick to this?
# Tucker’s Congruency Coefficients to see the similiarity of the splits
tcc_sh_table <- splithalf_tcc(split_half3)

tcc_sh_table

# Extra model validation
corcondia <- eempf_corcondia(pf4[[1]], eem_list_ex)

eemqual <- eempf_eemqual(pf4[[1]], eem_list_ex, split_half, cores = cores)

# Component importance
varimp <- eempf_varimp(pf4[[1]], eem_list_ex, cores = cores)

## Export the normalized components, denormalized components and sample loadings
#Openfluor export
eempf_openfluor(pf4[[1]], file = "4comp_openfluor-NonNormalized.txt")
eempf_openfluor(pf4p[[1]], file = "4comp_openfluor-Normalized.txt", Fmax=TRUE)
eempf_report(pf4[[1]], export = "parafac_report.html", eem_list = eem_list_ex, shmodel = splithalf, performance = TRUE)
eempf_export(pf4[[1]], export="4comp-NonNormalized.txt")
eempf_export(pf4p[[1]], export="4comp-Normalized.txt", Fmax=T)

