library(staRdom) # folder
library(data.table)
cores <- detectCores(logical = FALSE)

setwd("..")
Samples <- eem_read("Horiba", recursive = TRUE, import_function = "aqualog")

for (i in seq_along(Samples)) {
  Samples[[i]][["x"]]<- Samples[[i]][["x"]][,41:1]
}

#eem_overview_plot(Samples, spp=9, contour = TRUE)
Samples[[1]][["x"]]

# THis takes the name of the folders in the Hitachi folder, then I pasted it as the path of the absorbance files
BiblioDir = list.dirs(path = "Hitachi", full.names = TRUE, recursive = FALSE)
abs=list.files(path=BiblioDir, pattern=".DX")

write.csv(eem_metatemplate(Samples), file="metatable.csv", row.names=T)
meta <- read.table("metatable.csv", header = TRUE, sep = ",", dec = ".", row.names=1)

abs=sub(paste0(meta[,2],"/", meta[,1], ".dx"), pattern="Horiba", replacement="Hitachi")
absorbance <- absorbance_read(abs, package="staRdom") # load csv or txt tables in folder
names(absorbance)[-1] = c(substr(colnames(absorbance[-1]), 1, stop=nchar(colnames(absorbance[-1]))-3)) 

#Auto write a metatable with file locations
write.csv(eem_metatemplate(Samples, absorbance), file="metatable.csv", row.names=T)
meta <- read.table("metatable.csv", header = TRUE, sep = ",", dec = ".", row.names=1)

### Check for problematic samples
problem <- eem_checkdata(eem_list=Samples, absorbance=absorbance, metadata=meta, error=FALSE)

#Absorbance baseline correction
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

eem_list <- eem_extend2largest(Samples, interpolation = 1, extend = FALSE)
eem_raman_area(eem_list, blanks_only = T, average=F)
# There is 1 blank with a very low area, so we discard that
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("S13_Blank_002")
)
eem_list <- eem_exclude(eem_list, exclude)
# Remove the rest of the blankc from the list
eem_list <- eem_remove_blank(eem_list)
eem_overview_plot(eem_list[1:9], spp=9, contour = TRUE)

eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 2)

eem_list<- eem_raman_normalisation2(eem_list, blank = "blank")
eem_overview_plot(eem_list[1:9], spp=9, contour = TRUE)

#remove balnks from the sample list
eem_list <- eem_extract(eem_list, c("blank"),ignore_case = TRUE)
absorbance <- dplyr::select(absorbance, -matches("nano|miliq|milliq|mq|blank", ignore.case = TRUE))

#remove and interpolate scatter
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(40,10,10,10)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)
eem_overview_plot(eem_list[1:9], spp=9, contour = TRUE)

eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)

eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)
eem_overview_plot(eem_list[1], spp=9, contour = TRUE, cores=cores)

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

slope_parms <- abs_parms(absorbance, cuvl = 1, cores = cores)
abs_params=merge(indices_peaks, slope_parms, by="sample")

list(abs_params, abs_params$sample_date)


ggplot(abs_params)+
  facet_grid(~sample_date)+
  geom_line(aes(x=col_no, y=bix, group=replicate))

##### PARAFAC analysis #####
#### Simple models ####

eem_list=eem_red2smallest (eem_list)
#eem_list %>% 
# eem_extract(sample = "", keep = TRUE) %>%
#ggeem(contour = TRUE)

dim_min <- 3
dim_max <- 6

nstart <- 50 # number of similar models from which best is chosen
maxit = 5000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-6 # tolerance in PARAFAC analysis

# calculating PARAFAC models, one for each number of components
pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("uncons", "uncons", "uncons"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf1, contour = TRUE)


# same model but using non-negative constraints
pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf1n, contour = TRUE)

# Checkin the correlatino between components (there shouldn't be any components with high correlation)
eempf_cortable(pf1n[[2]])
eempf_corplot(pf1n[[4]], progress = FALSE, normalisation = FALSE)

# Calculations with normalized sample data
pf2n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
pf2n <- lapply(pf2n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf2n, contour = TRUE)

#check the correlation
#normalization can stay False when chekcing the correlation (see the manual for details)
eempf_corplot(pf2n[[2]], progress = FALSE, normalisation = F)

#### The best model seems to be model 3 (5 component model)

#### Finding and Excluding outliers ####
# calculate leverage
cpl <- eempf_leverage(pf2n[[3]])

# plot leverage (nice plot)
eempf_leverage_plot(cpl, qlabel=0.1)
#exclude <- eempf_leverage_ident(cpl,qlabel=0.1)

# samples, excitation and emission wavelengths to exclude, makes sense after calculation of leverage
exclude <- list("ex" = c(),
                "em" = c(),
                "sample" = c("S14_I_C1", "S13_I_C3")
)

# exclude outliers if neccessary. if so, restart analysis
eem_list_ex <- eem_exclude(eem_list, exclude)

#New model withough the outlilers
pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("uncons", "uncons", "uncons"),maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf3, contour = TRUE)

pf3n <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
pf3n <- lapply(pf3n, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf3n, contour = TRUE)

#Examine Resuduals
# When examining outliers, keep the excluded samples within the loop
eempf_residuals_plot(pf3n[[3]], eem_list, residuals_only = TRUE,  spp = 9, cores = cores, contour = TRUE)

# recaluclating the model with higher accuracy. This takes a lot of time, so try to narrow the moedl comp range as much as possible
####
dim_min <- 3
dim_max <- 6

ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
nstart = 50 # number of random starts
maxit = 10000 # increase number of maximum interations

pf4 <- eem_parafac(eem_list_ex, comps = 5, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, output = "all", cores = cores, strictly_converging = TRUE)
pf4p <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")
eempf_compare(pf4p, contour = TRUE)

eempf_leverage_plot(eempf_leverage(pf4p[[1]])) # [[1]] means the 4th model in the list, 6 component model in that case
eempf_corplot(pf4p[[1]], progress = FALSE)

eempf_comp_load_plot(pf4[[1]], contour = TRUE)
eempf_comps3D(pf4p[[1]])

# Final residuals and sample plots
eempf_residuals_plot(pf4[[1]], eem_list, select = eem_names(eem_list)[10:14], cores = cores, contour = TRUE)

# Split half analysis
split_half <- splithalf(eem_list_ex, comp=4, normalise = TRUE, rand = FALSE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)
split_half2 <- splithalf(eem_list_ex, comp=6, normalise = TRUE, rand = TRUE, cores = cores, nstart = nstart, strictly_converging = TRUE, maxit = maxit, ctol = ctol)

splithalf_plot(split_half)

# 4 component model is better than 5 so better stick to this?
# Tucker’s Congruency Coefficients to see the similiarity of the splits
tcc_sh_table <- splithalf_tcc(split_half)

tcc_sh_table

# Extra model validation
corcondia <- eempf_corcondia(pf4[[1]], eem_list_ex)

eemqual <- eempf_eemqual(pf4[[1]], eem_list_ex, split_half, cores = cores)

# Component importance
varimp <- eempf_varimp(pf4[[1]], eem_list_ex, cores = cores)

