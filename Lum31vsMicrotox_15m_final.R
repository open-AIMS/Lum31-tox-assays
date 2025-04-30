###################### Lum31 vs Microtox using measured copper and zinc values #######################

library(tidyverse)
library(readxl)
library(bayesnec)
library(ggpubr)
library(future.apply)

###---- load data: copper ----####

# three replicates for both strains

CuDat1 <- read_excel("Lum31 Microtox_CuZn_AIMSrepository.xlsx", sheet="Cu") %>%
  data.frame() 

CuDat1_Manip <- CuDat1 %>% 
  mutate(x = as.numeric(as.character(x_measuredCu)),
         log_x = log(x),
         Lum1_15_qtox = (1-(median(CuDat1[1:4,3])- Lum1_15)/median(CuDat1[1:4,3])), 
         Lum2_15_qtox = (1-(median(CuDat1[1:4,5])-Lum2_15)/median(CuDat1[1:4,5])),
         Lum3_15_qtox = (1-(median(CuDat1[1:4,7])- Lum3_15)/median(CuDat1[1:4,7])),
         Microtox1_15_qtox = (1-(median(CuDat1[1:4,9])-Microtox1_15)/median(CuDat1[1:4,9])),
         Microtox2_15_qtox = (1-(median(CuDat1[1:4,11])-Microtox2_15)/median(CuDat1[1:4,11])),
         Microtox3_15_qtox = (1-(median(CuDat1[1:4,13])-Microtox3_15)/median(CuDat1[1:4,13])),
         Lum1_15_qtoxB = Lum1_15_qtox/max(Lum1_15_qtox),
         Lum2_15_qtoxB = Lum2_15_qtox/max(Lum2_15_qtox),
         Lum3_15_qtoxB = Lum3_15_qtox/max(Lum3_15_qtox),
         Microtox1_15_qtoxB = Microtox1_15_qtox/max(Microtox1_15_qtox),
         Microtox2_15_qtoxB = Microtox2_15_qtox/max(Microtox2_15_qtox),
         Microtox3_15_qtoxB = Microtox3_15_qtox/max(Microtox3_15_qtox))

# select only log_x and qtoxB data to use as input for bayesnec

CuDat1_Manip2 <- CuDat1_Manip |> select(16,23:28)

# create dataset to run models at once

dat1 <- lapply(2:ncol(CuDat1_Manip2), FUN = function(i) {
  dat.i <- CuDat1_Manip2[,c("log_x", colnames(CuDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat1) <- colnames(CuDat1_Manip2)[-1]

save(dat1, file = "CuData_15m.RData")


######## HPC start ##########

##### need to install future apply as it is not in .sif Yun Kit made with bayesnec
#- module load singularity
#- singularity shell rstudio_bayesnec_20231205_02.sif
#- R
#-install.packages("~/Project/future.apply_1.11.2.tar.gz", repos = NULL)

# R script on HPC called: Lum31mods.R, run.slurm.bayesnec

library(bayesnec)
library(future.apply)

load(file = "CuData_15m.RData")

plan(multisession)
fits <- future_lapply(dat1, FUN = function(i) {
  
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fits) <- names(dat1)
save(fits, file = "Lum31MicrotoxCuFits_All.RData")

######## HPC end ##########

load("Lum31MicrotoxCuFits_All.RData")

####---- copper plots ----####

Cu_all_plots <- lapply(fits, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = FALSE, ecx = TRUE, ecx_val = 50) +
    scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Cu <- ggpubr::ggarrange(
  plotlist = Cu_all_plots, nrow = 2, ncol = 3, labels = names(Cu_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Cu, left = text_grob(expression("Proportion of Luminescence (inversed, scaled)"), rot = 90),
  bottom = text_grob(expression("Copper Concentration (mg/L)"))
)

save(Cu_all_plots, file = "plotData_Cu15m.RData")

####---- copper EC50s ----####

CuEC50_15m <- lapply(fits, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(CuEC50_15m, file = "CuEC50s_15m.csv")
save(CuEC50_15m, file = "CuEC50s_15m.RData")

###---- load data: zinc ----####

ZnDat1 <- read_excel("Lum31 Microtox_CuZn_AIMSrepository.xlsx", sheet="Zn")  |> 
  data.frame() 

ZnDat1 <- ZnDat1 |> 
  drop_na(Microtox3_15)

ZnDat1_Manip <- ZnDat1 %>% 
  mutate(x = as.numeric(as.character(x_measuredZn)),
         log_x = log(x),
         Lum1_15_qtox = (1-(median(ZnDat1[1:4,3])- Lum1_15)/median(ZnDat1[1:4,3])), 
         Lum2_15_qtox = (1-(median(ZnDat1[1:4,5])-Lum2_15)/median(ZnDat1[1:4,5])),
         Lum3_15_qtox = (1-(median(ZnDat1[1:4,7])- Lum3_15)/median(ZnDat1[1:4,7])),
         Microtox1_15_qtox = (1-(median(ZnDat1[1:4,9])-Microtox1_15)/median(ZnDat1[1:4,9])),
         Microtox2_15_qtox = (1-(median(ZnDat1[1:4,11])-Microtox2_15)/median(ZnDat1[1:4,11])),
         Microtox3_15_qtox = (1-(median(ZnDat1[1:4,13])-Microtox3_15)/median(ZnDat1[1:4,13])),
         Lum1_15_qtoxB = Lum1_15_qtox/max(Lum1_15_qtox),
         Lum2_15_qtoxB = Lum2_15_qtox/max(Lum2_15_qtox),
         Lum3_15_qtoxB = Lum3_15_qtox/max(Lum3_15_qtox),
         Microtox1_15_qtoxB = Microtox1_15_qtox/max(Microtox1_15_qtox),
         Microtox2_15_qtoxB = Microtox2_15_qtox/max(Microtox2_15_qtox),
         Microtox3_15_qtoxB = Microtox3_15_qtox/max(Microtox3_15_qtox))

# select only log_x and qtoxB data to use as input for bayesnec

ZnDat1_Manip2 <- ZnDat1_Manip |> select(16,23:28)

# create dataset to run models at once

datZn1 <- lapply(2:ncol(ZnDat1_Manip2), FUN = function(i) {
  dat.i <- ZnDat1_Manip2[,c("log_x", colnames(ZnDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(datZn1) <- colnames(ZnDat1_Manip2)[-1]

save(datZn1, file = "ZnData_15m.RData")

####---- Zinc fits HPC ----#####

plan(multisession)
fitsZn <- future_lapply(datZn1, FUN = function(i) {
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsZn) <- names(datZn1)
save(fitsZn, file = "Lum31MicrotoxCuFits_All.RData")

####---- Zinc plots ----####

load("Lum31MicrotoxZnFits_All.RData")

Zn_all_plots <- lapply(fitsZn, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = FALSE, ecx = TRUE, ecx_val = 50) +
    scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Zn <- ggpubr::ggarrange(
  plotlist = Zn_all_plots, nrow = 2, ncol = 3, labels = names(Zn_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Zn, left = text_grob(expression("Proportion of Luminescence (inversed, scaled)"), rot = 90),
  bottom = text_grob(expression("Zinc Concentration (mg/L)"))
)

save(Zn_all_plots, file = "plotData_Zn15m.RData")

####---- zinc EC50s ----####

ZnEC50_15m <- lapply(fitsZn, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(ZnEC50_15m, file = "ZnEC50s_15m.csv")
save(ZnEC50_15m, file = "ZnEC50s_15m.RData")


### cycle through model outputs
summary(fitsZn$Lum1_15_qtoxB)
check_chains(fitsZn$Lum1_15_qtoxB)
rhat(fitsZn$Lum1_15_qtoxB, rhat_cutoff = 1.05)
