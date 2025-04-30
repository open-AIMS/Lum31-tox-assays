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
         Lum1_30_qtox = (1-(median(CuDat1[1:4,4])- Lum1_30)/median(CuDat1[1:4,4])), 
         Lum2_30_qtox = (1-(median(CuDat1[1:4,6])-Lum2_30)/median(CuDat1[1:4,6])),
         Lum3_30_qtox = (1-(median(CuDat1[1:4,8])- Lum3_30)/median(CuDat1[1:4,8])),
         Microtox1_30_qtox = (1-(median(CuDat1[1:4,10])-Microtox1_30)/median(CuDat1[1:4,10])),
         Microtox2_30_qtox = (1-(median(CuDat1[1:4,12])-Microtox2_30)/median(CuDat1[1:4,12])),
         Microtox3_30_qtox = (1-(median(CuDat1[1:4,14])-Microtox3_30)/median(CuDat1[1:4,14])),
         Lum1_30_qtoxB = Lum1_30_qtox/max(Lum1_30_qtox),
         Lum2_30_qtoxB = Lum2_30_qtox/max(Lum2_30_qtox),
         Lum3_30_qtoxB = Lum3_30_qtox/max(Lum3_30_qtox),
         Microtox1_30_qtoxB = Microtox1_30_qtox/max(Microtox1_30_qtox),
         Microtox2_30_qtoxB = Microtox2_30_qtox/max(Microtox2_30_qtox),
         Microtox3_30_qtoxB = Microtox3_30_qtox/max(Microtox3_30_qtox))

# select only log_x and qtoxB data to use as input for bayesnec

CuDat1_Manip2 <- CuDat1_Manip |> select(16,23:28)

# create dataset to run models at once

dat1 <- lapply(2:ncol(CuDat1_Manip2), FUN = function(i) {
  dat.i <- CuDat1_Manip2[,c("log_x", colnames(CuDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat1) <- colnames(CuDat1_Manip2)[-1]

save(dat1, file = "CuData_30m.RData")


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

load("Lum31MicrotoxCuFits_All_30m.RData")

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

save(Cu_all_plots, file = "plotData_Cu30m.RData")

####---- copper EC50s ----####

CuEC50_30m <- lapply(fits, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(CuEC50_30m, file = "CuEC50s_30m.csv")
save(CuEC50_30m, file = "CuEC50s_30m.RData")

###---- load data: zinc ----####

ZnDat1 <- read_excel("Lum31 Microtox_CuZn_AIMSrepository.xlsx", sheet="Zn")  |> 
  data.frame() 

ZnDat1 <- ZnDat1 |> 
  drop_na(Microtox3_30)

ZnDat1_Manip <- ZnDat1 %>% 
  mutate(x = as.numeric(as.character(x_measuredZn)),
         log_x = log(x),
         Lum1_30_qtox = (1-(median(ZnDat1[1:4,4])- Lum1_30)/median(ZnDat1[1:4,4])), 
         Lum2_30_qtox = (1-(median(ZnDat1[1:4,6])-Lum2_30)/median(ZnDat1[1:4,6])),
         Lum3_30_qtox = (1-(median(ZnDat1[1:4,8])- Lum3_30)/median(ZnDat1[1:4,8])),
         Microtox1_30_qtox = (1-(median(ZnDat1[1:4,10])-Microtox1_30)/median(ZnDat1[1:4,10])),
         Microtox2_30_qtox = (1-(median(ZnDat1[1:4,12])-Microtox2_30)/median(ZnDat1[1:4,12])),
         Microtox3_30_qtox = (1-(median(ZnDat1[1:4,14])-Microtox3_30)/median(ZnDat1[1:4,14])),
         Lum1_30_qtoxB = Lum1_30_qtox/max(Lum1_30_qtox),
         Lum2_30_qtoxB = Lum2_30_qtox/max(Lum2_30_qtox),
         Lum3_30_qtoxB = Lum3_30_qtox/max(Lum3_30_qtox),
         Microtox1_30_qtoxB = Microtox1_30_qtox/max(Microtox1_30_qtox),
         Microtox2_30_qtoxB = Microtox2_30_qtox/max(Microtox2_30_qtox),
         Microtox3_30_qtoxB = Microtox3_30_qtox/max(Microtox3_30_qtox))

# select only log_x and qtoxB data to use as input for bayesnec

ZnDat1_Manip2 <- ZnDat1_Manip |> select(16,23:28)

# create dataset to run models at once

datZn1 <- lapply(2:ncol(ZnDat1_Manip2), FUN = function(i) {
  dat.i <- ZnDat1_Manip2[,c("log_x", colnames(ZnDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(datZn1) <- colnames(ZnDat1_Manip2)[-1]

save(datZn1, file = "ZnData_30m.RData")

####---- Zinc fits HPC ----#####

plan(multisession)
fitsZn <- future_lapply(datZn1, FUN = function(i) {
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsZn) <- names(datZn1)
save(fitsZn, file = "Lum31MicrotoxZnFits_All_30m.RData")

####---- Zinc plots ----####

load("Lum31MicrotoxZnFits_All_30m.RData") 

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

save(Zn_all_plots, file = "plotData_Zn30m.RData")

####---- zinc EC50s ----####

ZnEC50_30m <- lapply(fitsZn, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(ZnEC50_30m, file = "ZnEC50s_30m.csv")
save(ZnEC50_30m, file = "ZnEC50s_30m.RData")


### cycle through model outputs
summary(fitsZn$Lum1_15_qtoxB)
check_chains(fitsZn$Lum1_15_qtoxB)
rhat(fitsZn$Lum1_15_qtoxB, rhat_cutoff = 1.05)


############# Plots: Combined EC50 ############

combECs <- read_excel("EC50s_combined.xlsx", sheet="4R")  |> 
  data.frame() 
ECs15mZn <- combECs |> filter(exposure == "15m", toxicant == "Zn") 
ECs15mCu <- combECs |> filter(exposure == "15m" , toxicant == "Cu") 
ECs30mZn <- combECs |> filter(exposure == "30m" , toxicant == "Zn")
ECs30mCu <- combECs |> filter(exposure == "30m" , toxicant == "Cu")


ggplot(ECs15mZn, aes(x = plateID, y = Q50, color = strain)) +
  geom_point(size = 4) +
  geom_errorbar(aes(x = plateID, ymax = Q97.5, ymin = Q2.5), color = "black") 


Zn15 <- ggplot(ECs15mZn, aes(x = plate_rep, y = Q50, color = strain)) +
  geom_point(size = 5) +
  #scale_color_manual(values = c("blue", "green")) +
  scale_color_manual(values = c("#4da6ff", "#8ECF69")) +
  geom_errorbar(aes(x = plate_rep, ymax = Q97.5, ymin = Q2.5), color = "black", width = 0) +
  scale_y_continuous(limits = c(0, 7)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
Zn15

Zn30 <- ggplot(ECs30mZn, aes(x = plate_rep, y = Q50, color = strain)) +
  geom_point(size = 5) +
  #scale_color_manual(values = c("blue", "green")) +
  scale_color_manual(values = c("#4da6ff", "#8ECF69")) +
  geom_errorbar(aes(x = plate_rep, ymax = Q97.5, ymin = Q2.5), color = "black", width = 0) +
  scale_y_continuous(limits = c(0, 7)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
Zn30

ZnComb <- ggarrange(Zn15, Zn30, ncol = 1, nrow = 2, 
                    common.legend = TRUE, legend = "right") 
ZnComb
  annotate_figure(ZnComb,
                  left = text_grob("Modelled EC50 (mg/L Zinc)", color = "black", rot = 90),
                  bottom = text_grob(""),
                  top = text_grob("Zinc"))

Cu15 <- ggplot(ECs15mCu, aes(x = plate_rep, y = Q50, color = strain)) +
    geom_point(size = 5) +
    #scale_color_manual(values = c("blue", "green")) +
    scale_color_manual(values = c("#4da6ff", "#8ECF69")) +
    geom_errorbar(aes(x = plate_rep, ymax = Q97.5, ymin = Q2.5), color = "black", width = 0) +
    scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
Cu15
  
Cu30 <- ggplot(ECs30mCu, aes(x = plate_rep, y = Q50, color = strain)) +
    geom_point(size = 5) +
    #scale_color_manual(values = c("blue", "green")) +
    scale_color_manual(values = c("#4da6ff", "#8ECF69")) +
    geom_errorbar(aes(x = plate_rep, ymax = Q97.5, ymin = Q2.5), color = "black", width = 0) +
    scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
Cu30

CuComb <- ggarrange(Cu15, Cu30, ncol = 1, nrow = 2, 
                    common.legend = TRUE, legend = "right") 
CuComb
annotate_figure(CuComb,
                left = text_grob("Modelled EC50 (mg/L Copper)", color = "black", rot = 90),
                bottom = text_grob(""),
                top = text_grob("Copper"))

ZnCu <- ggarrange(Cu15, Cu30, Zn15, Zn30, ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "bottom", align = "hv") 
ZnCu
  
library(cowplot)
  
ZnCu_final <- ggdraw(ZnCu) +
  draw_label("Modelled EC50 (mg/L Zinc)", x = 0.01, y = 0.32, angle = 90, size = 12) +
  draw_label("Modelled EC50 (mg/L Copper)", x = 0.01, y = 0.78, angle = 90, size = 12) +
  draw_label("15 min", x = 0.27, y = 1, size = 12, fontface = "bold") +
  draw_label("30 min", x = 0.77, y = 1, size = 12, fontface = "bold") 
ZnCu_final
               
ggsave("Plot_CombinedEC50_v3.pdf", plot = ZnCu_final, device = "pdf", width = 8, height = 6)
ggsave("Plot_CombinedEC50_v3.tiff", plot = ZnCu_final, device = "tiff", width = 8, height = 6)
