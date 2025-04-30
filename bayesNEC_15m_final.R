
library(tidyverse)
library(readxl)
library(bayesnec)
library(ggpubr)
#library(future.apply)

###---- load data: copper ----####

# five lots of data corresponding to different chemistry sampling

CuDat1 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="19Mar24") %>%
  data.frame() 

CuDat1_Manip <- CuDat1 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y15_qtox = (1-(median(CuDat1[1:4,3])-Rep1B_y15)/median(CuDat1[1:4,3])), 
         Rep2B_y15_qtox = (1-(median(CuDat1[1:4,6])-Rep2B_y15)/median(CuDat1[1:4,6])),
         Rep3A_y15_qtox = (1-(median(CuDat1[1:4,9])-Rep3A_y15)/median(CuDat1[1:4,9])),
         Rep4A_y15_qtox = (1-(median(CuDat1[1:4,12])-Rep4A_y15)/median(CuDat1[1:4,12])),
         Rep1B_y15_qtoxB = Rep1B_y15_qtox/max(Rep1B_y15_qtox),
         Rep2B_y15_qtoxB = Rep2B_y15_qtox/max(Rep2B_y15_qtox),
         Rep3A_y15_qtoxB = Rep3A_y15_qtox/max(Rep3A_y15_qtox),
         Rep4A_y15_qtoxB = Rep4A_y15_qtox/max(Rep4A_y15_qtox))

CuDat2 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="28Mar24") %>%
  data.frame() 

CuDat2_Manip <- CuDat2 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y15_qtox = (1-(median(CuDat2[1:4,3])-Rep1B_y15)/median(CuDat2[1:4,3])), 
         Rep2B_y15_qtox = (1-(median(CuDat2[1:4,6])-Rep2B_y15)/median(CuDat2[1:4,6])),
         Rep1B_y15_qtoxB = Rep1B_y15_qtox/max(Rep1B_y15_qtox),
         Rep2B_y15_qtoxB = Rep2B_y15_qtox/max(Rep2B_y15_qtox))
        
CuDat3 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="5Apr24") %>%
  data.frame() 

CuDat3_Manip <- CuDat3 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep3B_y15_qtox = (1-(median(CuDat3[1:4,3])-Rep3B_y15)/median(CuDat3[1:4,3])), 
         Rep4B_y15_qtox = (1-(median(CuDat3[1:4,6])-Rep4B_y15)/median(CuDat3[1:4,6])),
         Rep5B_y15_qtox = (1-(median(CuDat3[1:4,9])-Rep5B_y15)/median(CuDat3[1:4,9])),
         Rep3B_y15_qtoxB = Rep3B_y15_qtox/max(Rep3B_y15_qtox),
         Rep4B_y15_qtoxB = Rep4B_y15_qtox/max(Rep4B_y15_qtox),
         Rep5B_y15_qtoxB = Rep5B_y15_qtox/max(Rep5B_y15_qtox))

CuDat4 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="23Apr24") %>%
  data.frame() 

CuDat4_Manip <- CuDat4 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep6B_y15_qtox = (1-(median(CuDat4[1:4,3])-Rep6B_y15)/median(CuDat4[1:4,3])), 
         Rep7B_y15_qtox = (1-(median(CuDat4[1:4,6])-Rep7B_y15)/median(CuDat4[1:4,6])),
         Rep8B_y15_qtox = (1-(median(CuDat4[1:4,9])-Rep8B_y15)/median(CuDat4[1:4,9])),
         Rep9B_y15_qtox = (1-(median(CuDat4[1:4,12])-Rep9B_y15)/median(CuDat4[1:4,12])),
         Rep10B_y15_qtox = (1-(median(CuDat4[1:4,15])-Rep10B_y15)/median(CuDat4[1:4,15])),
         Rep6B_y15_qtoxB = Rep6B_y15_qtox/max(Rep6B_y15_qtox),
         Rep7B_y15_qtoxB = Rep7B_y15_qtox/max(Rep7B_y15_qtox),
         Rep8B_y15_qtoxB = Rep8B_y15_qtox/max(Rep8B_y15_qtox),
         Rep9B_y15_qtoxB = Rep9B_y15_qtox/max(Rep9B_y15_qtox),
         Rep10B_y15_qtoxB = Rep10B_y15_qtox/max(Rep10B_y15_qtox))

CuDat5 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="1Oct24") %>%
  data.frame() 

CuDat5_Manip <- CuDat5 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Batch6_A1_y15_qtox = (1-(median(CuDat5[1:4,3])-Batch6_A1_y15)/median(CuDat5[1:4,3])), 
         Batch6_A2_y15_qtox = (1-(median(CuDat5[1:4,6])-Batch6_A2_y15)/median(CuDat5[1:4,6])),
         Batch6_A1_y15_qtoxB = Batch6_A1_y15_qtox/max(Batch6_A1_y15_qtox),
         Batch6_A2_y15_qtoxB = Batch6_A2_y15_qtox/max(Batch6_A2_y15_qtox))
         
# select only log_x and qtoxB data to use as input for bayesnec

CuDat1_Manip2 <- CuDat1_Manip |> select(14,19:22)
CuDat2_Manip2 <- CuDat2_Manip |> select(8,11:12)
CuDat3_Manip2 <- CuDat3_Manip |> select(11,15:17)
CuDat4_Manip2 <- CuDat4_Manip |> select(17,23:27)
CuDat5_Manip2 <- CuDat5_Manip |> select(8,11:12)

# create datasets to run models at once

dat1 <- lapply(2:ncol(CuDat1_Manip2), FUN = function(i) {
  dat.i <- CuDat1_Manip2[,c("log_x", colnames(CuDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat1) <- colnames(CuDat1_Manip2)[-1]

dat2 <- lapply(2:ncol(CuDat2_Manip2), FUN = function(i) {
  dat.i <- CuDat2_Manip2[,c("log_x", colnames(CuDat2_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat2) <- colnames(CuDat2_Manip2)[-1]

dat3 <- lapply(2:ncol(CuDat3_Manip2), FUN = function(i) {
  dat.i <- CuDat3_Manip2[,c("log_x", colnames(CuDat3_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat3) <- colnames(CuDat3_Manip2)[-1]

dat4 <- lapply(2:ncol(CuDat4_Manip2), FUN = function(i) {
  dat.i <- CuDat4_Manip2[,c("log_x", colnames(CuDat4_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat4) <- colnames(CuDat4_Manip2)[-1]

dat5 <- lapply(2:ncol(CuDat5_Manip2), FUN = function(i) {
  dat.i <- CuDat5_Manip2[,c("log_x", colnames(CuDat5_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat5) <- colnames(CuDat5_Manip2)[-1]

datCu_list <- c(dat1, dat2, dat3, dat4, dat5)

#save(datCu_list, file = "CuData_all.RData")


######## HPC start ##########

library(bayesnec)
library(future.apply)

load(file = "CuData_all.RData")

plan(multisession)
fitsCu <- future_lapply(datCu_list, FUN = function(i) {

  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsCu) <- names(datCu_list)
save(fitsCu, file = "CuFits_All_15m.RData")

######## HPC end ##########

load("CuFits_All_15m.RData")

#---- copper plots ----#

Cu_all_plots <- lapply(fitsCu, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = FALSE, ecx = TRUE, ecx_val = 50) +
    #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Cu <- ggpubr::ggarrange(
  plotlist = Cu_all_plots, nrow = 5, ncol = 4, labels = names(Cu_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Cu, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Copper concentration (mg/L)"))
)

####---- rename plate replicates ----####

names(Cu_all_plots) <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5", "Rep6", "Rep7", "Rep8", "Rep9", "Rep10",
                         "Rep11", "Rep12", "Rep13", "Rep14", "Rep15", "Rep16")

save(Cu_all_plots, file = "PlotData_Cu_15m_log10.RData")


  
####---- copper EC50s ----####

CuEC50_15m <- lapply(fitsCu, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(CuEC50_15m, file = "CuEC50s_15m.csv")
save(CuEC50_15m, file = "CuEC50s_15m.RData")

       
###---- load data: zinc ----####

ZnDat1 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="19Mar24") %>%
  data.frame() 

ZnDat1_Manip <- ZnDat1 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y15_qtox = (1-(median(ZnDat1[1:4,3])-Rep1B_y15)/median(ZnDat1[1:4,3])), 
         Rep2B_y15_qtox = (1-(median(ZnDat1[1:4,6])-Rep2B_y15)/median(ZnDat1[1:4,6])),
         Rep3A_y15_qtox = (1-(median(ZnDat1[1:4,9])-Rep3A_y15)/median(ZnDat1[1:4,9])),
         Rep4A_y15_qtox = (1-(median(ZnDat1[1:4,12])-Rep4A_y15)/median(ZnDat1[1:4,12])),
         Rep1B_y15_qtoxB = Rep1B_y15_qtox/max(Rep1B_y15_qtox),
         Rep2B_y15_qtoxB = Rep2B_y15_qtox/max(Rep2B_y15_qtox),
         Rep3A_y15_qtoxB = Rep3A_y15_qtox/max(Rep3A_y15_qtox),
         Rep4A_y15_qtoxB = Rep4A_y15_qtox/max(Rep4A_y15_qtox))

ZnDat2 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="28Mar24") %>%
  data.frame() 

ZnDat2_Manip <- ZnDat2 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep5B_y15_qtox = (1-(median(ZnDat2[1:4,3])-Rep5B_y15)/median(ZnDat2[1:4,3])),
         Rep1B_y15_qtox = (1-(median(ZnDat2[1:4,6])-Rep1B_y15)/median(ZnDat2[1:4,6])), 
         Rep2B_y15_qtox = (1-(median(ZnDat2[1:4,9])-Rep2B_y15)/median(ZnDat2[1:4,9])),
         Rep5B_y15_qtoxB = Rep5B_y15_qtox/max(Rep5B_y15_qtox),
         Rep1B_y15_qtoxB = Rep1B_y15_qtox/max(Rep1B_y15_qtox),
         Rep2B_y15_qtoxB = Rep2B_y15_qtox/max(Rep2B_y15_qtox))

ZnDat3 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="5Apr24") %>%
  data.frame() 

ZnDat3_Manip <- ZnDat3 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep3B_y15_qtox = (1-(median(ZnDat3[1:4,3])-Rep3B_y15)/median(ZnDat3[1:4,3])), 
         Rep4B_y15_qtox = (1-(median(ZnDat3[1:4,6])-Rep4B_y15)/median(ZnDat3[1:4,6])),
         Rep5B_y15_qtox = (1-(median(ZnDat3[1:4,9])-Rep5B_y15)/median(ZnDat3[1:4,9])),
         Rep3B_y15_qtoxB = Rep3B_y15_qtox/max(Rep3B_y15_qtox),
         Rep4B_y15_qtoxB = Rep4B_y15_qtox/max(Rep4B_y15_qtox),
         Rep5B_y15_qtoxB = Rep5B_y15_qtox/max(Rep5B_y15_qtox))

ZnDat4 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="23Apr24") %>%
  data.frame() 

ZnDat4_Manip <- ZnDat4 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep6B_y15_qtox = (1-(median(ZnDat4[1:4,3])-Rep6B_y15)/median(ZnDat4[1:4,3])), 
         Rep7B_y15_qtox = (1-(median(ZnDat4[1:4,6])-Rep7B_y15)/median(ZnDat4[1:4,6])),
         Rep8B_y15_qtox = (1-(median(ZnDat4[1:4,9])-Rep8B_y15)/median(ZnDat4[1:4,9])),
         Rep9B_y15_qtox = (1-(median(ZnDat4[1:4,12])-Rep9B_y15)/median(ZnDat4[1:4,12])),
         Rep10B_y15_qtox = (1-(median(ZnDat4[1:4,15])-Rep10B_y15)/median(ZnDat4[1:4,15])),
         Rep6B_y15_qtoxB = Rep6B_y15_qtox/max(Rep6B_y15_qtox),
         Rep7B_y15_qtoxB = Rep7B_y15_qtox/max(Rep7B_y15_qtox),
         Rep8B_y15_qtoxB = Rep8B_y15_qtox/max(Rep8B_y15_qtox),
         Rep9B_y15_qtoxB = Rep9B_y15_qtox/max(Rep9B_y15_qtox),
         Rep10B_y15_qtoxB = Rep10B_y15_qtox/max(Rep10B_y15_qtox))

ZnDat5 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="1Oct24") %>%
  data.frame() 

ZnDat5_Manip <- ZnDat5 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Batch6_A1_y15_qtox = (1-(median(ZnDat5[1:4,3])-Batch6_A1_y15)/median(ZnDat5[1:4,3])), 
         Batch6_A2_y15_qtox = (1-(median(ZnDat5[1:4,6])-Batch6_A2_y15)/median(ZnDat5[1:4,6])),
         Batch6_A1_y15_qtoxB = Batch6_A1_y15_qtox/max(Batch6_A1_y15_qtox),
         Batch6_A2_y15_qtoxB = Batch6_A2_y15_qtox/max(Batch6_A2_y15_qtox))
        


# select only log_x and qtoxB data to use as input for bayesnec

ZnDat1_Manip2 <- ZnDat1_Manip |> select(14,19:22)
ZnDat2_Manip2 <- ZnDat2_Manip |> select(11,15:17)
ZnDat3_Manip2 <- ZnDat3_Manip |> select(11,15:17)
ZnDat4_Manip2 <- ZnDat4_Manip |> select(17,23:27)
ZnDat5_Manip2 <- ZnDat5_Manip |> select(8,11:12)

dat1_Zn <- lapply(2:ncol(ZnDat1_Manip2), FUN = function(i) {
  dat.i <- ZnDat1_Manip2[,c("log_x", colnames(ZnDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat1_Zn) <- colnames(ZnDat1_Manip2)[-1]

dat2_Zn <- lapply(2:ncol(ZnDat2_Manip2), FUN = function(i) {
  dat.i <- ZnDat2_Manip2[,c("log_x", colnames(ZnDat2_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat2_Zn) <- colnames(ZnDat2_Manip2)[-1]

dat3_Zn <- lapply(2:ncol(ZnDat3_Manip2), FUN = function(i) {
  dat.i <- ZnDat3_Manip2[,c("log_x", colnames(ZnDat3_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat3_Zn) <- colnames(ZnDat3_Manip2)[-1]

dat4_Zn <- lapply(2:ncol(ZnDat4_Manip2), FUN = function(i) {
  dat.i <- ZnDat4_Manip2[,c("log_x", colnames(ZnDat4_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat4_Zn) <- colnames(ZnDat4_Manip2)[-1]

dat5_Zn <- lapply(2:ncol(ZnDat5_Manip2), FUN = function(i) {
  dat.i <- ZnDat5_Manip2[,c("log_x", colnames(ZnDat5_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat5_Zn) <- colnames(ZnDat5_Manip2)[-1]

dat_list_Zn <- c(dat1_Zn, dat2_Zn, dat3_Zn, dat4_Zn, dat5_Zn)

save(dat_list_Zn, file = "ZnData_all.RData")

####---- HPC Zinc fits ----#####

plan(multisession)
fitsZn <- future_lapply(dat_list_Zn, FUN = function(i) {
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsZn) <- names(dat_list_Zn)
save(fitsZn, file = "ZnFits_All.RData")

####---- Zinc plots ----####

load("ZnFits_All.RData")

Zn_all_plots <- lapply(fitsZn, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = FALSE, ecx = TRUE, ecx_val = 50) +
    #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})

####---- rename plate replicates ----####

names(Zn_all_plots) <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5", "Rep6", "Rep7", "Rep8", "Rep9", "Rep10",
                           "Rep11", "Rep12", "Rep13", "Rep14", "Rep15", "Rep16", "Rep17")

save(Zn_all_plots, file = "PlotData_Zn_15m_log10.RData")

figure_Zn <- ggpubr::ggarrange(
  plotlist = Zn_all_plots, nrow = 5, ncol = 4, labels = names(Zn_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Zn, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Zinc concentration (mg/L)"))
)

####---- zinc EC50s ----####

ZnEC50 <- lapply(fitsZn, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(ZnEC50, file = "ZnEC50s.csv")
save(ZnEC50, file = "ZnEC50s_15m.RData")
