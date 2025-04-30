
library(tidyverse)
library(readxl)
library(bayesnec)
library(ggpubr)
#library(future.apply) #on HPC 

###---- load data: copper ----####

# five lots of data corresponding to different chemistry sampling

CuDat1 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="19Mar24") %>%
  data.frame() 

CuDat1_Manip <- CuDat1 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y30_qtox = (1-(median(CuDat1[1:4,4])-Rep1B_y30)/median(CuDat1[1:4,4])), 
         Rep2B_y30_qtox = (1-(median(CuDat1[1:4,7])-Rep2B_y30)/median(CuDat1[1:4,7])),
         Rep3A_y30_qtox = (1-(median(CuDat1[1:4,10])-Rep3A_y30)/median(CuDat1[1:4,10])),
         Rep4A_y30_qtox = (1-(median(CuDat1[1:4,13])-Rep4A_y30)/median(CuDat1[1:4,13])),
         Rep1B_y30_qtoxB = Rep1B_y30_qtox/max(Rep1B_y30_qtox),
         Rep2B_y30_qtoxB = Rep2B_y30_qtox/max(Rep2B_y30_qtox),
         Rep3A_y30_qtoxB = Rep3A_y30_qtox/max(Rep3A_y30_qtox),
         Rep4A_y30_qtoxB = Rep4A_y30_qtox/max(Rep4A_y30_qtox))

CuDat2 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="28Mar24") %>%
  data.frame() 

CuDat2_Manip <- CuDat2 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y30_qtox = (1-(median(CuDat2[1:4,4])-Rep1B_y30)/median(CuDat2[1:4,4])), 
         Rep2B_y30_qtox = (1-(median(CuDat2[1:4,7])-Rep2B_y30)/median(CuDat2[1:4,7])),
         Rep1B_y30_qtoxB = Rep1B_y30_qtox/max(Rep1B_y30_qtox),
         Rep2B_y30_qtoxB = Rep2B_y30_qtox/max(Rep2B_y30_qtox))

CuDat3 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="5Apr24") %>%
  data.frame() 

CuDat3_Manip <- CuDat3 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep3B_y30_qtox = (1-(median(CuDat3[1:4,7])-Rep3B_y30)/median(CuDat3[1:4,4])), 
         Rep4B_y30_qtox = (1-(median(CuDat3[1:4,7])-Rep4B_y30)/median(CuDat3[1:4,7])),
         Rep5B_y30_qtox = (1-(median(CuDat3[1:4,10])-Rep5B_y30)/median(CuDat3[1:4,10])),
         Rep3B_y30_qtoxB = Rep3B_y30_qtox/max(Rep3B_y30_qtox),
         Rep4B_y30_qtoxB = Rep4B_y30_qtox/max(Rep4B_y30_qtox),
         Rep5B_y30_qtoxB = Rep5B_y30_qtox/max(Rep5B_y30_qtox))

CuDat4 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="23Apr24") %>%
  data.frame() 

CuDat4_Manip <- CuDat4 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep6B_y30_qtox = (1-(median(CuDat4[1:4,4])-Rep6B_y30)/median(CuDat4[1:4,4])), 
         Rep7B_y30_qtox = (1-(median(CuDat4[1:4,7])-Rep7B_y30)/median(CuDat4[1:4,7])),
         Rep8B_y30_qtox = (1-(median(CuDat4[1:4,10])-Rep8B_y30)/median(CuDat4[1:4,10])),
         Rep9B_y30_qtox = (1-(median(CuDat4[1:4,13])-Rep9B_y30)/median(CuDat4[1:4,13])),
         Rep10B_y30_qtox = (1-(median(CuDat4[1:4,16])-Rep10B_y30)/median(CuDat4[1:4,16])),
         Rep6B_y30_qtoxB = Rep6B_y30_qtox/max(Rep6B_y30_qtox),
         Rep7B_y30_qtoxB = Rep7B_y30_qtox/max(Rep7B_y30_qtox),
         Rep8B_y30_qtoxB = Rep8B_y30_qtox/max(Rep8B_y30_qtox),
         Rep9B_y30_qtoxB = Rep9B_y30_qtox/max(Rep9B_y30_qtox),
         Rep10B_y30_qtoxB = Rep10B_y30_qtox/max(Rep10B_y30_qtox))

CuDat5 <- read_excel("Cu Acute tests_measured_AIMSrepository.xlsx", sheet="1Oct24") %>%
  data.frame() 

CuDat5_Manip <- CuDat5 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Batch6_A1_y30_qtox = (1-(median(CuDat5[1:4,4])-Batch6_A1_y30)/median(CuDat5[1:4,4])), 
         Batch6_A2_y30_qtox = (1-(median(CuDat5[1:4,7])-Batch6_A2_y30)/median(CuDat5[1:4,7])),
         Batch6_A1_y30_qtoxB = Batch6_A1_y30_qtox/max(Batch6_A1_y30_qtox),
         Batch6_A2_y30_qtoxB = Batch6_A2_y30_qtox/max(Batch6_A2_y30_qtox))
        

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

datCu30_list <- c(dat1, dat2, dat3, dat4, dat5)

save(datCu30_list, file = "CuData30m_all.RData")

########## HPC START ########
fitsCu30 <- future_lapply(datCu30_list, FUN = function(i) {
  
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsCu30) <- names(datCu30_list)
save(fitsCu30, file = "CuFits_All_30m.RData")
########## HPC END ########

load("CuFits_All_30m.RData")

#---- copper plots ----#

Cu30_all_plots <- lapply(fitsCu30, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = FALSE, ecx = TRUE, ecx_val = 50) +
    #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Cu30 <- ggpubr::ggarrange(
  plotlist = Cu30_all_plots, nrow = 5, ncol = 4, labels = names(Cu30_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Cu30, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Copper concentration (mg/L)"))
)

####---- rename plate replicates ----####

names(Cu30_all_plots) <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5", "Rep6", "Rep7", "Rep8", "Rep9", "Rep10",
                         "Rep11", "Rep12", "Rep13", "Rep14", "Rep15", "Rep16")

save(Cu30_all_plots, file = "PlotData_Cu_30m_log10.RData")

####---- copper EC50s ----####

CuEC50_30m <- lapply(fitsCu30, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(CuEC50_30m, file = "CuEC50s_30m.csv")
save(CuEC50_30m, file = "CuEC50s_30m.RData")

###---- load data: zinc ----####

ZnDat1 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="19Mar24") %>%
  data.frame() 

ZnDat1_Manip <- ZnDat1 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep1B_y30_qtox = (1-(median(ZnDat1[1:4,4])-Rep1B_y30)/median(ZnDat1[1:4,4])), 
         Rep2B_y30_qtox = (1-(median(ZnDat1[1:4,7])-Rep2B_y30)/median(ZnDat1[1:4,7])),
         Rep3A_y30_qtox = (1-(median(ZnDat1[1:4,10])-Rep3A_y30)/median(ZnDat1[1:4,10])),
         Rep4A_y30_qtox = (1-(median(ZnDat1[1:4,13])-Rep4A_y30)/median(ZnDat1[1:4,13])),
         Rep1B_y30_qtoxB = Rep1B_y30_qtox/max(Rep1B_y30_qtox),
         Rep2B_y30_qtoxB = Rep2B_y30_qtox/max(Rep2B_y30_qtox),
         Rep3A_y30_qtoxB = Rep3A_y30_qtox/max(Rep3A_y30_qtox),
         Rep4A_y30_qtoxB = Rep4A_y30_qtox/max(Rep4A_y30_qtox))

ZnDat2 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="28Mar24") %>%
  data.frame() 

ZnDat2_Manip <- ZnDat2 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep5B_y30_qtox = (1-(median(ZnDat2[1:4,3])-Rep5B_y30)/median(ZnDat2[1:4,3])),
         Rep1B_y30_qtox = (1-(median(ZnDat2[1:4,6])-Rep1B_y30)/median(ZnDat2[1:4,6])), 
         Rep2B_y30_qtox = (1-(median(ZnDat2[1:4,9])-Rep2B_y30)/median(ZnDat2[1:4,9])),
         Rep5B_y30_qtoxB = Rep5B_y30_qtox/max(Rep5B_y30_qtox),
         Rep1B_y30_qtoxB = Rep1B_y30_qtox/max(Rep1B_y30_qtox),
         Rep2B_y30_qtoxB = Rep2B_y30_qtox/max(Rep2B_y30_qtox))

ZnDat3 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="5Apr24") %>%
  data.frame() 

ZnDat3_Manip <- ZnDat3 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep3B_y30_qtox = (1-(median(ZnDat3[1:4,7])-Rep3B_y30)/median(ZnDat3[1:4,4])), 
         Rep4B_y30_qtox = (1-(median(ZnDat3[1:4,7])-Rep4B_y30)/median(ZnDat3[1:4,7])),
         Rep5B_y30_qtox = (1-(median(ZnDat3[1:4,10])-Rep5B_y30)/median(ZnDat3[1:4,10])),
         Rep3B_y30_qtoxB = Rep3B_y30_qtox/max(Rep3B_y30_qtox),
         Rep4B_y30_qtoxB = Rep4B_y30_qtox/max(Rep4B_y30_qtox),
         Rep5B_y30_qtoxB = Rep5B_y30_qtox/max(Rep5B_y30_qtox))

ZnDat4 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="23Apr24") %>%
  data.frame() 

ZnDat4_Manip <- ZnDat4 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Rep6B_y30_qtox = (1-(median(ZnDat4[1:4,4])-Rep6B_y30)/median(ZnDat4[1:4,4])), 
         Rep7B_y30_qtox = (1-(median(ZnDat4[1:4,7])-Rep7B_y30)/median(ZnDat4[1:4,7])),
         Rep8B_y30_qtox = (1-(median(ZnDat4[1:4,10])-Rep8B_y30)/median(ZnDat4[1:4,10])),
         Rep9B_y30_qtox = (1-(median(ZnDat4[1:4,13])-Rep9B_y30)/median(ZnDat4[1:4,13])),
         Rep10B_y30_qtox = (1-(median(ZnDat4[1:4,16])-Rep10B_y30)/median(ZnDat4[1:4,16])),
         Rep6B_y30_qtoxB = Rep6B_y30_qtox/max(Rep6B_y30_qtox),
         Rep7B_y30_qtoxB = Rep7B_y30_qtox/max(Rep7B_y30_qtox),
         Rep8B_y30_qtoxB = Rep8B_y30_qtox/max(Rep8B_y30_qtox),
         Rep9B_y30_qtoxB = Rep9B_y30_qtox/max(Rep9B_y30_qtox),
         Rep10B_y30_qtoxB = Rep10B_y30_qtox/max(Rep10B_y30_qtox))

ZnDat5 <- read_excel("Zn Acute tests_measured_AIMSrepository.xlsx", sheet="1Oct24") %>%
  data.frame() 

ZnDat5_Manip <- ZnDat5 %>% 
  mutate(x = as.numeric(as.character(x)),
         log_x = log(x),
         Batch6_A1_y30_qtox = (1-(median(ZnDat5[1:4,4])-Batch6_A1_y30)/median(ZnDat5[1:4,4])), 
         Batch6_A2_y30_qtox = (1-(median(ZnDat5[1:4,7])-Batch6_A2_y30)/median(ZnDat5[1:4,7])),
         Batch6_A1_y30_qtoxB = Batch6_A1_y30_qtox/max(Batch6_A1_y30_qtox),
         Batch6_A2_y30_qtoxB = Batch6_A2_y30_qtox/max(Batch6_A2_y30_qtox))
       
# select only log_x and qtoxB data to use as input for bayesnec

ZnDat1_Manip2 <- ZnDat1_Manip |> select(14,19:22)
ZnDat2_Manip2 <- ZnDat2_Manip |> select(11,15:17)
ZnDat3_Manip2 <- ZnDat3_Manip |> select(11,15:17)
ZnDat4_Manip2 <- ZnDat4_Manip |> select(17,23:27)
ZnDat5_Manip2 <- ZnDat5_Manip |> select(8,11:12)

# create datasets to run models at once

dat1 <- lapply(2:ncol(ZnDat1_Manip2), FUN = function(i) {
  dat.i <- ZnDat1_Manip2[,c("log_x", colnames(ZnDat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat1) <- colnames(ZnDat1_Manip2)[-1]

dat2 <- lapply(2:ncol(ZnDat2_Manip2), FUN = function(i) {
  dat.i <- ZnDat2_Manip2[,c("log_x", colnames(ZnDat2_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat2) <- colnames(ZnDat2_Manip2)[-1]

dat3 <- lapply(2:ncol(ZnDat3_Manip2), FUN = function(i) {
  dat.i <- ZnDat3_Manip2[,c("log_x", colnames(ZnDat3_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat3) <- colnames(ZnDat3_Manip2)[-1]

dat4 <- lapply(2:ncol(ZnDat4_Manip2), FUN = function(i) {
  dat.i <- ZnDat4_Manip2[,c("log_x", colnames(ZnDat4_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat4) <- colnames(ZnDat4_Manip2)[-1]

dat5 <- lapply(2:ncol(ZnDat5_Manip2), FUN = function(i) {
  dat.i <- ZnDat5_Manip2[,c("log_x", colnames(ZnDat5_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(dat5) <- colnames(ZnDat5_Manip2)[-1]

datZn30_list <- c(dat1, dat2, dat3, dat4, dat5)

save(datZn30_list, file = "ZnData30m_all.RData")


####---- Zinc plots ----####

load("ZnFits_All_30m.RData")

Zn30_all_plots <- lapply(fitsZn30, function(x) {
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

names(Zn30_all_plots) <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5", "Rep6", "Rep7", "Rep8", "Rep9", "Rep10",
                         "Rep11", "Rep12", "Rep13", "Rep14", "Rep15", "Rep16", "Rep17")

save(Zn30_all_plots, file = "PlotData_Zn_30m_log10.RData")

figure_Zn30 <- ggpubr::ggarrange(
  plotlist = Zn30_all_plots, nrow = 5, ncol = 4, labels = names(Zn30_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Zn30, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Zinc concentration (mg/L)"))
)

####---- zinc EC50s ----####

ZnEC50_30m <- lapply(fitsZn30, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(ZnEC50_30m, file = "ZnEC50s_30m.csv")
save(ZnEC50_30m, file = "ZnEC50s_30m.RData")
