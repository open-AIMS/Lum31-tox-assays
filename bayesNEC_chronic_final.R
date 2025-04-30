##### Chronic Exp using fresh cells; measured metal data #######

library(tidyverse)
library(readxl)
library(bayesnec)
library(ggpubr)

############### plate 1 Cu ################

dat1 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Cu_Plate1") %>%
  data.frame() 

dat1 <- dat1 |> 
  mutate(x = as.numeric(as.character(x)),     
         log_x = log(x))

dat1_Manip <- dat1 |>  
  mutate(y24h_qtox = (1-(median(dat1[1:4,1])- y24h)/median(dat1[1:4,1])),
         cells_qtox = (1-(median(dat1[1:4,2])-cells_ml)/median(dat1[1:4,2])),
         sgr = ((log(dat1[ ,2]) - log(5.72e5))/24), #average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

dat1_Manip <- dat1_Manip[-15, ] #remove outlier identified as pipetting error

ggplot(dat1_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log copper conc", y="SGR", title = "")

ggplot(dat1_Manip, aes(x=log_x, y=y24h_qtox))+
  geom_point(size=2)+
  labs(x="log copper conc", y="proprotion of luminescence", title = "")

ggplot(dat1_Manip, aes(x=cells_ml, y=y24h))+
  geom_point(size=2)+
  labs(x="cells/mL", y="luminescence (RLU)", title = "")

############### plate 2 Cu ################

dat2 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Cu_Plate2") %>%
  data.frame() 

dat2 <- dat2 |> 
  mutate(x = as.numeric(as.character(x)), 
         log_x = log(x))

dat2_Manip <- dat2 |>  
  mutate(y24h_qtox = (1-(median(dat2[1:4,1])- y24h)/median(dat2[1:4,1])),
         cells_qtox = (1-(median(dat2[1:4,2])-cells_ml)/median(dat2[1:4,2])),
         sgr = ((log(dat2[ ,2]) - log(4.20E+05))/24), #average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

ggplot(dat2_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log copper conc", y="SGR", title = "")

ggplot(dat2_Manip, aes(x=log_x, y=y24h_qtoxB))+
  geom_point(size=2)+
  labs(x="log copper conc", y="proprotion of luminescence", title = "")

############### plate 3 Cu ################

dat3 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Cu_Plate3") %>%
  data.frame() 

dat3 <- dat3 |> 
  mutate(x = as.numeric(as.character(x)), 
         log_x = log(x))

dat3_Manip <- dat3 |>  
  mutate(y24h_qtox = (1-(median(dat3[1:4,1])- y24h)/median(dat3[1:4,1])),
         cells_qtox = (1-(median(dat3[1:4,2])-cells_ml)/median(dat3[1:4,2])),
         sgr = ((log(dat3[ ,2]) - log(3.52E+05))/24), #average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

ggplot(dat3_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log copper conc", y="SGR", title = "")

ggplot(dat3_Manip, aes(x=log_x, y=y24h_qtoxB))+
  geom_point(size=2)+
  labs(x="log copper conc", y="proprotion of luminescence", title = "")

####---- combined Cu data for HPC ----####

# log_x and all qtoxB data for chronic endpoints

dat1_Manip2 <- dat1_Manip |> select(3, 7:9)
dat1_Manip2$sgrB[dat1_Manip2$sgrB < 0] <- 0 #change negative sgr's to 0

dat2_Manip2 <- dat2_Manip |> select(3, 7:9)
dat2_Manip2$sgrB[dat2_Manip2$sgrB < 0] <- 0

dat3_Manip2 <- dat3_Manip |> select(3, 7:9)
dat3_Manip2$sgrB[dat3_Manip2$sgrB < 0] <- 0

names(dat1_Manip2) <- c("log_x", "y24h_qtoxB_p1",  "cells_qtoxB_p1" ,"sgrB_p1")
names(dat2_Manip2) <- c("log_x", "y24h_qtoxB_p2",  "cells_qtoxB_p2" ,"sgrB_p2")
names(dat3_Manip2) <- c("log_x", "y24h_qtoxB_p3",  "cells_qtoxB_p3" ,"sgrB_p3")

Dat1 <- lapply(2:ncol(dat1_Manip2), FUN = function(i) {
  dat.i <- dat1_Manip2[,c("log_x", colnames(dat1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(Dat1) <- colnames(dat1_Manip2)[-1]

Dat2 <- lapply(2:ncol(dat2_Manip2), FUN = function(i) {
  dat.i <- dat2_Manip2[,c("log_x", colnames(dat2_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(Dat2) <- colnames(dat2_Manip2)[-1]

Dat3 <- lapply(2:ncol(dat3_Manip2), FUN = function(i) {
  dat.i <- dat3_Manip2[,c("log_x", colnames(dat3_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(Dat3) <- colnames(dat3_Manip2)[-1]

datCu_list <- c(Dat1, Dat2, Dat3)

#save(datCu_list, file = "Cu_All_Chronic_Data.RData")

################# HPC Cu ##################
library(bayesnec)
library(future.apply)

load(file = "Cu_All_Chronic_Data.RData")

plan(multisession)
fitsCu_C <- future_lapply(datCu_list, FUN = function(i) {
  
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsCu_C) <- names(datCu_list)
save(fitsCu_C, file = "CuFits_All_chronic.RData")

########################################

####---- Copper plots ----####

load("CuFits_All_chronic.RData")

Cu_all_plots <- lapply(fitsCu_C, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = TRUE, ecx = TRUE, ecx_val = 50) +
    scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Cu <- ggpubr::ggarrange(
  plotlist = Cu_all_plots, nrow = 3, ncol = 3, labels = names(Cu_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Cu, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Copper Concentration (mg/L)"))
)

######################## just copper SGR plots ################

Cu_plates <- list(Copper_plate1 = fitsCu_C$sgrB_p1, Copper_plate2 = fitsCu_C$sgrB_p2, Copper_plate3 = fitsCu_C$sgrB_p3)

Cu_all_plots2 <- lapply(Cu_plates, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = TRUE, ecx = TRUE, ecx_val = 50) +
    #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Cu2 <- ggpubr::ggarrange(
  plotlist = Cu_all_plots2, nrow = 3, ncol = 1,
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Cu2, left = text_grob(expression("Specific Growth Rate (scaled)"), rot = 90),
  bottom = text_grob(expression("Copper Concentration (mg/L)"))
)

### cycle through model outputs
summary(fitsCu_C$y24h_qtoxB)
check_chains(fitsCu_C$y24h_qtoxB)
rhat(fitsCu_C$y24h_qtoxB, rhat_cutoff = 1.05)

####---- Copper EC50s ----####

CuEC50s <- lapply(fitsCu_C, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(CuEC50s, file = "CuEC50s_AllChronicEndPts.csv")
save(CuEC50s, file = "CuEC50s_AllChronicEndPts.RData")

############### plate 1 Zn ################

datZn1 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Zn_Plate1") %>%
  data.frame() 

datZn1 <- datZn1 |> 
  mutate(x = as.numeric(as.character(x)), 
         log_x = log(x))

datZn1_Manip <- datZn1 |>  
  mutate(y24h_qtox = (1-(median(datZn1[1:4,1])- y24h)/median(datZn1[1:4,1])),
         cells_qtox = (1-(median(datZn1[1:4,2])-cells_ml)/median(datZn1[1:4,2])),
         sgr = ((log(datZn1[ ,2]) - log(5.72e5))/24), #used average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

ggplot(datZn1_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log zinc sulfate conc", y="SGR", title = "")

####---- plate 2 Zn ----####

datZn2 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Zn_Plate2") %>%
  data.frame() 

datZn2 <- datZn2 |> 
  mutate(x = as.numeric(as.character(x)), 
         log_x = log(x))

datZn2_Manip <- datZn2 |>  
  mutate(y24h_qtox = (1-(median(datZn2[1:4,1])- y24h)/median(datZn2[1:4,1])),
         cells_qtox = (1-(median(datZn2[1:4,2])-cells_ml)/median(datZn2[1:4,2])),
         sgr = ((log(datZn2[ ,2]) - log(4.20E+05))/24), #average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

ggplot(datZn2_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log zinc sulfate conc", y="SGR", title = "")

ggplot(datZn2_Manip, aes(x=log_x, y=y24h_qtox))+
  geom_point(size=2)+
  labs(x="log zinc sulfate conc", y="proprotion of luminescence", title = "")

####---- plate 3 Zn ----####

datZn3 <- read_excel("Cu_Zn_Chronic_measured_AIMSrepository.xlsx", sheet="Zn_Plate3") %>%
  data.frame() 

datZn3 <- datZn3 |> 
  mutate(x = as.numeric(as.character(x)), 
         log_x = log(x))

datZn3_Manip <- datZn3 |>  
  mutate(y24h_qtox = (1-(median(datZn3[1:4,1])- y24h)/median(datZn3[1:4,1])),
         cells_qtox = (1-(median(datZn3[1:4,2])-cells_ml)/median(datZn3[1:4,2])),
         sgr = ((log(datZn3[ ,2]) - log(3.52E+05))/24), #average T0 cell count for 11 cells
         y24h_qtoxB = y24h_qtox/max(y24h_qtox),
         cells_qtoxB = cells_qtox/max(cells_qtox),
         sgrB = sgr/max(sgr)) #bound between 0-1 to model as beta like all others

ggplot(datZn3_Manip, aes(x=log_x, y=sgr))+
  geom_point(size=2)+
  labs(x="log zinc sulfate conc", y="SGR", title = "")

ggplot(datZn3_Manip, aes(x=log_x, y=y24h_qtox))+
  geom_point(size=2)+
  labs(x="log zinc sulfate conc", y="proprotion of luminescence", title = "")

####---- combined Zn data for HPC ----####

# select only log_x and qtoxB data 

datZn1_Manip2 <- datZn1_Manip |> select(3, 7:9)
datZn1_Manip2$sgrB[datZn1_Manip2$sgrB < 0] <- 0 #change negative sgr's to 0

datZn2_Manip2 <- datZn2_Manip |> select(3, 7:9)
datZn2_Manip2$sgrB[datZnC_Manip2$sgrB < 0] <- 0

datZn3_Manip2 <- datZn3_Manip |> select(3, 7:9)
datZn3_Manip2$sgrB[datZn3_Manip2$sgrB < 0] <- 0

names(datZn1_Manip2) <- c("log_x", "y24h_qtoxB_p1",  "cells_qtoxB_p1" ,"sgrB_p1")
names(datZn2_Manip2) <- c("log_x", "y24h_qtoxB_p2",  "cells_qtoxB_p2" ,"sgrB_p2")
names(datZn3_Manip2) <- c("log_x", "y24h_qtoxB_p3",  "cells_qtoxB_p3" ,"sgrB_p3")


DatZn1 <- lapply(2:ncol(datZn1_Manip2), FUN = function(i) {
  dat.i <- datZn1_Manip2[,c("log_x", colnames(datZn1_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(DatZn1) <- colnames(datZn1_Manip2)[-1]

DatZn2 <- lapply(2:ncol(datZn2_Manip2), FUN = function(i) {
  dat.i <- datZn2_Manip2[,c("log_x", colnames(datZn2_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(DatZn2) <- colnames(datZn2_Manip2)[-1]

DatZn3 <- lapply(2:ncol(datZn3_Manip2), FUN = function(i) {
  dat.i <- datZn3_Manip2[,c("log_x", colnames(datZn3_Manip2)[i])]
  colnames(dat.i) <- c("log_x", "y") 
  dat.i
})
names(DatZn3) <- colnames(datZn3_Manip2)[-1]

datZn_list <- c(DatZn1, DatZn2, DatZn3)

#save(datZn_list, file = "Zn_All_Chronic_Data.RData")

################# HPC Zn ##################
library(bayesnec)
library(future.apply)

load(file = "Zn_All_Chronic_Data.RData")

plan(multisession)
fitsZn_C <- future_lapply(datZn_list, FUN = function(i) {
  
  bnec(y~crf(log_x, model = "all"),
       data = i, seed = 333, control = list(adapt_delta = 0.99),
       iter = 20000, warmup = 19000)
})

names(fitsZn_C) <- names(datZn_list)
save(fitsZn_C, file = "ZnFits_All_chronic.RData")
######################################

####---- Zinc plots ----####

load("ZnFits_All_Chronic.RData")

Zn_all_plots <- lapply(fitsZn_C, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = TRUE, ecx = TRUE, ecx_val = 50) +
    scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Zn <- ggpubr::ggarrange(
  plotlist = Zn_all_plots, nrow = 3, ncol = 3, labels = names(Zn_all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Zn, left = text_grob(expression("Proportion of Luminescence (scaled)"), rot = 90),
  bottom = text_grob(expression("Zinc Concentration (mg/L)"))
)

############# Just Zinc SGR plots ###############

Zn_plates <- list(Zinc_plate1 = fitsZn_C$sgrB_p1, Zinc_plate2 = fitsZn_C$sgrB_p2, Zinc_plate3 = fitsZn_C$sgrB_p3)

Zn_all_plots2 <- lapply(Zn_plates, function(x) {
  autoplot(x, xform=function(x){exp(x)}, nec = TRUE, ecx = TRUE, ecx_val = 50) +
    #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) + 
    scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure_Zn2 <- ggpubr::ggarrange(
  plotlist = Zn_all_plots2, nrow = 3, ncol = 1,
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure_Zn2, left = text_grob(expression("Specific Growth Rate (scaled)"), rot = 90),
  bottom = text_grob(expression("Zinc Concentration (mg/L)"))
)

### cycle through model outputs
summary(fitsZn_C$sgrB_p1)
check_chains(fitsZn_C$sgrB_p1)
rhat(fitsZn_C$sgrB_p1, rhat_cutoff = 1.05)

####---- Zinc EC50s ----####

ZnEC50s <- lapply(fitsZn_C, function(x) {
  ecx(x, ecx_val = 50, xform = function(x){exp(x)})
})
write.csv(ZnEC50s, file = "ZnEC50s_AllChronicEndPts.csv")
save(ZnEC50s, file = "ZnEC50s_AllChronicEndPts.RData")


####---- combined Cu/Zn SGR plot ----####

comb <- ggpubr::ggarrange(figure_Cu2, figure_Zn2, ncol = 2)

ggpubr::annotate_figure(
  comb, left = text_grob(expression("Specific Growth Rate (normalised)"), rot = 90),
  bottom = text_grob(expression("Concentration (mg/L)"))
)

#ggsave("combPlot_CuZnSGRv3.pdf", plot = comb, device = "pdf", width = 6.5, height = 7)

#library(cowplot)

#combF <- ggdraw(comb) +
#  # draw_label("Specific Growth Rate (normalised)", x = 0.01, y = 0.50, angle = 90, size = 12) +
#  draw_label("Copper (mg/L)", x = 0.27, y = 0.01, size = 12) +
#  draw_label("Zinc (mg/L)", x = 0.77, y = 0.01, size = 12) 
#combF


####---- posterior comparisons of EC50s ----####

CuPostComp <- compare_posterior(Cu_plates,
                                comparison = "ecx", ecx_val = 50)

ggplot(data = CuPostComp$posterior_data, mapping = aes(x = value)) + 
  geom_density(mapping = aes(group = model, colour = model, fill = model),
               alpha = 0.3) +
  theme_classic()

ZnPostComp <- compare_posterior(Zn_plates,
                                comparison = "ecx", ecx_val = 50)

ggplot(data = ZnPostComp$posterior_data, mapping = aes(x = value)) + 
  geom_density(mapping = aes(group = model, colour = model, fill = model),
               alpha = 0.3) +
  theme_classic()
