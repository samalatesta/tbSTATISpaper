library(dplyr)
library(data.table)
library(lubridate)
library(arsenal)
library(flextable)
library(ggplot2)
library(arsenal)
library(dplyr)
library(data.table)
library(combinat)
library(ggpubr)

#### TRUST
#redcap data set file name
input_data <- "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/TRUST/data/TRUST_DATA_2024-02-21_0840.cleaned.wide_noreenroll.csv"

#functions
source("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/functions.R")

#path name to data file
data <- read.csv(input_data, stringsAsFactors = F)

#subset 
data$tot_samps <- rowSums(is.na(data[,c(paste0("culture_conversion_sputum_specimen_", 1:12))])==F)

data$pos <- ifelse(data$culture_conversion_sputum_specimen_1 %in% c("tb_positive", "tb_positive_contaminated")|data$culture_conversion_sputum_specimen_2 %in% c("tb_positive", "tb_positive_contaminated"),1,0)
data <- data %>% dplyr::filter(pos==1 & tot_samps>=3) 

#### data prep for model
vars <- data %>% dplyr::select(pid, bl_cough, bl_sweat, bl_fever, bl_hemop, bl_wtloss, culture_conversion_sputum_specimen_1, s_concafb_sputum_specimen_1, cxr_cavity_chest_radiograph_1, cxr_finding_chest_radiograph_1, bl_hiv, cxr_infiltrate_chest_radiograph_1, cxr_infiltrate_chest_radiograph_1 )

#symptoms
vars$cough <- ifelse(vars$bl_cough > 0, 1, 0)
vars$cough[is.na(vars$cough)] <- 0

vars$sweat <- ifelse(vars$bl_sweat > 0, 1, 0)
vars$sweat[is.na(vars$sweat)] <- 0

vars$fever <- ifelse(vars$bl_fever > 0, 1, 0)
vars$fever[is.na(vars$fever)] <- 0

vars$wtloss <- ifelse(vars$bl_wtloss > 0,1,0)
vars$wtloss[is.na(vars$wtloss)] <- 0

vars$fever <- ifelse(vars$bl_fever > 0,1,0)
vars$fever[is.na(vars$fever)] <- 0

#smear
vars$smear_pos <- as.numeric(factor(vars$s_concafb_sputum_specimen_1, levels = c(0,4,1,2,3), labels = c("0", "1", "1", "1", "1")))-1

vars$smear_3 <- as.numeric(factor(vars$s_concafb_sputum_specimen_1, levels = c(0,4,1,2,3), labels = c("0", "0", "0", "0", "1")))-1


#cavity
vars$cavity <- ifelse(vars$cxr_cavity_chest_radiograph_1 ==1, 1, ifelse(vars$cxr_cavity_chest_radiograph_1 ==0, 0, NA))

#infiltrates
vars$infiltrates <- ifelse(vars$cxr_infiltrate_chest_radiograph_1 %in% c(2), 1, ifelse(vars$cxr_infiltrate_chest_radiograph_1 %in% c(1, 0), 0, NA))

vars$unilatinfiltrates <- ifelse(vars$cxr_infiltrate_chest_radiograph_1 %in% c(1,2), 1, ifelse(vars$cxr_infiltrate_chest_radiograph_1 %in% c(0), 0, NA))

#binary vars
vars_bin <- vars %>% dplyr::select( pid,cough,fever, sweat,wtloss, cavity, smear_pos, smear_3,unilatinfiltrates, infiltrates) %>% na.omit()

#run model
D=data.frame(bio=c("cavity", "wtloss", "sweat", "smear","smear", "infiltrates", "infiltrates", "cough", "fever"), events=c(1,1,1, 1,2,1,2,1, 1), var_name=c("cavity", "wtloss", "sweat", "smear_pos", "smear_3", "unilatinfiltrates", "infiltrates", "cough", "fever"))

measures <- vars_bin %>% dplyr::select(cavity, wtloss, sweat, smear_pos, smear_3, unilatinfiltrates, infiltrates, cough, fever)

clinical_info=D
p_vec=c(.95,.85,.85, .95,.95, .95, .95, .85, .85)
nstart=20
initial_iter=5000
trust <- fit_STATIS(data=measures, p_vec, D, nstart, initial_iter)

#trust[[4]]

#get results and output

st <- get_stage(measures, data.frame(trust[[4]]), p_vec)
#table(st$pred_stage)

#output results
ml=data.frame(trust[[4]])
write.csv(ml, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/TRUST/results/trust_ml_seq.csv", row.names=F)

likes=trust[[5]]
write.csv(likes, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/TRUST/results/trust_ml_likes.csv", row.names = F)

write.csv(st, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/TRUST/results/trust_class.csv", row.names = F)

#library(tbSTATIS)
#plot_likes(trust[[5]])
