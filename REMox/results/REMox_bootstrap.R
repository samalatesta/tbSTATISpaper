rm(list = ls())
options(scipen=999)


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


####Batch Job####

#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
####

#number of times to run sim 
z=1

#redcap data set file name
data <-  read.csv("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/data/REMox_data.csv")

#functions
source("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/functions.R")

measures <- data %>% dplyr::select( CAV, SMEAR1, SMEAR2, SMEAR3, SMEAR4,SWEATS, WTLOSS, COUGH, FEVER)

D=data.frame(bio=c( "cav", "smear", "smear", "smear", "smear","sweats", "wtloss", "cough", "fever"), events=c(1,  1,2, 3,4,1,1,1,1), var_name=c("CAV", "SMEAR1", "SMEAR2", "SMEAR3", "SMEAR4","SWEATS",  "WTLOSS", "COUGH", "FEVER"))
clinical_info=D
p_vec=c(.95,.95,.95, .95, .95, .85, .85, .85, .85)
nstart=20
initial_iter=5000

#bootstrap
boot_ml = data.frame()
for(i in 1:z){
boot = measures[sample(nrow(measures),nrow(measures), T),]
remox <- fit_STATIS(data=measures, p_vec, D, nstart, initial_iter)
ml <- t(data.frame(remox[[4]]) %>% dplyr::arrange(pos) %>% select(sub))
boot_ml = rbind(boot_ml, ml)
}

write.csv(boot_ml, paste0("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/results/REMox_boot_ml_",iTask,".csv"), row.names = F)

