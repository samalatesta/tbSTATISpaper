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


#redcap data set file name
data <- read.csv("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/data/REMox_data.csv")

#functions
source("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/functions.R")

measures <- data %>% dplyr::select( CAV, SMEAR1, SMEAR2, SMEAR3, SMEAR4,SWEATS, WTLOSS, COUGH, FEVER)

D=data.frame(bio=c( "cav", "smear", "smear", "smear", "smear","sweats", "wtloss", "cough", "fever"), events=c(1,  1,2, 3,4,1,1,1,1), var_name=c("CAV", "SMEAR1", "SMEAR2", "SMEAR3", "SMEAR4","SWEATS",  "WTLOSS", "COUGH", "FEVER"))
clinical_info=D
p_vec=c(.95,.95,.95, .95, .95, .85, .85, .85, .85)
nstart=20
initial_iter=5000
remox <- fit_STATIS(data=measures, p_vec, D, nstart, initial_iter)
#remox[[4]]

st <- get_stage(measures, data.frame(remox[[4]]), p_vec)
#table(st$pred_stage)

#output results
ml=data.frame(remox[[4]])
write.csv(ml, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/results/remox_ml_seq.csv", row.names=F)

likes=remox[[5]]
write.csv(likes, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/results/remox_ml_likes.csv", row.names = F)

write.csv(st, "/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox/results/remox_class.csv", row.names = F)
#library(tbSTATIS)
#plot_likes(remox[[5]])
