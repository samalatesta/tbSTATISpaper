##################################################
## Project:tbSTATISpaper1
## Script purpose: This program combines all simulated data into two unique data sets per simulated setting. 
#                  One data set for individual-level data and one for all sequence-level data. 
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################
library(dplyr)
library(ggplot2)
library(tidyverse)

#main simulation results

setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/sim_results")

#file names
files_data <-list.files(path="/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/sim_results",pattern = "sim_data*")
files_seq <- list.files(path="/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/sim_results",pattern = "sim_seq*")

#data and stages (N=4)
dat1 <- files_data[str_detect(files_data,"_4_")==T & str_detect(files_data,"_4_")==T]  %>% map(readRDS)
dat1 <- data.frame(Reduce(rbind, dat1))

dat1 <- dat1 %>% na.omit()
dat1$id=paste(dat1$z, "-", dat1$run)
dat1 <- dat1 %>% filter(run <= 50)
write.csv(dat1, paste0("res_data_4_", Sys.Date(), ".csv"))

#sequences (N=4)
seqs1 <-  files_seq[str_detect(files_seq,"_4_")==T] %>% map(readRDS) 
seqs1 <- data.frame(Reduce(rbind, seqs1)) %>% na.omit()
seqs1 <- seqs1  %>% filter(run <= 50)
write.csv(seqs1, paste0("res_seqs_4_", Sys.Date(), ".csv"), row.names=F)

#data and stages (N=8)
dat2 <- files_data[str_detect(files_data,"_8_")==T] %>% map(readRDS)
dat2 <- data.frame(Reduce(rbind, dat2))
dat2 <- dat2 %>% na.omit()
dat2$id=paste(dat2$z, "-", dat2$run)
write.csv(dat2, paste0("res_data_8_", Sys.Date(), ".csv"))

#sequences (N=8)
seqs2 <-  files_seq[str_detect(files_seq,"_8_")==T] %>% map(readRDS) 
seqs2 <- data.frame(Reduce(rbind, seqs2))
write.csv(seqs2, paste0("res_seqs_8_", Sys.Date(), ".csv"), row.names=F)

#data and stages (N=12)
dat3 <- files_data[str_detect(files_data,"_12_")==T] %>% map(readRDS)
dat3 <- data.frame(Reduce(rbind, dat3))
dat3 <- dat3 %>% na.omit()
dat3$id=paste(dat3$z, "-", dat3$run)
write.csv(dat3, paste0("res_data_12_", Sys.Date(), ".csv"))


#sequences (N=12)
seqs3 <-  files_seq[str_detect(files_seq,"_12_")==T] %>% map(readRDS) 
seqs3 <- data.frame(Reduce(rbind, seqs3))
write.csv(seqs3, paste0("res_seqs_12_", Sys.Date(), ".csv"), row.names=F)


#sensitivity analysis results

setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/sensitivity/sim_results")

#file names
files_data <-list.files(path="/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/sensitivity/sim_results",pattern = "sim_data*")
files_seq <- list.files(path="/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/sensitivity/sim_results",pattern = "sim_seq*")

#data and stages (N=4)
dat1 <- files_data[str_detect(files_data,"_4_")==T & str_detect(files_data,"_4_")==T]  %>% map(readRDS)
dat1 <- data.frame(Reduce(rbind, dat1))

dat1 <- dat1 %>% na.omit()
dat1$id=paste(dat1$z, "-", dat1$run)
write.csv(dat1, paste0("res_data_4_sens_", Sys.Date(), ".csv"))

#sequences (N=4)
seqs1 <-  files_seq[str_detect(files_seq,"_4_")==T] %>% map(readRDS) 
seqs1 <- data.frame(Reduce(rbind, seqs1)) %>% na.omit()
write.csv(seqs1, paste0("res_seqs_4_sens_", Sys.Date(), ".csv"), row.names=F)

#data and stages (N=8)
dat2 <- files_data[str_detect(files_data,"_8_")==T] %>% map(readRDS)
dat2 <- data.frame(Reduce(rbind, dat2))
dat2 <- dat2 %>% na.omit()
dat2$id=paste(dat2$z, "-", dat2$run)
write.csv(dat2, paste0("res_data_8_sens_", Sys.Date(), ".csv"))

#sequences (N=8)
seqs2 <-  files_seq[str_detect(files_seq,"_8_")==T] %>% map(readRDS) 
seqs2 <- data.frame(Reduce(rbind, seqs2))
write.csv(seqs2, paste0("res_seqs_8_sens_", Sys.Date(), ".csv"), row.names = F)

#data and stages (N=12)
dat3 <- files_data[str_detect(files_data,"_12_")==T] %>% map(readRDS)
dat3 <- data.frame(Reduce(rbind, dat3))
dat3 <- dat3 %>% na.omit()
dat3$id=paste(dat3$z, "-", dat3$run)
write.csv(dat3, paste0("res_data_12_sens_", Sys.Date(), ".csv"))

#sequences (N=12)
seqs3 <-  files_seq[str_detect(files_seq,"_12_")==T] %>% map(readRDS) 
seqs3 <- data.frame(Reduce(rbind, seqs3))
write.csv(seqs3, paste0("res_seqs_12_sens_", Sys.Date(), ".csv"), row.names=F)
