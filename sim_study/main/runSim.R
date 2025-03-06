##################################################
## Project:tbSTATISpaper1
## Script purpose: This program contains the function to simulate data and estimate TB-STATIS.  
#                  z iterations are run and the sequence and data are output as .rds files. 
#                  Set up to run as a batch job using the iTask variable. 
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################


rm(list = ls())
options(scipen=999)

library(dplyr)
library(data.table)
library(combinat)

setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main")

####Batch Job####

#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#Get sim parameters from arguments 
args <- commandArgs(trailingOnly = TRUE)

#arg1 sample size
M <- as.numeric(args[1])

#arg2 p
p <- as.numeric(args[2])

#arg3 nstart
nstart <- as.numeric(args[3])

#arg4 iteratations 
initial_iter <- as.numeric(args[4])

#arg5 transitions 
N <- as.numeric(args[5])
####

#number of times to run sim 
z=20 

#source functions
source("functions.R")

#N=12
#M=500
#p=.95
#nstart=1
#initial_iter=2
#z=20
#iTask=45


sim <- function(M, N,p, z, nstart, initial_iter){

  #set up objects to store results
  seqs <- matrix(ncol=N)
  true_seqs <- matrix(ncol=N)
  ml_seqs <- matrix(ncol=N)
  all_data <- data.frame(matrix(ncol=N+11))
  seed_vec <- vector( length=z)
  likes <- matrix(ncol=2, nrow=z)
  p_vec= rep(p, N)
  
  for(i in 1:z){
    #set seed that is unique to each data set
    seed=as.numeric(paste(iTask,M/10,p*100, N/2,i, sep=""))
    #print(seed)
    set.seed(seed)
    seed_vec[i] <- seed
    
    #generate number of clinical measures and states per measure
    D=make_D(N)
    info=D
    bio_info_long<- data.frame(bio=NA, event=NA)
    for(t in 1: dim(info)[1]){
    add <- data.frame(bio=rep(info[t,1]), event = c(1:info[t,2]))
    bio_info_long <- rbind(bio_info_long, add)
    }
    bio_info_long <- na.omit(bio_info_long)
    bio_info_long$var_name <- paste0("V", 1:dim(bio_info_long)[1])

    #make data set
    dat <- make_data(bio_info_long, p, M)
    data<- dat[[1]]
    #save true disease sequence
    true = data.frame(dat[[2]])[,6]
    #save likelihood for true S
    l <- get_likelihood(data[3:ncol(data)], dat[[2]], p_vec)[[4]]
    
    #run TB-STATIS
    ml <- fit_STATIS(data[,3:ncol(data)],  p_vec, bio_info_long, nstart, initial_iter )
    
    #save maximum likelihood sequence and add to seqs df
    ml_seq=t(data.frame(ml[[4]]) %>% arrange(pos) %>% select(sub))
    ml_seqs <- rbind(ml_seqs, ml_seq)
    true_seqs=rbind(true_seqs, true)
    
    #save likelihood value for estimated seq
    m <- get_likelihood(data[3:ncol(data)], data.frame(ml[[4]]), p_vec)[[4]]
    
    #get predicted disease class 
    stage=get_stage(data[,3:ncol(data)],S=data.frame(ml[[4]]),p=p_vec)
    
    #update likes matrix with true and estimated likelihood
    likes[i,1] <- l
    likes[i,2] <- m
    
    #additional info to save for each run
    data$M=M
    data$N=N
    data$z=i
    data$run = iTask
    data$p <- p
    data$nbio <- dim(D)[1]
    data$nevents <- max(D$events)
    data$pred_stage=stage$pred_stage
    data$seed=seed
    
    colnames(all_data) <- colnames(data)
    all_data <- rbind(all_data, data)
    
  }
    #organize results
    likes <- data.frame(likes)
    ml_seqs <- data.frame(ml_seqs) %>% na.omit()
    true_seqs <- data.frame(true_seqs) %>% na.omit()
    colnames(likes) <- c("true_like", "ml_like")
    colnames(ml_seqs) <- paste0("ml_", colnames(ml_seqs))
    colnames(true_seqs) <- paste0("true_", colnames(true_seqs))
    res <- cbind(likes, ml_seqs, true_seqs)
    res$M=M
    res$N=N
    res$z=c(1:z)
    res$run=iTask
    res$seed=seed_vec
    res$p <- p
    res$nbio <- dim(D)[1]
    res$nevents <- max(D$events)
    all_data$stage = all_data$stage-1
    
  #return res (all sequences) and all_data (all simulated data sets)
  return(list(res, all_data))
}

#run sim 
sim1 <- sim(M, N,p, z, nstart, initial_iter)

#save results 
saveRDS(sim1[[1]], paste0("./sim_results/sim_seq_",M,"_", p, "_", N,"_", iTask, ".rds"))
saveRDS(sim1[[2]], paste0("./sim_results/sim_data_",M,"_", p, "_", N,"_", iTask, ".rds"))



