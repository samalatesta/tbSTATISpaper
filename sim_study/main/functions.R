##################################################
## Project:tbSTATISpaper1
## Script purpose: This program contains all functions to run TB-STATIS and generate data for a simulation study.
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################

library(dplyr)

#' Propose simultaneous events
#'
#' @param info A data frame.
#' @param all_seqs List of all possible sequences.
#' @return A ggplot object.
#' @keywords internal
#'

get_stage <- function(data=data.frame(), S=data.frame(), p=vector()){
  likelihood=get_likelihood(data, S, p )
  
  stage_probs <- likelihood[[2]]
  
  
  stage_probs=data.frame(stage_probs)
  colnames(stage_probs) = c(0:max(S$sub))
  
  pred_stage=as.numeric(colnames(stage_probs)[apply(stage_probs,1,which.max)])
  
  stagedf <- cbind(data,pred_stage)
  
  return(stagedf)
}


make_D <- function(N){
  
  nbio <- sample(2:N, 1)
  D <- data.frame(bio=paste0("bio", 1:nbio), events=rep(1, nbio))
  
  repeat{
    if(sum(D$events)==N){
      break
    }
    rowadd <- sample(1:nbio, 1)
    
    D$events[rowadd] <- D$events[rowadd]+1
    
  }
  
  return(D)
}


get_group <- function( info=data.frame(), all_seqs){

  info <- info %>% dplyr::arrange(order)
  N=dim(info)[1]
  
  all_seqs2=all_seqs
  all_seqs=data.frame(t(apply(all_seqs[,1:N], 1, function(i) paste(i, info$bio) )))
  all_seqs$unique= apply(all_seqs[,1:N], 1, function(x) length(unique(x)))==N
  
  final <-  cbind(all_seqs2, data.frame(unique=all_seqs$unique))#%>% dplyr::filter(unique==T)
  select <- sample(1:dim(final)[1],1)
  
  info$sub<- as.numeric(final[select,1:N])

  
  # print(info)
  return(info)
}


#' Propose sequence as start point for maximum likelihood estimation.
#'
#' @param info A data frame.
#' @return A data frame.
#' @keywords internal
#'
#### Initialize sequence events
get_seq <- function(info=data.frame()){
  #info=make_D(5)
  #transform D to long form
  #bio_info_long<- data.frame(bio=NA, event=NA)
  #for(i in 1: dim(info)[1]){
   #add <- data.frame(bio=rep(info[i,1]), event = c(1:info[i,2]))
   #bio_info_long <- rbind(bio_info_long, add)
  
  #}
  
  #bio_info_long <- bio_info_long %>% dplyr::filter(is.na(event)==F)
  bio_info_long <- info
  colnames(bio_info_long) <- c("bio", "event", "var_name")
  #bio_info_long <- ml[[2]][[1]]
  N=dim(bio_info_long)[1]
  bio_info_long$pos = c(1:N)
  
  #sample(bio_info_long$bio, N)
  #all_perms <- combinat::permn(bio_info_long$bio)
  #r <- sample(1:length(all_perms),1)
  #temp <- all_perms[[r]]
  
  temp<-sample(bio_info_long$bio, N)
  
  bio_info_long_temp <- bio_info_long
  bio_info_long_temp$temp <- temp
  bio_info_long_temp <- bio_info_long_temp %>% dplyr::arrange(temp) 
  bio_info_long_temp$order <- c(1:dim(bio_info_long_temp)[1])
  #bio_info_long_temp
  bio_info_long_new <- bio_info_long #%>% dplyr::arrange(pos) %>% dplyr::select(-temp)
  bio_info_long_new$order <- bio_info_long_temp$pos
  
  return(bio_info_long_new)
}




#' Generate data frame of all possible simultaneous sequences.
#'
#' @param N Number of events.
#' @return A data frame.
#' @keywords internal
#'
#


possible_seqs <- function(N=numeric()){
  
  filter_rows <- function(g, x) {
    ok  <- function(z) all(diff(z) %in% 0:1)
    out <- g[apply(g, 1, ok), ]
    replace(out, TRUE, lapply(out, \(i) x[i]))
  }
  f2 <- function(x = c(1:N), n=N+1, n1=2) {
    data.frame(as.list(rep(1, n1)),
               gtools::combinations(length(x), n-n1, repeats.allowed = TRUE)) |>
      filter_rows(x)
  }
  
  # test runs
  all=f2() # as per question
  all <- all[,-1]
  return(all)
}

#D<- make_D(4)
#p=.95
#M=100
#S=seq
get_likelihood <- function(data=data.frame(), S=data.frame(), p_vec=vector()){

  bio <- data.frame(data[,1:ncol(data)])
  
  # number of biomarkers
  N = as.numeric(dim(bio)[2])
  
  #individuals in data set
  M = dim(bio)[1]
  
  new_bio <- bio
  
  #multiply each p by each corresponding variable
  for(i in 1:dim(new_bio)[2]){
    #x<- S$bio[i]
    #n <- sum(S$bio==x)
    new_bio[,i] <- ifelse(new_bio[,i]==1, dbinom(1,1,p_vec[i]), dbinom(0,1,p_vec[i]))
  }
  
  
  #empty matrix to store individual stage probabilities
  p_perm_k = matrix(NA,M,N+1)
  
  #order based on sequence S and length
  #order = as.numeric(S$pos)
  
  group = as.numeric(S$sub)
  
  #reorder biomarker columns based on S
  bio_order = new_bio
  
  #add alpha to make more continuous
  #alpha=.1
  
  #bio_order=bio_order+alpha
  normal <- 1-bio_order
  
  #special case stage 1
  normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
  #normal_prob[which(normal_prob==0)] <- .1
  
  tot_prob_stage = normal_prob[,1]
  
  p_perm_k[,1] <- tot_prob_stage
  S2=S
  
  for (i in 1:length(unique(group))) {
    #print(i)
    
    if(i < max(group)){
      #abnormal biomarkers assuming S is true sequence
      abnormal_cols <- S2$var_name[S2$sub<=i]
      
      abnormal <- data.frame(bio_order[,abnormal_cols])
      
      #normal data assuming S is true sequence
      normal_cols = S2 %>% filter(!(var_name %in% abnormal_cols)) %>% select(var_name) 
      
      normal <- 1-data.frame(bio_order[,normal_cols[,1]])
      normal[normal==-1]<- 1
      #abnormal_n <- apply(abnormal, 1, sum)
      #normal_n <- apply(normal, 1, function(x) sum(x==0))
      
      
      
      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #abnormal_prob[which(abnormal_prob==0)] <- .1
      
      normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #normal_prob[which(normal_prob==0)] <- .1
      #fill in empty matrix
      #abnormal_prob <- ifelse(sum(abnormal_n,normal_n)>0 & abnormal_n==0, .5*normal_prob, abnormal_prob)
      #normal_prob <- ifelse(sum(abnormal_n,normal_n)>0 & normal_n==0, .5*abnormal_prob, normal_prob)
      
      
      tot_prob_stage = abnormal_prob[,1]*normal_prob[,1]
      #tot_prob_stage= (abnormal+normal)/dim(bio_order)[1]
      
      
      p_perm_k[,i+1] <-   tot_prob_stage
      
    }
    
    if(i == max(group)){
      #special case stage 1
      abnormal <- bio_order
      
      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #abnormal_prob[which(abnormal_prob==0)] <- .1
      
      
      p_perm_k[,i+1] <- abnormal_prob[,1]
    }
    
    
    
    
  }
  #p_perm_k <-  ifelse(p_perm_k==0,y,p_perm_k)
  
  prob_subj = p_perm_k*(1/(max(group)+1))
  #prob_subj <- prob_subj+0.01
  total_prob_subj = apply(prob_subj,1,sum, na.rm=T)
  
  loglike = sum(log(total_prob_subj + 1e-250))
  
  
  return(list(p_perm_k, prob_subj, total_prob_subj, loglike))
}


fit_STATIS <- function(data,p_vec,clinical_info, nstart, initial_iter){
  # clinical_info=D
  colnames(clinical_info) <- c("clinical_measure", "event_number", "var_name")
  
  info <- clinical_info %>% dplyr::group_by(clinical_measure) %>% dplyr::summarise(events=dplyr::n())
  info <- data.frame(info) 
  colnames(info) <- c("bio", "events")
  #save likelihoods for prelim sequence search
  prelim_like <- vector(mode='list', length=nstart)
  
  #save sequences for prelim sequence search
  prelim_seq <- vector(mode='list', length=nstart)
  all_likes <- data.frame(start=NA, iter=NA,  like=NA)
  add <- possible_seqs(sum(info$events))
  N=dim(info)[1]
  max_like <- rep(NA, nstart)
  for(i in 1:nstart){
    
    all_seqs <- vector(mode='list', length=nstart)
    
    start_seq <- get_seq(clinical_info)
    
    start_group <- get_group(start_seq, add) 
    #start_group <- start_group %>% arrange(pos)
    current_likelihood <- get_likelihood(data, start_group,p_vec)[[4]]
    
    prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
    prelim_like_sub[1,1] <- current_likelihood
    prelim_like_sub[,2] <- i
    
    
    current_seq <- start_group
    
    
    for(j in 1:initial_iter){
      
      #print(j)
      
      bio_info_long <- dplyr::arrange(current_seq, pos)
      
      selected_pos = sample(x = 1:dim(bio_info_long)[1], size  = 1, replace=F)
      
      elig_pos <- sample(x = setdiff(c(1:dim(bio_info_long)[1]), c(selected_pos)), size  = dim(bio_info_long)[1]-1, replace=F)
      
      
      l=0
      for(k in elig_pos){
        l=l+1
        #k=6
        #print(k)
        bio_info_temp = bio_info_long
        possible_pos <- k
        temp <- bio_info_temp[k,5]
        bio_info_temp[k,5]<- bio_info_temp[selected_pos,5]
        bio_info_temp[selected_pos,5]<- temp
        #print(bio_info_temp)
        
        #check if eligible sequence based on events within each bio
        check <- bio_info_temp %>% dplyr::arrange(order) %>% dplyr::group_by(bio) %>% dplyr::summarize(Result = all(diff(event) == 1)) %>% dplyr::ungroup()
        #print(check)
        if(all(check$Result)){
          bio_info_long2 <- bio_info_temp %>% dplyr::arrange(order)
          break
        }
        if(l==length(elig_pos)){
          bio_info_long2 <- bio_info_long
        }
        
      }
      
      #start_group
      #bio_info_long2
      #get new group
      new <- get_group(bio_info_long2, add) #%>% dplyr::arrange(pos)
      
      #print(new$sub)
      all_seqs[[j]] <- new
      temp_likelihood = get_likelihood(data, new,p_vec)[[4]]
      
      
      prelim_like_sub[j,1] <- current_likelihood
      all_likes=rbind(all_likes, data.frame(start=i, iter=j, like=current_likelihood))
      #print(new)
      #if current sequence improves likelihood, update current likelihood and sequence
      if(temp_likelihood > current_likelihood){
        
        prelim_like_sub[j,1] <- temp_likelihood
        
        current_likelihood <-  temp_likelihood
        
        current_seq <- new
        
      }
      
      
    }
    #max_pos <- which.max(prelim_like_sub[,1])
    prelim_like[[i]] <- current_likelihood
    prelim_seq[[i]] <- current_seq
    max_like[i] <- current_likelihood
    
    
    
  }
  
  ml_seq <- prelim_seq[[which.max(max_like)]]
  
  #ml_seq <- cbind(ml_seq, event_name=data.frame(event_name=clinical_info$var_name))
  #ml_seq$est_seq <- ml_seq$sub
  #ml_seq <- ml_seq %>% dplyr::select(-pos,-order, -sub)
  
  all_likes <- na.omit(all_likes)
  #find maximum likelihood and corresponding sequence
  
  return(list(prelim_like=prelim_like, prelim_seq=prelim_seq, ml=max_like, ml_seq=ml_seq, loglikes=all_likes))
  
}



make_data <- function(D, p, M){
  seq <- get_seq(D)
  all_seqs  <- possible_seqs(dim(D)[1])
  N= dim(D)[1]
  
  
  group <- get_group(seq, all_seqs)
  
  #save true order and arrange data frame to be bio, event
  S <- group %>% dplyr::arrange(bio,event)
  true_order <- S$order
  
  
  #generate stage for individuals 1:M, randomly sample from all possible stages with replacement
  dat <- data.frame(Index=c(1:M),stage=sample(1:(max(S$sub)+1), M, replace=T))
  
  #empty matrix of M rows, K-1 columns to store event variables
  bio <- as.data.frame(matrix(nrow=M,ncol=max(S$pos)))
  
  #combine true stage and empty matrix to fill in
  dat2 <- cbind(dat, bio)
  
  check <- matrix(nrow=M, ncol=N+1)
  check[,1] <- dat2$stage
  #loop across 1:M
  for(i in 1:M){
    
    #special case when stage=1, when p=1, all events should = 0
    if(dat2$stage[i]==1){
      
      for(j in 1:ncol(bio)){
        x<- S$bio[j]
        n <- sum(S$bio==x)
        dat2[i,j+2] <- rbinom(1,1,(1-p))
        check[i,j+1] <- (1-p)
      }
      
    }
    
    
    if(dat2$stage[i]!=1){
      
      S_temp <- S %>% dplyr::filter(sub <= (dat2$stage[i])-1)
      to_fill <- dim(S_temp)[1]
      for(j in 1:(to_fill)){
        dat2[i,j+2] <- rbinom(1,1,p)
        check[i,j+1] <- p
      }
      
    }
    
    
    for(j in 1:ncol(bio)){
      x<- S$bio[j]
      n <- sum(S$bio==x)
      if(is.na(dat2[i,j+2])){
        dat2[i,j+2] <- rbinom(1,1,1-p)
        check[i,j+1] <- (1-p)
      }
      
    }
    
  }
  
  #fill in events that have not occurred
  
  dat2 <- dat2[,c(1,2,2+as.numeric(true_order))]
  colnames(dat2)[3:ncol(dat2)] <- paste0("V", 1:(ncol(dat2)-2))
  
  return(list(dat2, S))
  
}

