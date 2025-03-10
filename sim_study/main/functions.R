##################################################
## Project:tbSTATISpaper1
## Script purpose: This program contains all functions to run TB-STATIS and generate data for a simulation study.
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################

library(dplyr)

#' Function: get_stage
#' Purpose: Obtain predicted disease class for individuals after estimating TB-STATIS
#'

get_stage <- function(data=data.frame(), S=data.frame(), p=vector()){
  
  #get likelhood 
  likelihood=get_likelihood(data, S, p )
  
  #save class probs
  stage_probs <- likelihood[[2]]
  stage_probs=data.frame(stage_probs)
  #update colnames to class numbers
  colnames(stage_probs) = c(0:max(S$sub))
  
  #assign disease class with max prob
  pred_stage=as.numeric(colnames(stage_probs)[apply(stage_probs,1,which.max)])
  
  #cbind to original data
  stagedf <- cbind(data,pred_stage)
  
  return(stagedf)
}

#' Function: make_D
#' Purpose: Generate set of clinical measures and # of states per measure given N total clinical states 
#'
make_D <- function(N){
 
  #select number of clinical measures
  nbio <- sample(2:N, 1)
  D <- data.frame(bio=paste0("bio", 1:nbio), events=rep(1, nbio))
  
  #generate number of states per measure
  repeat{
    if(sum(D$events)==N){
      break
    }
    rowadd <- sample(1:nbio, 1)
    D$events[rowadd] <- D$events[rowadd]+1
  }
  
  return(D)
}

#' Function: get_group
#' Purpose: Generate simulataneous disease sequence given sequence of clinical states. 
#'
get_group <- function( info=data.frame(), all_seqs){
  #arrange sequence of states
  info <- info %>% dplyr::arrange(order)
  N=dim(info)[1]
  
  #select disease sequence allowing multiple states per disease class
  #code allows for multiple states in one measure to occur at the same class
  #uncomment code below if clinical measures should be unique per class
  all_seqs2=all_seqs
  all_seqs=data.frame(t(apply(all_seqs[,1:N], 1, function(i) paste(i, info$bio) )))
  all_seqs$unique= apply(all_seqs[,1:N], 1, function(x) length(unique(x)))==N
  
  #all possible seqs after criteria above
  final <-  cbind(all_seqs2, data.frame(unique=all_seqs$unique))#%>% dplyr::filter(unique==T)
  #randomly select one sequence to use
  select <- sample(1:dim(final)[1],1)
  
  #add to info df
  info$sub<- as.numeric(final[select,1:N])

  # print(info)
  return(info)
}

#' Function: get_seq
#' Purpose: Propose sequence as start point for maximum likelihood estimation.
#'

get_seq <- function(info=data.frame()){
  #save info df to work with
  bio_info_long <- info
  #update col names
  colnames(bio_info_long) <- c("bio", "event", "var_name")

  N=dim(bio_info_long)[1]
  #save original order of states
  bio_info_long$pos = c(1:N)
  
  #sample from clinical measures w/o replacement as new order and save as temp df
  temp<-sample(bio_info_long$bio, N)
  
  #update order in df
  bio_info_long_temp <- bio_info_long
  bio_info_long_temp$temp <- temp
  bio_info_long_temp <- bio_info_long_temp %>% dplyr::arrange(temp) 
  bio_info_long_temp$order <- c(1:dim(bio_info_long_temp)[1])
  #bio_info_long_temp
  bio_info_long_new <- bio_info_long #%>% dplyr::arrange(pos) %>% dplyr::select(-temp)
  bio_info_long_new$order <- bio_info_long_temp$pos
  
  return(bio_info_long_new)
}




#' Function: possible_seqs
#' Purpose: Generate data frame of all possible sequences given N total clinical states
#'

possible_seqs <- function(N=numeric()){
  #filter such that first disease class is 1 and are monotonically non-decreasing
  filter_rows <- function(g, x) {
    ok  <- function(z) all(diff(z) %in% 0:1)
    out <- g[apply(g, 1, ok), ]
    replace(out, TRUE, lapply(out, \(i) x[i]))
  }
  #generate all seqs
  f2 <- function(x = c(1:N), n=N+1, n1=2) {
    data.frame(as.list(rep(1, n1)),
               gtools::combinations(length(x), n-n1, repeats.allowed = TRUE)) |>
      filter_rows(x)
  }
  #generate all possible sequences of size N and filter using filter_rows and f2 functions
  all=f2() 
  all <- all[,-1]
  return(all)
}


#' Function: get_likelihood
#' Purpose: Calculate likelhood given observed data, sequence S, vector of clinical measure accuracies p_vec
#'
get_likelihood <- function(data=data.frame(), S=data.frame(), p_vec=vector()){
  #input data
  bio <- data.frame(data[,1:ncol(data)])
  
  # number of biomarkers
  N = as.numeric(dim(bio)[2])
  
  #individuals in data set
  M = dim(bio)[1]
  
  #new df to become probabilities from 0/1 data
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
  group = as.numeric(S$sub)
  
  #reorder biomarker columns based on S
  bio_order = new_bio
  
  #bio_order=bio_order+alpha
  normal <- 1-bio_order
  
  #special case class 1
  normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
  
  #save probs for class 1
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
    
      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)

      normal_prob = normal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
    
      tot_prob_stage = abnormal_prob[,1]*normal_prob[,1]
      #tot_prob_stage= (abnormal+normal)/dim(bio_order)[1]
      
      #save probs
      p_perm_k[,i+1] <-   tot_prob_stage
      
    }
    #special case max class, all states should have occurred
    if(i == max(group)){
      #special case stage 1
      abnormal <- bio_order
      
      abnormal_prob = abnormal %>% dplyr::mutate(prod = apply(., 1, prod, na.rm=T)) %>% dplyr::select(prod)
      #abnormal_prob[which(abnormal_prob==0)] <- .1
      
      
      p_perm_k[,i+1] <- abnormal_prob[,1]
    }
    
    
    
    
  }
  #divide probs by total classess
  prob_subj = p_perm_k*(1/(max(group)+1))
  
  #multiply across rows to get individual likelihood for class k
  total_prob_subj = apply(prob_subj,1,sum, na.rm=T)
  
  #full likelhood by summing log likes for each ind. 
  loglike = sum(log(total_prob_subj + 1e-250))
  
  
  return(list(p_perm_k, prob_subj, total_prob_subj, loglike))
}

#' Function: fit_STATIS
#' Purpose: Estimate TB-STATIS
#'
fit_STATIS <- function(data,p_vec,clinical_info, nstart, initial_iter){

  colnames(clinical_info) <- c("clinical_measure", "event_number", "var_name")
  
  info <- clinical_info %>% dplyr::group_by(clinical_measure) %>% dplyr::summarise(events=dplyr::n())
  info <- data.frame(info)
  colnames(info) <- c("bio", "events")
  
  
  #save likelihoods for prelim sequence search
  prelim_like <- vector(mode='list', length=nstart)
  
  #save sequences for prelim sequence search
  prelim_seq <- vector(mode='list', length=nstart)
  all_likes <- data.frame(start=NA, iter=NA,  like=NA)
  #df of all possible seqs to sample from during likelihood ascent
  add <- possible_seqs(sum(info$events))
  N=dim(info)[1]
  max_like <- rep(NA, nstart)
  
  #for each nstart, generate an initial sequence then search for mlseq 
  for(i in 1:nstart){
    
    all_seqs <- vector(mode='list', length=nstart)
    #get initial seq to start search
    start_seq <- get_seq(clinical_info)
    
    start_group <- get_group(start_seq, add) 
    #calc first likelihood 
    current_likelihood <- get_likelihood(data, start_group,p_vec)[[4]]
    
    prelim_like_sub <- matrix(nrow=initial_iter, ncol=2)
    prelim_like_sub[1,1] <- current_likelihood
    prelim_like_sub[,2] <- i
    
    
    current_seq <- start_group
    
    #for each iteration generate new seq, calculate likelihood and compare to current likelihood, update if new likelihood is greater
    for(j in 1:initial_iter){
      
      #print(j)
      
      bio_info_long <- dplyr::arrange(current_seq, pos)
      #select state to move
      selected_pos = sample(x = 1:dim(bio_info_long)[1], size  = 1, replace=F)
      #elig states to swap with
      elig_pos <- sample(x = setdiff(c(1:dim(bio_info_long)[1]), c(selected_pos)), size  = dim(bio_info_long)[1]-1, replace=F)
      
      #loop through eligible states, select first one that works
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
      

      #get new group after swapping states
      new <- get_group(bio_info_long2, add) #%>% dplyr::arrange(pos)
      
      #get likelihood 
      all_seqs[[j]] <- new
      temp_likelihood = get_likelihood(data, new,p_vec)[[4]]
      
      
      prelim_like_sub[j,1] <- current_likelihood
      all_likes=rbind(all_likes, data.frame(start=i, iter=j, like=current_likelihood))
 
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
  
  #save ml seq
  ml_seq <- prelim_seq[[which.max(max_like)]]
  
  #save all likes to plot ascent later
  all_likes <- na.omit(all_likes)

  return(list(prelim_like=prelim_like, prelim_seq=prelim_seq, ml=max_like, ml_seq=ml_seq, loglikes=all_likes))
  
}


#' Function: make_data
#' Purpose: Generate simulated data set with sample size M, clinical measures in D, and clinical measure accuracies in p.
#'
make_data <- function(D, p, M){
  #get true sequence to generate data from
  seq <- get_seq(D)
  all_seqs  <- possible_seqs(dim(D)[1])
  N= dim(D)[1]
  
  #get true simultaneous sequence 
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
    
    #special case when class=1, when p=1, all events should = 0
    if(dat2$stage[i]==1){
      
      for(j in 1:ncol(bio)){
        x<- S$bio[j]
        n <- sum(S$bio==x)
        dat2[i,j+2] <- rbinom(1,1,(1-p))
        check[i,j+1] <- (1-p)
      }
      
    }
    
    #fill in states that have occurred
    if(dat2$stage[i]!=1){
      
      S_temp <- S %>% dplyr::filter(sub <= (dat2$stage[i])-1)
      to_fill <- dim(S_temp)[1]
      for(j in 1:(to_fill)){
        dat2[i,j+2] <- rbinom(1,1,p)
        check[i,j+1] <- p
      }
      
    }
    
    #fill in events that have not occurred
    for(j in 1:ncol(bio)){
      x<- S$bio[j]
      n <- sum(S$bio==x)
      if(is.na(dat2[i,j+2])){
        dat2[i,j+2] <- rbinom(1,1,1-p)
        check[i,j+1] <- (1-p)
      }
      
    }
    
  }
  

  #order columns and update col names
  dat2 <- dat2[,c(1,2,2+as.numeric(true_order))]
  colnames(dat2)[3:ncol(dat2)] <- paste0("V", 1:(ncol(dat2)-2))
  
  return(list(dat2, S))
  
}

