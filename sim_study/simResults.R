##################################################
## Project:tbSTATISpaper1
## Script purpose: This program contains functions to summarize the data and sequences produced in simulation.   
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################

library(dplyr)
library(ggplot2)
library(ggpubr)

#function to check if ml sequence and estimated sequence are identical, generate plots 
same_seq <- function(ml, N, label){
#ml=seqs1
#N=4
#label=NA

ml=ml %>% na.omit()

same <- rep(NA, dim(ml)[1])

#compare ml seq and true seq by looping over columns
for(i in 1:dim(ml)[1]){
  
  x <- ml[i,3:(3+N-1)]
  colnames(x) <- c(1:N)
  y <- ml[i,(3+N):(3+2*N-1)]
  colnames(y) <- c(1:N)
  same[i] = identical(x,y)
  
}

dat <- cbind(ml, same)

#summarize into group proportions
same_seq <-dat %>%  dplyr::group_by(M,p, same) %>% dplyr::summarise(n=n()) %>% dplyr:::filter(same==T)
denom <- dat %>%  dplyr::group_by(M,p) %>% dplyr::summarise(denom=n()) 
same_seq <- cbind(same_seq, denom=denom$denom)
same_seq$prop=same_seq$n/same_seq$denom

#make bar plot 
p1=ggplot(data=same_seq) + geom_bar(aes(x=factor(p), y=prop, fill=factor(M)), stat = "identity", position="dodge")+ 
  ggtitle(paste0(label,") ",N," Clinical States"))+xlab("p") +
  scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + ylim(0,1) + 
  theme_bw() + ylab("Proportion\nCorrect Sequence")  + theme(legend.position = "bottom", text=element_text(size=16, color="black"), legend.text = element_text(size=16), axis.text = element_text(color="black"))

return(list(same_seq, p1))
}


#generate table for proportion correct sequence
proptab <- function(data){
  t1=ggplot(data) + facet_grid(~p+M, labeller =label_both, scales="free_y")+geom_text(data,mapping=aes(y=factor(0), x=factor(0), label=paste0(round(prop, digits=2)), hjust=.7, vjust=0.3))+
    ggtitle("Proportion Kendall's Tau Distance Equals 0") + theme_pubr()  +theme(
      axis.line = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.y  = element_blank(),
      axis.title.x  = element_blank(),
      axis.text.x = element_text(color="white"),
      axis.text.y = element_text(color="white"),
      plot.title = element_text(size=14)
    ) +theme(
      strip.text.x = element_text(
        size = 12, color = "black", face = "bold"
      ),strip.background = element_rect(
        color="black", fill="lightgrey", size=1, linetype="solid"
      ),panel.spacing = unit(0,'lines'),panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +theme(plot.margin = margin(0,0,0,1, "cm"))
  
  return(t1)
}


#function to calculate Kendall Tau Distance 
seq_dist <- function(data, N, label){
all_combos = data.frame(t(combn(N,2)))

#data=seq4
data <- data %>% na.omit()

#1 if a < b
#2 if a > b
#3 if a = b

#get diff vector for estimated sequences
data <- na.omit(data)

all_props <- vector(mode="numeric", length=dim(data)[1])

ml<- (data[,3:(3+N-1)])
t <- (data[,(3+N):(3+2*N-1)])

for(j in 1:dim(data)[1]){
  est=vector(mode="numeric", length = dim(all_combos)[1])
  true=vector(mode="numeric", length = dim(all_combos)[1])
  
  for(i in 1:length(est)){
    x=all_combos[i,1]
    y=all_combos[i,2]

    a <- t[j,x]
    b <- t[j,y]
    est[i] <- ifelse(a < b, 1, ifelse(a>b, 2, 3))
    
    a <- ml[j,x]
    b <- ml[j,y]
    true[i] <- ifelse(a < b, 1, ifelse(a>b, 2, 3))
  }
  
  #get diff vector for true sequences
  diff = data.frame(est, true)
  diff$disc <- diff$true!=diff$est
  
  prop = sum(diff$disc) /dim(diff)[1]
  all_props[j] <- prop
}

props_info <- data %>% dplyr::select(p,M)
props_info$dist=all_props

propsinfo_no0 <- props_info %>% dplyr::filter(dist!=0)

propsinfo0 <- props_info %>%  filter(dist==0)

max=props_info  %>% dplyr::group_by(p, M) %>% summarise( 
  n=n())

propsinfo0_summary <- propsinfo0 %>% dplyr::group_by(p, M) %>% summarise( 
  n=n()) 

propsinfo0_summary$denom =max$n
propsinfo0_summary$prop <- propsinfo0_summary$n/propsinfo0_summary$denom

props_summary <- propsinfo_no0 %>% dplyr::group_by(p, M) %>% summarise( 
  n=n(),
  mean=mean(dist),
  sd=sd(dist)
) %>%
  mutate( se=sd/sqrt(n))

closeto1 <- propsinfo0_summary %>% filter(prop >.995)

propsinfo_no0<- propsinfo_no0 %>% filter(!(paste0(M,p) %in% paste0(closeto1$M, closeto1$p)))

for(k in 1:dim(propsinfo0_summary)[1]){
if(propsinfo0_summary$prop[k]>.994){
 add <- propsinfo0_summary[k,1:2]
 add$dist=0
 propsinfo_no0=rbind(propsinfo_no0, add)
}
}

p1=ggplot(data = props_summary, 
          aes(x=factor(p),
              y= mean, 
              ymin=mean, 
              ymax=mean+se,       
              fill=factor(M))) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar( position = position_dodge(), colour="black") + ylim(0,1) + 
  scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + theme_bw() + ylab("Average Kendall's\nTau Distance") + xlab("p") +
  theme(legend.position = "top", text=element_text(size=14, color="black"), legend.text = element_text(size=14), axis.text = element_text(color="black")) + ggtitle("")


p2= ggplot(propsinfo_no0) + geom_boxplot(aes(y=dist, x=factor(p), fill=factor(M)), outlier.shape=NA) + ylab("Kendall's Tau Distance") +
  ggtitle(paste0(label,") ",N," Clinical States"))+xlab("p") + scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + 
  theme_bw()   + theme(legend.position = "bottom", text=element_text(size=16, color="black"), legend.text = element_text(size=16), axis.text = element_text(color="black")) +ylim(0,1)

p3= ggplot(propsinfo_no0) + geom_violin(aes(y=dist, x=factor(p), fill=factor(M))) + ylab("Kendall's Tau Distance") +
  ggtitle(paste0(label,") ",N," Clinical States"))+xlab("p") + scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + 
  theme_bw()   + theme(legend.position = "bottom", text=element_text(size=16, color="black"), legend.text = element_text(size=16), axis.text = element_text(color="black")) +ylim(0,1)

p4= ggplot(propsinfo_no0) + geom_density(aes( x=dist), fill="lightblue") + facet_grid(~p+M, labeller =label_both) + coord_flip() +xlab("Kendall's Tau Distance > 0") +
  ggtitle(paste0("   ",label,") ",N," Clinical States")) + ylab("Density") +
  theme_bw() + xlim(0,1) + theme(text=element_text(size=14, color="black")) +theme(
    strip.text.x = element_text(
      size = 12, color = "black", face = "bold"
    ),strip.background = element_rect(
      color="black", fill="lightgrey", size=1, linetype="solid"
    ),panel.spacing = unit(0,'lines'),panel.border = element_rect(color = "black", fill = NA, size = 1), plot.title = element_text(hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank()
  )

return(list(props_summary,p1, p2, p3, p4))
}



#function to compare predicted and true disease class, generate box plot 
class_diff <- function(data, N, label){
#data=dat4
data <- data %>% na.omit()

data$diff = (data$stage-data$pred_stage)
data$diffn <- abs(data$diff/(N+1))

#data=dat4
kendall.corr =data %>% dplyr::group_by(id,M,p) %>% dplyr::summarise(cor=cor(x=stage, y=pred_stage, method="kendall"))

#kendall.corr.mean = kendall.corr %>% dplyr::group_by(M, p) %>% dplyr::summarise(kt.avg=mean(cor))
p1= ggplot(kendall.corr, aes(y=cor, x=factor(p), fill=factor(M))) + geom_boxplot() + ylab("Kendall's Tau Correlation") +
  ggtitle(paste0(label,") ",N," Clinical States"))+xlab("p") + scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + 
  theme_bw()   + theme(legend.position = "top", text=element_text(size=14, color="black"), legend.text = element_text(size=14), axis.text = element_text(color="black"))+ylim(0,1)
  
p2= ggplot(data, aes(y=diffn, x=factor(p), fill=factor(M))) + geom_boxplot() + ylab("Distance Between\nTrue and Predicted Class") +
 ggtitle("")+xlab("p") + scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + 
 theme_bw()   + theme(legend.position = "bottom", text=element_text(size=16, color="black"), legend.text = element_text(size=16), axis.text = element_text(color="black"))+ ylim(0,1)

#p2= ggplot(data, aes(y=diffn, x=factor(p), fill=factor(M))) + geom_violin(adjust=40, width=1) + ylab("Average Distance Between\nTrue and Predicted Class") +
 # ggtitle(paste0(label,") ",N," Clinical States"))+xlab("Probability Observed Value\nEquals True Value") + scale_fill_manual(name="Sample size", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + 
 # theme_bw()   + theme(legend.position = "bottom", text=element_text(size=16, color="black"), legend.text = element_text(size=16), axis.text = element_text(color="black"))+ylim(0,1)

return(list(data, p1, p2))
}

