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
library(tbSTATIS)

setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/REMox")

source("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/functions.R")


#data
data <- read.csv("./data/REMox_data.csv")

mlseq <- read.csv("./results/remox_ml_seq.csv")

mllikes <- read.csv("./results/remox_ml_likes.csv")

class <- read.csv("./results/remox_class.csv")



plot_likes(mllikes) + theme(legend.position="none")

ggsave("REMox_likes.png", height=5, width=7, units="in", dpi=500)

obs <- data %>% dplyr::select( CAV, SMEAR1, SMEAR2, SMEAR3, SMEAR4,SWEATS, WTLOSS, COUGH, FEVER)

#make expected data given ml seq and predicted disease class
exp <- data.frame(matrix(ncol=dim(obs)[2], nrow=dim(obs)[1]))
colnames(exp) <- colnames(obs)
exp$class <- class$pred_stage

for(i in 1:dim(exp)[1]){
  if(exp$class[i]==0){
    exp[i,1:(ncol(exp)-1)] <- 0
  }
  
  if(exp$class[i]==max(mlseq$sub)){
    exp[i,1:(ncol(exp)-1)] <- 1
  }
  
  if(!(exp$class[i] %in% c(0, max(mlseq$sub)))){
    abnormal <- mlseq$var_name[mlseq$sub <= exp$class[i]]
    normal <-  mlseq$var_name[mlseq$sub > exp$class[i]]
    exp[i,normal] <- 0
    exp[i,abnormal] <- 1
  }
}

dh <- data.frame(matrix(ncol=dim(obs)[2], nrow=dim(obs)[1]))
colnames(dh) <- colnames(dh)
for(j in 1:ncol(dh)){
  dh[,j] <- obs[,j] != exp[,j]
  
}

#dh for each clinical state
dhi <- data.frame(dhi=apply(dh, 2, sum)/nrow(obs))
rownames(dhi) <- colnames(obs)

#overall fit
dhbar <- mean(dhi$dhi)


classdf=class %>% dplyr::select(pred_stage) %>% dplyr::group_by(pred_stage) %>% dplyr::summarise(n=n())

classdf$perc=classdf$n/sum(classdf$n)

stageplot=ggplot(data=classdf) + geom_bar(aes(x=pred_stage, y=perc), stat="identity", fill="#3A68AB") +  geom_text(aes(x=pred_stage, y=perc,label=paste0("N=",n)), vjust=-.1, size=4) + theme_bw()+ ylim(0,.35) + xlab("Disease Class")+ scale_y_continuous(labels = scales::percent_format()) + ylab("Participants (%)") + scale_x_continuous(breaks=c(0:5))+ theme(text=element_text(size=14)) +ggtitle("A)") 

new <- cbind(data, pred_class = class$pred_stage) #%>% dplyr::filter(ARM=="2EHRZ/4HR")


prop.table(table(new$pred_class, new$PRIME),1)

prop.table(table(new$pred_class, new$TCC<.15), margin=1)



#association with outcomes
cult <- new 
cult$status <- 1
cult$time <- as.numeric(cult$TCC)*52
cult$time2 <- ifelse(cult$time >=8, 8, cult$time)
cult$status2 <- ifelse(cult$time >=8, 0, cult$status)

log_model = glm(1-status2 ~ factor(pred_class) + AGECAT +SEX+HIV + ARM, data = cult, family = "binomial")
summary(log_model)
exp(coefficients(log_model))
exp(confint.default(log_model))
tccors <- data.frame(or=(coefficients(log_model)), lb=(confint.default(log_model))[,1], ub=(confint.default(log_model))[,2])
tccors <- tccors[2:6,]
tccors$class <- c(1,2,3,4,5)
tccors$out="CC"



cult <- new 
cult$status <-ifelse(cult$PRIME=="Favorable", 0, ifelse(cult$PRIME == "Unfavorable", 1, NA))
cult$time <- as.numeric(cult$TIME_PRIME)*52
cult$time2 <- ifelse(cult$time >=78, 78, cult$time)
cult$status2 <- ifelse(cult$time >= 78, 0, cult$status)

log_model = glm(status ~ factor(pred_class) + AGECAT +SEX+HIV + ARM, data = cult, family = "binomial")
summary(log_model)
exp(coefficients(log_model))
exp(confint.default(log_model))

favors <- data.frame(or=(coefficients(log_model)), lb=(confint.default(log_model))[,1], ub=(confint.default(log_model))[,2])
favors <- favors[2:6,]
favors$class <- c(1,2,3,4,5)
favors$out="FAV"

f=rbind(tccors, favors)
stageplot
f$out<- factor(f$out,levels=c("CC", "FAV"), labels=c("Outcome: culture positive at week 8", "Outcome: unfavorable TB outcome"))
c=ggplot(data=f) + geom_point(aes(x=class, y=or)) + geom_errorbar(aes(x=class, ymin=lb, ymax=ub), width=.2) + facet_wrap(~out, scales="free") + theme_bw()+ ylab("log(Odds Ratio)") + xlab("Disease Class") + theme(text=element_text(size=14)) +ggtitle("B)")

ggarrange(stageplot, c, nrow=2)
ggsave("./plots/remox_res.png", width = 7, height=8, unit="in", dpi=500)    


#make pvd
setwd("./results")
filenames <- list.files(pattern = 'REMox_boot_ml.*.csv', recursive = TRUE)


boot <- purrr::map_df(filenames, read.csv, .id = 'id', header=F)
ex <- seq(1,300,2)*-1
boot<- boot[ex,]

print(mlseq)

#obs <- data %>% dplyr::select( CAV, SMEAR1, SMEAR2, SMEAR3, SMEAR4,SWEATS, WTLOSS, COUGH, FEVER)

pvd <- data.frame(WTLOSS=c(51,13,0,35,0,0,0, 0, 0),FEVER=c(85,12,1,0,0,0,1, 0, 0),SMEAR1=c(99, 0, 0, 0, 0, 0, 0, 0, 0) ,COUGH=c(98,0,0,1,0,0,0, 0, 0),
                  CAV=c(98,0, 0, 0, 1, 0, 0, 0, 0), SMEAR2=c(12,87,0, 0, 0, 0, 0, 0, 0), SMEAR3=c(0,1,98,0,0,0,0, 0, 0), SMEAR4=c(0,0,0,63,35,1,0, 0, 0), SWEATS=c(1,12,0,35,50,0,1, 0, 0))
pvd<-pvd/100
pvd$rowname=rownames(pvd)

#wtloss
table(boot$V7)
#fever
table(boot$V9)
#smear1
table(boot$V2)
#cough
table(boot$V8)
#cav
table(boot$V1)
#smear2
table(boot$V3)
#smear3
table(boot$V4)
#smear4
table(boot$V5)
#sweats
table(boot$V6)


dt2 <- pvd %>%
  gather(colname, value, -rowname)
dt2$colname <- factor(dt2$colname, levels=c("WTLOSS", "FEVER","SMEAR1","COUGH", "CAV", "SMEAR2",   "SMEAR3", "SMEAR4", "SWEATS"), labels=c("Unexplained\nweight loss","Fever", "Smear: 1+", "Cough", "Cavity", "Smear: 2+", "Smear: 3+", "Smear: 4+", "Night sweats" ))

d <- data.frame(Y=unique(dt2$colname), X=mlseq$sub)

pvd=ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() + xlab("Disease Class") + ylab("Clinical State")+ geom_tile(data=d, aes(X,Y), fill="transparent", colour="black", linewidth=1.2) + scale_fill_gradient(name="Proportion of\nbootstrap samples   ",low = "white", high = "#3A68AB") + theme_bw() + theme(text=element_text(size=16, color="black"), axis.text = element_text(color="black"), legend.position="top", legend.key.size = unit(.75, 'cm'), legend.title = element_text(size=10) , legend.text = element_text(size=10))

ggsave("../plots/pvd_remox.png", height = 5, width = 5)
