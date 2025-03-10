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
library(tidyr)

setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/TRUST")

### TRUST
#redcap data set file name
input_data <- "./data//TRUST_DATA_2024-02-21_0840.cleaned.wide_noreenroll.csv"

mlseq = read.csv("./results/trust_ml_seq.csv")
likes = read.csv("./results/trust_ml_likes.csv")
class = read.csv("./results/trust_class.csv")

#plot class dist
classdf=class %>% dplyr::select(pred_stage) %>% dplyr::group_by(pred_stage) %>% dplyr::summarise(n=n())
classdf$perc=classdf$n/sum(classdf$n)
dist=ggplot(data=classdf) + geom_bar(aes(x=pred_stage, y=perc), stat="identity", fill="#3A68AB") +  geom_text(aes(x=pred_stage, y=perc,label=paste0("N=",n)), vjust=-.1, size=5) + theme_bw() + xlab("Disease Class")+ scale_y_continuous(labels = scales::percent_format()) + ylab("Participants (%)") + scale_x_continuous(breaks=c(0:6))+ theme(text=element_text(size=14))
ggsave("./plots/TRUST dist.png",width=5, height=4 )

#plot likes
library(tbSTATIS)
plot_likes(likes) + theme(legend.position="none")
ggsave("./plots/TRUST likes.png",width=5, height=4 )


#outcomes
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

obs <- vars_bin %>% dplyr::select(-pid)

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

long <- vars_bin %>% melt()%>% dplyr::arrange(value)

#cough, weightloss,smear pos, cavity, sweat,infiltrates, smear 3, hemop

dem <- data %>% dplyr::filter(pid %in% vars_bin$pid)

relabel_binary.fn = function(x) {
  if(class(x) == "logical") {
    x <- factor(x, levels = c(TRUE, FALSE), labels = c("Yes", "No"))
  }else {
    x <- factor(x, levels = c(1,0), labels = c("Yes", "No"))
  }
  return(x)
}


dem$cough <- ifelse(dem$bl_cough > 0, 1, 0)
dem$cough[is.na(dem$cough)] <- 0

dem$sweat <- ifelse(dem$bl_sweat > 0, 1, 0)
dem$sweat[is.na(dem$sweat)] <- 0

dem$hemop <- ifelse(dem$bl_hemop > 0, 1, 0)
dem$hemop[is.na(dem$hemop)] <- 0

dem$fever <- ifelse(dem$bl_fever > 0, 1, 0)
dem$fever[is.na(dem$fever)] <- 0

dem$wtloss <- ifelse(dem$bl_wtloss > 0,1,0)
dem$wtloss[is.na(dem$wtloss)] <- 0


#cavity
dem$cavity <- ifelse(dem$cxr_cavity_chest_radiograph_1 ==1, 1, ifelse(dem$cxr_cavity_chest_radiograph_1 ==0, 0, NA))

dem$infiltrates <- ifelse(dem$cxr_infiltrate_chest_radiograph_1 %in% c(2), 1, ifelse(dem$cxr_infiltrate_chest_radiograph_1 %in% c(1, 0), 0, NA))
#cxr finding

#baseline demographics/clinical 
dem$screen_sex <- factor(dem$screen_sex, levels=c(1,2), labels=c( "Male", "Female"))
dem$education_baseline <- factor(dem$education_baseline, levels = c(T, F), labels = c("< grade 9", ">= grade 9"))
dem$bl_prison <- factor(dem$bl_prison, levels = c(1, 0), labels = c("Yes", "No"))
dem$prevtb_outcome <- factor(dem$prevtb_outcome, levels = c(1, 2, 3), labels = c("Cured", "Treatment Completed", "Treatment Defaulted"))
dem$bmi <- factor(dem$bmi, levels = c( "Underweight", "Normal Weight","Overweight and Obese"), labels=c( "Underweight", "Normal\nWeight","Overweight\nand Obese"))
dem$diabetes2 <- factor(dem$diabetes2, levels = c("Normal", "Pre-diabetes", "Diabetes"))
dem$mixed_ancestry_race <- factor(dem$mixed_ancestry_race, levels = c(1, 0), labels = c("Yes", "No"))
dem$screen_race <- factor(dem$screen_race, levels = c(1,2,3,4,5), labels = c("Mixed ancestry", "Black African", "White", "Indian/Asian", "Other"))
dem$unemployment_baseline <- factor(dem$unemployment_baseline, levels = c(T, F), labels=c("Unemployed", "Employed"))
dem$smoked_substance_use <- factor(dem$smoked_substance_use, levels = c(T, F), labels = c( "Smoked Substance Use", "No Smoked Substance Use"))
dem$age_cat <- factor(dem$age_cat, levels = c(1,2,3,4,5,6), labels = c("<30", "<30", "30-39", "40-49", "50+", "50+"))

dem$bl_hiv <- factor(dem$bl_hiv, levels=c(1,0),labels=c("Positive", "Negative"))

#these variables are TRUE/FALSE -- for table purposes, reassign as categorical
dem = dem %>% mutate_at(c("fstrom1_baseline",
                          "bl_amphetamine",
                          "dudit4d_baseline",
                          "dudit4f_baseline",
                          "dudit4g_baseline",
                          "bl_prevtb",
                          "bl_art", "cannabis_use", "meth_use", "mandrax_use", "bl_inh_monoresistant", "cough", "sweat", "hemop", "fever", "wtloss", "cavity", "infiltrates"), 
                        relabel_binary.fn)

dem$cxr_cavity_chest_radiograph_1 <- factor(dem$cxr_cavity_chest_radiograph_1, levels = c(1,0,2), c("Yes", "No", "Unknown"))
dem$smear_pos <-factor(dem$s_concafb_sputum_specimen_1, levels = c(0,4,1,2,3), labels = c("No AFB", "Scanty-++", "Scanty-++", "Scanty-++", "+++"))
dem$alc= factor(dem$tlfb_currdrink_baseline==1, levels=c(T,F), labels=c("Yes", "No"))


tableby(~ age_cat + screen_sex + bl_hiv + bmi + bl_prison  + unemployment_baseline + education_baseline + fstrom1_baseline + smoked_substance_use+ alc +bl_prevtb + cough + sweat + fever+ wtloss+ cavity+ infiltrates + smear_pos, data=dem) %>% summary()

dem$outcome<- factor(dem$to_programto_treatment_outcome, levels = c(1,2,3, 4, 5, 6, 7), labels = c("Cured/completed treatment", "Cured/completed treatment", "LTFU", "Treatment failure", "Died", "Moved/transferred", "Moved/transferred") )

dem$favorable=ifelse(dem$outcome %in% "Cured/completed treatment", "Favorable", "Unfavorable")
dem$favorable[dem$outcome=="Moved/transferred"] <- NA
data1 <- cbind(vars_bin, data.frame(stage=class$pred_stage))
data1 <- data1%>% dplyr::select(pid, stage)

data3 <- dplyr::left_join(dem,data1, by="pid")
data3$culture_conversion_sputum_specimen_9[data3$culture_conversion_sputum_specimen_9=="tb_positive_contaminated"]="tb_positive"
data3$culture_conversion_sputum_specimen_10[data3$culture_conversion_sputum_specimen_10=="tb_positive_contaminated"]="tb_positive"

prop.table(table(data3$stage, data3$culture_conversion_sputum_specimen_9))

data3$pos <- ifelse(data3$culture_conversion_sputum_specimen_9=="tb_positive", 1, ifelse(is.na(data3$culture_conversion_sputum_specimen_9)==F, 0, NA))

classdf=class %>% dplyr::select(pred_stage) %>% dplyr::group_by(pred_stage) %>% dplyr::summarise(n=n())

classdf$perc=classdf$n/sum(classdf$n)

stageplot=ggplot(data=classdf) + geom_bar(aes(x=pred_stage, y=perc), stat="identity", fill="#3A68AB") +  geom_text(aes(x=pred_stage, y=perc,label=paste0("N=",n)), vjust=-.1, size=4) + theme_bw()+ ylim(0,.35) + xlab("Disease Class")+ scale_y_continuous(labels = scales::percent_format()) + ylab("Participants (%)") + scale_x_continuous(breaks=c(0:5))+ theme(text=element_text(size=14)) +ggtitle("A)") 

log_model=glm(pos ~ factor(stage, levels=c(0,1,2,3,4,5,6,7), labels=c(1,1,2,3,4,5,6,7)) + age_cat + screen_sex + bl_hiv ,
    data = data3,
    family = "binomial")

tccors <- data.frame(or=(coefficients(log_model)), lb=(confint.default(log_model))[,1], ub=(confint.default(log_model))[,2])
tccors <- tccors[2:7,]
tccors$class <- c(2,3,4,5, 6, 7)

c=ggplot(data=tccors) + geom_point(aes(x=class, y=or)) + geom_errorbar(aes(x=class, ymin=lb, ymax=ub), width=.2) + theme_bw()+ ylab("log(Odds Ratio)") + xlab("Disease Class") + theme(text=element_text(size=14)) +ggtitle("B)")


ggarrange(stageplot, c, nrow=2)
ggsave("./plots/trust_res.png", width = 7, height=8, unit="in", dpi=500)    


#covariates and disease stage
#chi square tests
data3 %>% dplyr::group_by(stage, pos) %>% dplyr::summarise(n=n())

#plots 

#age, sex, hiv, smoked drug use, alcohol use, tobacco use, BMI
data3$smoked_substance_use <- factor(as.character(data3$smoked_substance_use), levels=c("Smoked Substance Use", "No Smoked Substance Use"), labels=c("Yes", "No"))

#age
age=data3 %>% dplyr::group_by(stage, age_cat) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
ageplot=ggplot(data=age) + geom_bar(aes(x=stage,y=prop, fill=age_cat), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7)) +scale_fill_manual(name="Age", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + theme(plot.margin = margin(0,0,0,5),legend.position = "top",text=element_text(size=16),legend.text=element_text(size=12))+ylim(0,1)

hiv=data3 %>% dplyr::group_by(stage, bl_hiv) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
hivplot=ggplot(data=hiv) + geom_bar(aes(x=stage,y=prop, fill=bl_hiv), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7)) +scale_fill_manual(name="HIV", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + theme(plot.margin = margin(0,0,0,5),legend.position = "top",text=element_text(size=16), legend.text=element_text(size=12))+ylim(0,1)
drugs=data3 %>% dplyr::group_by(stage, smoked_substance_use) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
drugsplot=ggplot(data=drugs) + geom_bar(aes(x=stage,y=prop, fill=smoked_substance_use), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7)) +scale_fill_manual(name="Smoked\nSubstance Use ", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + theme(plot.margin = margin(0,0,0,5),text=element_text(size=16),legend.position = "top", legend.text=element_text(size=12))+ylim(0,1)

sex=data3 %>% dplyr::group_by(stage, screen_sex) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
sexplot=ggplot(data=sex) + geom_bar(aes(x=stage,y=prop, fill=screen_sex), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency")+scale_fill_manual(name="Sex", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + scale_x_continuous(breaks=c(0:7))  + theme(plot.margin = margin(0,0,0,5),legend.position = "top",text=element_text(size=16), legend.text=element_text(size=12))+ylim(0,1)

bmi=data3 %>% dplyr::group_by(stage, bmi) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
bmiplot=ggplot(data=bmi) + geom_bar(aes(x=stage,y=prop, fill=bmi), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency")+scale_fill_manual(name="BMI", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + scale_x_continuous(breaks=c(0:7))  + theme(plot.margin = margin(0,0,0,5),legend.position = "top",legend.text=element_text(size=10), text=element_text(size=16))+ylim(0,1)

tob=data3 %>% dplyr::group_by(stage, fstrom1_baseline) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
tobplot=ggplot(data=tob) + geom_bar(aes(x=stage,y=prop, fill=fstrom1_baseline), stat="identity", position="dodge")+scale_fill_manual(name="Tobacco\nUse", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7))  + theme(plot.margin = margin(0,0,0,5),text=element_text(size=16),legend.position = "top", legend.text=element_text(size=12))+ylim(0,1)
data3$alc <- factor(data3$tlfb_currdrink_baseline, levels=c(1, 0), labels=c("Yes", "No"))
alc=data3 %>% dplyr::group_by(stage, alc) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
alcplot=ggplot(data=alc) + geom_bar(aes(x=stage,y=prop, fill=alc), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7))+ scale_fill_manual(name="Alcohol \nUse", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3"))  + theme(plot.margin = margin(0,0,0,5),legend.position = "top",text=element_text(size=16), legend.text=element_text(size=12))+ylim(0,1)

prevtb=data3 %>% dplyr::group_by(stage, bl_prevtb) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
prevtbplot=ggplot(data=prevtb) + geom_bar(aes(x=stage,y=prop, fill=bl_prevtb), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class")+scale_fill_manual(name="History \nof TB", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3")) + ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7))  + theme(plot.margin = margin(0,0,0,5),legend.position = "top",text=element_text(size=16), legend.text=element_text(size=12)) +ylim(0,1)

data3$hhs <- factor(data3$household_hunger_bin_baseline, levels=c(1,0), labels=c("Moderate\nto Severe", "Little\nto None"))
hhs <- data3 %>% dplyr::group_by(stage, hhs) %>% dplyr::summarise(n=n()) %>% na.omit() %>% dplyr::mutate(prop=n/sum(n))
hhsplot=ggplot(data=hhs) + geom_bar(aes(x=stage,y=prop, fill=hhs), stat="identity", position="dodge") + theme_bw() + xlab("Disease Class") + scale_fill_manual(name="Household \nHunger", values=c("#00008B", "#6495EB", "#B0C4DE", "#D3D3D3"))+ylab("Relative Frequency") + scale_x_continuous(breaks=c(0:7))  + theme(legend.position = "top",plot.margin = margin(0,0,0,5),text=element_text(size=16), legend.text = element_text(size=12)) +ylim(0,1)

q=ggarrange(ageplot, sexplot, hivplot, bmiplot, prevtbplot, hhsplot, alcplot, tobplot, drugsplot, nrow=3, ncol=3)

ggsave("./plots/barplots.png", width=14,height=14)


#make pvd

#main analysis
setwd("./results")
filenames <- list.files(pattern = 'TRUST_boot_ml.*.csv', recursive = TRUE)
boot <- purrr::map_df(filenames, read.csv, .id = 'id', header=F)
ex <- seq(1,300,2)*-1
boot<- boot[ex,]
boot <- boot[1:100,]
print(mlseq)

#measures <- vars_bin %>% dplyr::select(cavity, wtloss, sweat, smear_pos, smear_3, unilatinfiltrates, infiltrates, cough, fever)

pvd <- data.frame(wtloss=c(99,1,0,0,0,0,0,0,0),cough=c(99,1,0,0,0,0,0,0,0), unilatinfiltrates=c(100, 0, 0, 0, 0, 0, 0, 0,0) ,smearpos=c(1,79,10,10,0,0,0,0,0),
                  cavity=c(0,13,83, 4, 0, 0, 0, 0,0), bilat=c(0,13,12,74, 1, 0, 0, 0,0), sweat=c(36,0,0,7,56,1,0,0,0),fever=c(0,0,0,5,38,56,7,0,0), smear3=c(0,0,0,0,5,38,56,1,0))
pvd<-pvd/100
pvd$rowname=rownames(pvd)

#wtloss
table(boot$V2)
#cough
table(boot$V8)
#uni
table(boot$V6)
#smearpos
table(boot$V4)
#cav
table(boot$V1)

#bilat
table(boot$V7)
#sweats
table(boot$V3)
#fever
table(boot$V9)
#smear3
table(boot$V5)

dt2 <- pvd %>%
  gather(colname, value, -rowname)
dt2$colname <- factor(dt2$colname, levels=c("wtloss", "cough","unilatinfiltrates","smearpos", "cavity", "bilat",   "sweat", "fever", "smear3"), labels=c("Unexplained\nweight loss","Cough", "Unilateral\ninfiltrates", "Smear:\nScanty-2+", "Cavity", "Bilateral\n infiltrates", "Night sweats", "Fever", "Smear: 3+" ))

d <- data.frame(Y=unique(dt2$colname), X=mlseq$sub)

pvd=ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() + xlab("Disease Class") + ylab("Clinical State")+ geom_tile(data=d, aes(X,Y), fill="transparent", colour="black", linewidth=1.2) + scale_fill_gradient(name="Proportion of\nbootstrap samples   ",low = "white", high = "#3A68AB") + theme_bw() + theme(text=element_text(size=16, color="black"), axis.text = element_text(color="black"), legend.position="top", legend.key.size = unit(.75, 'cm'), legend.title = element_text(size=10) , legend.text = element_text(size=10))

ggsave("../plots/pvd_trust.png", height = 5, width = 5)
