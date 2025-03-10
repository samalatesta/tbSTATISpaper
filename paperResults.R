##################################################
## Project:tbSTATISpaper1
## Script purpose: This program generates all plots to summarize results from the simulation study. 
## Date: 03/05/2025
## Author: Samantha Malatesta
##################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

#main simulation results
setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/main/sim_results")

source("../../simResults.R")

#read in simulated data
seq4 <- read.csv("res_seqs_4_2025-03-07.csv") 

dat4 <- read.csv("res_data_4_2025-03-07.csv")

seq8 <- read.csv("res_seqs_8_2025-03-07.csv")

dat8 <- read.csv("res_data_8_2025-03-07.csv")

seq12 <- read.csv("res_seqs_12_2025-03-08.csv")

dat12 <- read.csv("res_data_12_2025-03-08.csv")

#plot main simulation results 

#plot proportion correctly estimated sequences
x1 <- same_seq(seq4, 4, "A")
x2 <- same_seq(seq8, 8, "B")
x3 <- same_seq(seq12, 12, "C")

#arrange plots and save
xplot=ggarrange(x1[[2]], x2[[2]],x3[[2]], common.legend = T, nrow=3)
ggsave(paste0("./plots/same_seq_", Sys.Date(), ".png"), xplot, unit="in", width=8, height=9, res=500)

#plot difference between predicted and true disease class
y1 <- class_diff(dat4, 4, "A")[[2]]
y2 <- class_diff(dat8, 8, "B")[[2]]
y3 <- class_diff(dat12, 12, "C")[[2]]

#arrange plots and save
ggarrange(y1, y2,y3, common.legend = T)
yplot <- ggarrange(y1, y2,y3, common.legend = T, nrow=3)
ggsave(paste0("./plots/class_res_", Sys.Date(), ".png"), yplot, unit="in", width=7, height=14, dpi=500)

#plot Kendall Tau distance
z1 <- seq_dist(seq4, 4, "A")
z2 <- seq_dist(seq8, 8, "B")
z3 <- seq_dist(seq12, 12, "C")

#arrange plots and save
ggarrange(z1[[2]], z2[[2]], z3[[2]], common.legend = T, ncol=3)
zplot <- ggarrange(z1[[3]], z2[[3]], z3[[3]], common.legend = T, ncol=1)
ggsave(paste0("./plots/seq_dist_", Sys.Date(), ".png"), zplot)

#plot Kendall tau distance and add table for proportion correct seq below
f1=plot_grid(z1[[5]],proptab(x1[[1]]),nrow=2, rel_heights = c(2.5, 1))
f2=plot_grid(z2[[5]],proptab(x2[[1]]),nrow=2, rel_heights = c(2.5, 1))
f3=plot_grid(z3[[5]],proptab(x3[[1]]),nrow=2, rel_heights = c(2.5, 1))

f=ggarrange(f1,f2,f3,  nrow=2, ncol=2)

png(paste0("./plots/seq_res_", Sys.Date(), ".png"), unit="in", width=12, height=10, res=500)
f
dev.off()



#sensitivity analysis

#simulation results
setwd("/usr3/graduate/samalate/projectnb/cbs/samalate/tbSTATISpaper/sim_study/sensitivity/sim_results")

seq4sens <- read.csv("res_seqs_4_sens_2025-03-08.csv") 

dat4sens <- read.csv("res_data_4_sens_2025-03-08.csv")

seq8sens <- read.csv("res_seqs_8_sens_2025-03-07.csv")

dat8sens <- dat2#read.csv("res_data_8_sens_2025-03-09.csv")

seq12sens <- read.csv("res_seqs_12_sens_2025-03-09.csv")

dat12sens <- read.csv("res_data_12_sens_2025-03-09.csv")

#plot main simulation results 

#plot proportion correctly estimated sequences
x1 <- same_seq(seq4sens, 4, "A")
x2 <- same_seq(seq8sens, 8, "B")
x3 <- same_seq(seq12sens, 12, "C")

#arrange plots and save
xplot=ggarrange(x1[[2]], x2[[2]],x3[[2]], common.legend = T, nrow=3)
ggsave(paste0("./plots/same_seq_sens_", Sys.Date(), ".png"), xplot)

#plot difference between predicted and true disease class
y1 <- class_diff(dat4sens, 4, "A")[[2]]
y2 <- class_diff(dat8sens, 8, "B")[[2]]
y3 <- class_diff(dat12sens, 12, "C")[[2]]

#arrange plots and save
ggarrange(y1, y2,y3, common.legend = T)
yplot <- ggarrange(y1, y2,y3, common.legend = T, nrow=3)
ggsave(paste0("./plots/class_res_sens_", Sys.Date(), ".png"), yplot, unit="in", width=7, height=14, dpi=500)

#plot Kendall Tau distance
z1 <- seq_dist(seq4sens, 4, "A")
z2 <- seq_dist(seq8sens, 8, "B")
z3 <- seq_dist(seq12sens, 12, "C")

#arrange plots and save
ggarrange(z1[[2]], z2[[2]], z3[[2]], common.legend = T, ncol=3)
zplot <- ggarrange(z1[[3]], z2[[3]], z3[[3]], common.legend = T, ncol=1)
ggsave(paste0("./plots/seq_dist_sens_", Sys.Date(), ".png"), zplot)

#plot Kendall tau distance and add table for proportion correct seq below
f1=plot_grid(z1[[5]],proptab(x1[[1]]),nrow=2, rel_heights = c(2.5, 1))
f2=plot_grid(z2[[5]],proptab(x2[[1]]),nrow=2, rel_heights = c(2.5, 1))
f3=plot_grid(z3[[5]],proptab(x3[[1]]),nrow=2, rel_heights = c(2.5, 1))

f=ggarrange(f1,f2,f3, nrow=2, ncol=2)

png(paste0("./plots/seq_res_sens_", Sys.Date(), ".png"), unit="in", width=12, height=10, res=500)
f
dev.off()
