###This script details calculating how many wild samples it would take to 
##Capture different levels of diversity in wild samples. 
#

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)
library(tidyr)

#######################
#     Load files      #
#######################
#set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/Garden_Wild", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames", pattern = "_clean_df.csv")

#list out allele categories
list_sp_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#list of scenarios 
species_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#load in function to calculate allele frequency categories
source("../Analyses/RScripts/Fa_sample_funcs.R")

###################################
#     Allelic Resampling Code     #
###################################
#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(species_list)){
    
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/Garden_Wild/",sp_genind_list[[sp]]), ncode = 3)
  
  #load data frames 
  sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]])) 
  
  #organize genind object
  rownames(sp_genind_temp@tab) <- sp_df_temp[,1]
  levels(sp_genind_temp@pop) <- unique(sp_df_temp[,3])
  
  #run resampling code on all species 
  sp_resampling <- resampling(sp_genind_temp, species_list[[sp]])
  
}

########################################
#     Reporting Resampling Results     #
########################################
#create a list to store frequency summary results
all_mean_list <- list()

#create a data frame with mean  
all_mean_list[[1]] <- as.data.frame(apply(summ_results_tree_ndrop0[,,1:num_reps],c(1,2),mean,na.rm=T)*100)[-1,]
all_mean_list[[2]] <- as.data.frame(apply(summ_results_tree_ndrop2[,,1:num_reps],c(1,2),mean,na.rm=T)*100)[-1,]

##loop to plot resampling results including and removing very rare alleles  
for(ndrop in 1:length(all_mean_list)){
  
  #write PDF with name
  pdf(paste0("../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_all_resampling_ndrop", n_drop[[ndrop]],".pdf"))
  #add points
  plot(all_mean_list[[ndrop]][,1], col = "red", pch = 20, xlab = "Number of Individuals", 
       ylab = "Percent Diversity Capture", xlim = c(0,171), ylim = c(0,100), cex = 1.2,
       main = "Percent Diversity Capture (All Alleles Included)")
  points(all_mean_list[[ndrop]][,2], col = "firebrick", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,3], col = "darkorange3", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,4], col = "coral", pch = 20, cex = 1.2)
  points(all_mean_list[[ndrop]][,5], col = "deeppink4", pch = 20, cex = 1.2)
  
  dev.off()
}

##Create data frame of min sample size to sample 95% of diversity
#create a data frame to store results 
min_samp_95 <- matrix(nrow = length(n_drop), ncol = length(list_allele_cat))

##loop to calculate min sample size
for(ndrop in 1:length(n_drop)){ 
  for(col in 1:length(list_allele_cat)){
    
    min_samp_95[ndrop,col] <- which(all_mean_list[[ndrop]][,col] >= 95)[1]
    
  }
}
##name rows and columns of the matrix 
rownames(min_samp_95) <- c("N Ind for 95% (All Alleles)", 
                           "N Ind for 95% (Rare Alleles Dropped)")
colnames(min_samp_95) <- list_allele_cat

##write out data frame 
write.csv(min_samp_95, "../QUAC_analyses/Results/Wild_Garden_Comparison/QUAC_min_samp_95.csv")
