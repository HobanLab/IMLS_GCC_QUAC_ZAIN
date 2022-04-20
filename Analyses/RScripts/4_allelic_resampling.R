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
setwd("../Analyses/Results/Garden_Wild_Comparison")

resampling_list <- list.files(pattern = "resampling_df")
name_list <- c("QUAC_wK_ndrop0", "QUAC_wK_ndrop2", "QUAC_woK_ndrop0", "QUAC_woK_ndrop2",
               "ZAIN_og_ndrop0", "ZAIN_og_ndrop2", "ZAIN_rebinned_ndrop0", "ZAIN_rebinned_ndrop2")

for(sp in 1:length(resampling_list)){
  
  #load in data files 
  sp_resampling <- read.csv(resampling_list[[sp]])
  
  #export plots for each 
  sp_resampling_plot <- resampling_plot(sp_resampling, name_list[[sp]])
  
}

##Create data frame of min sample size to sample 95% of diversity
#create a data frame to store results 
all_sp_min_samp_95 <- matrix(nrow = length(n_drop), ncol = length(list_allele_cat))

#######minimum sample size loop   
resampling_list <- list.files(pattern = "resampling_df")
name_list <- c("QUAC_wK_ndrop0", "QUAC_wK_ndrop2", "QUAC_woK_ndrop0", "QUAC_woK_ndrop2",
               "ZAIN_og_ndrop0", "ZAIN_og_ndrop2", "ZAIN_rebinned_ndrop0", "ZAIN_rebinned_ndrop2")

#create data frame 
sp_min_sample_95 <- matrix(nrow = length(name_list), ncol = length(list_allele_cat))

for(sp in 1:length(resampling_list)){
  
  sp_resampling_df <- read.csv(resampling_list[[sp]])
  
  #clean up data frame 
  sp_resampling_df <- sp_resampling_df[-1,c(2:10)]
  
  for(all_cat in 1:length(list_allele_cat)){
    
    #line to store in a data frame 
    sp_min_sample_95[sp,all_cat] <- which(sp_resampling_df[,all_cat] >= 95)[1]
    
  }
  
}


##name rows and columns of the matrix 
rownames(sp_min_sample_95) <- name_list
colnames(sp_min_sample_95) <- list_allele_cat

##write out data frame 
write.csv(sp_min_sample_95, "sp_min_samp_95.csv")
