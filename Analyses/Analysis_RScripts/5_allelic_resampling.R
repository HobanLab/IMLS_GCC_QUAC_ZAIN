###This script calculates how many wild samples it would take to 
##capture different levels of diversity in wild samples of both ZAIN and QUAC. 
#We do this by resampling wild individuals and determining how many
#alleles we capture. We then determine the number of individuals  
#we would need to sample to capture 95% of diversity, a common target of 
#sampling analyses. We also do this by allele frequency category (rare - global)
#and use code to plot these graphs to visualize the improvement in sampling
#different numbers of individuals. 
#In this script we use genind objects that are referred to as "clean,"
#which means clones and individuals with too much missing data (>25%) 
#have been removed.

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)

#######################
#     Load files      #
#######################
#set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files", pattern = "_clean.gen")

#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#list of scenarios 
species_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#pop list 
pop_list <- list(c(1:17), c(1:17), c(1:10), c(1:10),
                 c(18:22), c(18:21), c(11:35), c(11:35))

#load in function to calculate allele frequency categories
source("../Analyses/Functions/Fa_sample_funcs.R")
source("../Analyses/Functions/resampling.R")

###################################
#     Allelic Resampling Code     #
###################################
#create a list to store arrays 
sp_resampling_list <- list()
all_mean_list <- list()
ndrop_list <- c(0,2)

#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(species_list)){
    
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  #separate wild individuals into a separate genind object 
  sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp+4]]])
  #rename to wild only 
  levels(sp_wild_genind@pop) <- rep("Wild", length(pop_list[[sp+4]]))
  
  #run resampling code on all species 
  for(n in 1:length(ndrop_list)){
   
    #set rep number
    num_reps <- 1000
    #include function for take the maximum value of a column
    colMax <- function(data) sapply(data, max, na.rm = TRUE)
  
    #create documents for allelic categorization code 
    n_total_indivs <- length(sp_wild_genind@tab[,1])
    n_ind_p_pop <- table(sp_wild_genind@pop)
    
    #calculate allele category
    allele_cat <- get.allele.cat(sp_wild_genind, region_makeup=NULL, 2, n_ind_p_pop,n_drop = ndrop_list[[n]], glob_only=T)
  
    #create summary results for allelic capture 
    summ_results_tree <- array(dim = c((nrow(sp_wild_genind@tab)-1), length(list_allele_cat), num_reps)) 
  
    #create a summary table 
    sum_results_df <- array(dim = c((nrow(sp_wild_genind@tab)-1), length(list_allele_cat), num_reps)) 
  
      #Repeat the resampling many times
      for(nrep in 1:num_reps) {
    
        #create empty matrix to store sampling code 
        alleles_samp <- matrix(nrow=nrow(sp_wild_genind@tab)-1,ncol=length(list_allele_cat))
    
          #This loop will sample trees from t = 2 to the total number of trees
          for (t in 2:(nrow(sp_wild_genind@tab)-1)){
      
            #create a sample of trees of length t, by using 'sample()' which randomly samples rows
            alleles_samp <- colSums(sp_wild_genind@tab[sample(1:nrow(sp_wild_genind@tab), t),],na.rm=T)
          
      
      #Then simply compare that sample to your wild population with allele_cat
      for (cat in 1:length(allele_cat)) summ_results_tree[t,cat,nrep] <- sum(alleles_samp[allele_cat[[cat]]]>0, na.rm=T)
      
      #Divide by the number of alleles
      sum_results_df[,,nrep] <- t(t(summ_results_tree[,,nrep])/summ_results_tree[length(summ_results_tree[,1,1]),,nrep])

       }
    #mean across reps using apply
    all_mean <- apply(sum_results_df[,,1:num_reps],c(1,2),mean,na.rm=T)*100
    
    write.csv(all_mean, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_resampling_df", ndrop_list[[n]], ".csv"))
    }
  }
}

########################################
#     Reporting Resampling Results     #
########################################
#set working directory to output results
setwd("../Analyses/Results/Garden_Wild_Comparison")

#collecting all of the resampling data frames
resampling_list <- list.files(pattern = "resampling_df")
#all names for results 
name_list <- c("QUAC_wK_ndrop0", "QUAC_wK_ndrop2", "QUAC_woK_ndrop0", "QUAC_woK_ndrop2",
               "ZAIN_og_ndrop0", "ZAIN_og_ndrop2", "ZAIN_rebinned_ndrop0", "ZAIN_rebinned_ndrop2")

#loop to create the resampling plots 
for(sp in 1:length(resampling_list)){
  
  #load in data files 
  sp_resampling <- read.csv(resampling_list[[sp]])
  
  #export plots for each 
  sp_resampling_plot <- resampling_plot(sp_resampling, name_list[[sp]])
  
}

##Create data frame of min sample size to sample 95% of diversity
#create a data frame to store results 
all_sp_min_samp_95 <- matrix(nrow = length(ndrop_list), ncol = length(list_allele_cat))

##minimum sample size to represent 95% of diversity ex situ   
#create 95% 
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

#name rows and columns of the matrix 
rownames(sp_min_sample_95) <- name_list
colnames(sp_min_sample_95) <- list_allele_cat

#write out data frame 
write.csv(sp_min_sample_95, "sp_min_samp_95.csv")

sessionInfo()
