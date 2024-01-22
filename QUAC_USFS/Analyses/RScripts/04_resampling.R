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

#load genind object 
QUAC_gen <- read.genepop("Genepop_Files/USFS_QUAC_gen_clean.gen", ncode = 2)

#load csv file as a data frame 
QUAC_df <- read.csv("CSV_Files/USFS_QUAC_rebinned_df.csv")

#remove individual with too much missing data
QUAC_df <- QUAC_df[QUAC_df$ID %in% rownames(QUAC_gen@tab),]

#create pop list name
pop_list <- unique(QUAC_df$Pop)

#reorganize genind object 
levels(QUAC_gen@pop) <- pop_list

#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com",
                   "glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#load in function to calculate allele frequency categories
source("../../Analyses/Functions/Fa_sample_funcs.R")
#include function for take the maximum value of a column
colMax <- function(data) sapply(data, max, na.rm = TRUE)

###########################
#     Run Resampling      #
###########################
#set rep number
num_reps <- 1000

#create documents for allelic categorization code 
n_total_indivs <- length(QUAC_wild_gen@tab[,1])
n_ind_p_pop <- table(QUAC_wild_gen@pop)

##repool
QUAC_scion_gen <- repool(seppop(QUAC_gen)[1:4])
levels(QUAC_scion_gen@pop) <- rep("Scion", 4)
QUAC_wild_gen <- repool(seppop(QUAC_gen)[5:8])
levels(QUAC_wild_gen@pop) <- rep("Wild", 4)

QUAC_scion_wild_gen <- repool(QUAC_scion_gen,
                              QUAC_wild_gen)

#convert to genpop 
QUAC_wild_genpop <- genind2genpop(seppop(QUAC_scion_wild_gen)[2]$Wild)

#calculate allelic categorization 
allele_cat <- get.allele.cat(QUAC_wild_genpop, region_makeup=NULL, 2, n_ind_p_pop,
                             n_drop = 0, glob_only=T)	

for(g in 1:length(levels(QUAC_wild_gen@pop))){
  
  #organize the genind objects 
  QUAC_scion_gen <- seppop(QUAC_gen)[[g]]
  QUAC_wild_gen <- seppop(QUAC_gen)[[g+4]]
 
  #This loop will sample trees from t = 2 to the total number of trees
  for (t in 2:(nrow(QUAC_wild_gen@tab)-1)){
    
    #sampling from alleles from wild to generate a matrix 
    #including scions 
    ex_situ_coll <- rbind(tab(QUAC_wild_gen)[sample(1:nrow(tab(QUAC_wild_gen)),
                                                            2),], tab(QUAC_scion_gen))
    #summary of alleles sampled 
    alleles_samp <- colSums(ex_situ_coll)
  
    #Repeat the resampling many times
    for(nrep in 1:num_reps){
    
      #Then simply compare that sample to your wild population with allele_cat
      for(cat in 1:length(allele_cat)) sum_results_tree[t,cat,nrep] <- sum(alleles_samp[allele_cat[[cat]]]>0, na.rm=T)

    #Divide by the number of alleles
    sum_results_df[,,nrep] <- t(t(sum_results_tree[,,nrep])/sum_results_tree[length(sum_results_tree[,1,1]),,nrep])
    
    #mean across reps using apply
    all_mean <- apply(sum_results_df[,,1:num_reps],c(1,2),mean,na.rm=T)*100

    }
    write.csv(all_mean, paste0("../Analyses/Results/Resampling/", pop_list[[g+4]], "_resampling_df.csv"))
  }
}


