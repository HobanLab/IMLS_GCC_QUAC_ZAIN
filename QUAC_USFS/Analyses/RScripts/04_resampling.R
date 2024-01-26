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
QUAC_gen <- read.genepop("USFS_QUAC_gen_clean.gen", ncode = 2)

#load csv file as a data frame - this is ONLY used to get population names
QUAC_df <- read.csv("USFS_QUAC_rebinned_df.csv")

#remove individual with too much missing data <-- TO DO: Note this does not remove them from QUAC_gen, which is where they need to be removed.  FIX
QUAC_df <- QUAC_df[QUAC_df$ID %in% rownames(QUAC_gen@tab),]

#create pop list name
pop_list <- unique(QUAC_df$Pop)

#organize genind object by giving the populations the correct names
levels(QUAC_gen@pop) <- pop_list
wild_pops <- pop_list[5:8]

#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com",
                   "glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#load in function to calculate allele frequency categories
source("../../Analyses/Functions/Fa_sample_funcs.R")
#create function to take the maximum value of a column
colMax <- function(data) sapply(data, max, na.rm = TRUE)

###########################
#     Run Resampling      #
###########################
#set rep number
num_reps <- 1000

#create documents for allelic categorization code 

##repool... this will create a WILD and a SCION (EX SITU) set of populations
QUAC_scion_gen <- repool(seppop(QUAC_gen)[1:4])
levels(QUAC_scion_gen@pop) <- rep("Scion", 4)
QUAC_wild_gen <- repool(seppop(QUAC_gen)[5:8])
levels(QUAC_wild_gen@pop) <- rep("Wild", 4)

#Get number of individuals, note this has to come after the lines above
n_total_indivs <- length(QUAC_wild_gen@tab[,1])
n_ind_p_pop <- table(QUAC_wild_gen@pop)

#Put them all back together
QUAC_scion_wild_gen <- repool(QUAC_scion_gen,
                              QUAC_wild_gen)

#convert to genpop 
QUAC_wild_genpop <- genind2genpop(seppop(QUAC_scion_wild_gen)[2]$Wild)

#calculate allelic categorization - TO DO I think this should be removed because I think you want allele_cat calculated for each population (within the g loop below).  Commented out for now
##allele_cat <- get.allele.cat(QUAC_wild_genpop, region_makeup=NULL, 2, n_ind_p_pop,
 ##                            n_drop = 0, glob_only=T)	

for(g in 1:length(wild_pops)){
  
  #organize the genind objects 
  QUAC_scion_gen <- seppop(QUAC_gen)[[g]]
  QUAC_wild_gen <- seppop(QUAC_gen)[[g+4]]
  
  n_indiv<-nrow(QUAC_wild_gen@tab)

#I think its correct to have allele_cat calculated within the g loop.  This is because each loop you only want to consider alleles in that wild population... it would be impossible to collect alleles found in other populations but not the focal population, so you should calculate allele_cat within each population.  It produces some wonky results if you have it outside the g loop (including overestimating sample size needed)

allele_cat <- get.allele.cat(genind2genpop(QUAC_wild_gen), region_makeup=NULL, 2, n_indiv,
                             n_drop = 0, glob_only=T)	
                             
  #create data frames to store results in - these are the results of the number of alleles conserved
  sum_results_tree <- array(dim = c((nrow(QUAC_wild_gen@tab)-1), length(list_allele_cat), num_reps)) 
  
  #create a summary table - this will be the proportions of alleles conserved
  sum_results_df <- array(dim = c((nrow(QUAC_wild_gen@tab)-1), length(list_allele_cat), num_reps)) 
  
  #create matrix 
  ##NOTE- not sure why this is here... all_samp is never used
  all_samp <- matrix(nrow=nrow(QUAC_wild_gen@tab)-1,ncol=length(list_allele_cat))
  
 #Loop over nreps.  Notes the nrep loop should be outside (placed above/ before) the t loop, because you want each sampling replicated many times
for(nrep in 1:num_reps){

  #This loop will sample trees from t = 2 to the total number of trees
     for (t in 2:(nrow(QUAC_wild_gen@tab)-1)){
    
    #sampling from alleles from wild to generate a matrix 
    #including scions 
    ex_situ_coll <- rbind(tab(QUAC_scion_gen), 
                          tab(QUAC_wild_gen)[sample(1:nrow(tab(QUAC_wild_gen)),t),])
    #summary of alleles sampled 
   #Note you need to have an na.rm=T.  Otherwise, anytime an individual has an NA, the whole colSum gets an NA, and since we later compare colSums to 0, it looks like the allele is not conserved- this gives the surprising and incorrect result of proportion of alleles conserved decreases with sample size! 
    alleles_samp <- colSums(ex_situ_coll,na.rm=T)
  
    #Repeat the resampling many times
     
      #Then simply compare that sample to your wild population with allele_cat
      for(cat in 1:length(allele_cat)) sum_results_tree[t,cat,nrep] <- sum(alleles_samp[allele_cat[[cat]]]>0, na.rm=T)

    }
      #Divide by the number of alleles-  needs to be outside the t loop.  Its dividing all the lines ("number of alleles sampled") by the last line, in each slice (the total number of alleles).  Its dividing all the lines in that matrix by the final line.  So it can't happen until that matrix is populated, e.g. after all the t's have been run e.g. outside the t loop.  

     sum_results_df[,,nrep] <- t(t(sum_results_tree[,,nrep])/sum_results_tree[length(sum_results_tree[,1,1]),,nrep])

  }

#mean across reps using apply
#note, all_mean needs to be outside the reps and the t loops.  Its calculating the means across all reps.  Its taking the t x allele cat matrix that is in 1000 slices and taking the mean across all those slices. 

    all_mean <- apply(sum_results_df[,,1:num_reps],c(1,2),mean,na.rm=T)*100
   
    write.csv(all_mean, paste0(pop_list[[g+4]], "_resampling_df.csv"))
  #print out the recommendation of minimum sample size, based on all alleles or just those present at 0.05 locally.
    print(paste(pop_list[[g+4]], which(all_mean[,3]>95)[1],which(all_mean[,1]>95)[1]))
  
}



#####################
#     Visualize     #
#####################

#read in all resampling data files 

resample_list <- list.files(path = "../Analyses/Results/Resampling")

for(n in 1:length(resample_list)){
  
  resample_df <- read.csv(paste0("../Analyses/Results/Resampling/",resample_list[[1]]))
  
  pdf(paste0("../Analyses/Results/Resampling/QUAC", ))
  plot(resample_df[,2], col = "red", pch = 20, xlab = "Number of Individuals", 
       ylab = "Percent Diversity Capture", xlim = c(0,length(rownames(resample_df))), 
       ylim = c(0,100), cex = 1.2,
       main = "Percent Diversity Capture (All Alleles Included)")
  points(resample_df[,4], col = "darkorange3", pch = 20, cex = 1.2)
  points(resample_df[,5], col = "coral", pch = 20, cex = 1.2)
  points(resample_df[,6], col = "deeppink4", pch = 20, cex = 1.2)
  
}

