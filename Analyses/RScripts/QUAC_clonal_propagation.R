###This script was created to run analyses for providing recommendations to 
##Mark Coggeshall and Ryan Russel for potentially clonally propagating 
#Q. acerifolia individuals from ex situ populations. Therefore we compare 
#diversity levels in ex situ and in situ populations as well as levels of 
#relatedness

#####################
#     Libraries     #
#####################

library(adegenet)
library(Demerelate)

#######################
#     Load Files      #
#######################
#set working directory 
setwd("../../Data_Files")

#load QUAC genind object
QUAC_genind <- read.genepop("Adegenet_Files/Garden_Wild/QUAC_woK_garden_wild_clean.gen", ncode = 3)
#load data frame
QUAC_df <- read.csv("Data_Frames/QUAC_woK_allpop_clean_df.csv")
rownames(QUAC_genind@tab) <- QUAC_df$Ind
levels(QUAC_genind@pop) <- unique(QUAC_df$Garden_Wild)

#load allele categorization code 
source("../Analyses/RScripts/Fa_sample_funcs.R")

###########################
#     Resampling code     #
###########################
#set number of replicates 
num_reps <- 1000

##run resampling code 
#botanic garden individuals 
sp_genind <- seppop(QUAC_genind)[[1]]
n_total_indivs <- length(sp_genind@tab[,1])
n_ind_p_pop <- table(sp_genind@pop)
#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
#calculate allele category
allele_cat <- get.allele.cat(sp_genind, region_makeup=NULL, 2, n_ind_p_pop,n_drop = 0, glob_only=T)

#create summary results for allelic capture 
garden_results_tree <- array(dim = c((nrow(sp_genind@tab)-1), length(list_allele_cat), num_reps)) 

#create a summary table 
garden_results_df <- array(dim = c((nrow(sp_genind@tab)-1), length(list_allele_cat), num_reps)) 

#Repeat the resampling many times
for (nrep in 1:num_reps) {
  
  #create empty matrix to store sampling code 
  alleles_samp <- matrix(nrow=nrow(sp_genind@tab)-1,ncol=length(list_allele_cat))
  
  #This loop will sample trees from t = 2 to the total number of trees
  for (t in 2:(nrow(sp_genind@tab)-1)){
    
    #create a sample of trees of length t, by using 'sample()' which randomly samples rows
    alleles_samp <- colSums(sp_genind@tab[sample(1:nrow(sp_genind@tab), t),],na.rm=T)
    
    #Then simply compare that sample to your wild population with allele_cat
    for (cat in 1:length(allele_cat)) garden_results_tree[t,cat,nrep] <- sum(alleles_samp[allele_cat[[cat]]]>0, na.rm=T)
    
    #Divide by the number of alleles
    garden_results_df[,,nrep] <- t(t(garden_results_tree[,,nrep])/garden_results_tree[length(garden_results_tree[,1,1]),,nrep])
    
  }
}

#mean across reps using apply
garden_all_mean <- apply(garden_results_df[,,1:num_reps],c(1,2),mean,na.rm=T)*100

#run on wild individuals 
sp_genind <- seppop(QUAC_genind)[[2]]
n_total_indivs <- length(sp_genind@tab[,1])
n_ind_p_pop <- table(sp_genind@pop)
#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
#calculate allele category
allele_cat <- get.allele.cat(sp_genind, region_makeup=NULL, 2, n_ind_p_pop,n_drop = 0, glob_only=T)

#create summary results for allelic capture 
wild_results_tree <- array(dim = c((nrow(sp_genind@tab)-1), length(list_allele_cat), num_reps)) 

#create a summary table 
wild_results_df <- array(dim = c((nrow(sp_genind@tab)-1), length(list_allele_cat), num_reps)) 

#Repeat the resampling many times
for (nrep in 1:num_reps) {
  
  #create empty matrix to store sampling code 
  alleles_samp <- matrix(nrow=nrow(sp_genind@tab)-1,ncol=length(list_allele_cat))
  
  #This loop will sample trees from t = 2 to the total number of trees
  for (t in 2:(nrow(sp_genind@tab)-1)){
    
    #create a sample of trees of length t, by using 'sample()' which randomly samples rows
    alleles_samp <- colSums(sp_genind@tab[sample(1:nrow(sp_genind@tab), t),],na.rm=T)
    
    #Then simply compare that sample to your wild population with allele_cat
    for (cat in 1:length(allele_cat)) wild_results_tree[t,cat,nrep] <- sum(alleles_samp[allele_cat[[cat]]]>0, na.rm=T)
    
    #Divide by the number of alleles
    wild_results_df[,,nrep] <- t(t(wild_results_tree[,,nrep])/wild_results_tree[length(wild_results_tree[,1,1]),,nrep])
    
  }
}

#mean across all wild individuals 
wild_all_mean <- apply(wild_results_df[,,1:num_reps],c(1,2),mean,na.rm=T)*100

#store # of individuals to objects
garden_min_95 <- which(garden_all_mean[,1] >= 95)[1]
#line to store in a data frame 
wild_min_95 <- which(wild_all_mean[,1] >= 95)[1]

##plot both lines 
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_allresampling_garden_wild.pdf", width = 8, height = 10)
plot(garden_all_mean[,1], col = "seagreen1", pch = 15, xlim = c(0,300), xlab = "Number of individuals",
     ylab = "Percent Allele Capture", ylim = c(0,100))
points(wild_all_mean[,1], col = "darkgreen", pch = 15)
legend('bottomright', legend = c("Garden", "Wild"), pch = 15, col = c("seagreen1", "darkgreen"), bty = 'n')
dev.off()
