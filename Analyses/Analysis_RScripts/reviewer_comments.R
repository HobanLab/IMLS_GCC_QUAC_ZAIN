#####################
#     Libraries     #
#####################

library(adegenet)
library(sjmisc)

###########################
#     Load Data Files     #
###########################
setwd("../../Data_Files")

#first try code on JUST QUAC without Kessler 
QUAC_woK_genind <- read.genepop("Adegenet_Files/QUAC_woK_allpop_clean.gen",
                                ncode = 3)

#load in fa sample functions
source("../Analyses/Functions/Fa_sample_funcs.R")

#allele categories list
all_cat_list <-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")

#list of duplicate reps
dup_reps <- c(0:9)

#################################################################
#     Calculate Individual Numbers That Contain Each Allele     #
#################################################################
##convert to garden/wild genind object 
#create garden genind
QUAC_garden_genind <- repool(seppop(QUAC_woK_genind)[1:17])
#rename pops
levels(QUAC_garden_genind@pop) <- rep("Garden",17)

#create wild genind object
QUAC_wild_genind <- repool(seppop(QUAC_woK_genind)[18:21])
#rename pops 
levels(QUAC_wild_genind@pop) <- rep("Wild",4)

#recombine into garden/wild genind object
QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

#convert to the wild genpop object
QUAC_wild_genpop <- genind2genpop(seppop(QUAC_garden_wild_genind)[[2]])

#create documents for comparison 
n_ind_W <- nrow(QUAC_wild_genind@tab);  n_ind_G <- nrow(QUAC_garden_genind@tab)

#calculate how alleles are represented ex situ
QUAC_all_rep <- colSums(seppop(QUAC_garden_wild_genind)[[1]]@tab,na.rm=T)

#calculate the allele categories in the wild populations
QUAC_all_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, n_ind_W, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
QUAC_all_cat <- QUAC_all_cat[1:5]

#create a dataframe of just allele representation counts
QUAC_all_rep_df <- as.data.frame(QUAC_all_rep_df)

###count how many individuals each alleles is represented by 
##global alleles
#create matrix 
QUAC_all_glob_df <- matrix(nrow = length(QUAC_all_cat[[1]]),
                           ncol = 3)
#name rows with allele
rownames(QUAC_all_glob_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[1]]]

#name columns
colnames(QUAC_all_glob_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                "Number_of_Ind")

#loop for global alleles represented  
for(g in 1:length(QUAC_all_cat[[1]])){
  
  ##very common alleles
  #homozygotes
  QUAC_all_glob_df[g,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] == 2, 
                               na.rm = TRUE)
  
  #heterozygotes
  QUAC_all_glob_df[g,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] == 1,
                                na.rm = TRUE)
  
  #total ind
  QUAC_all_glob_df[g,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] > 0,
                                na.rm = TRUE)
}

##very common alleles
#create matrix
QUAC_all_vcom_df <- matrix(nrow = length(QUAC_all_cat[[2]]),
                           ncol = 3)

#name rows with allele
rownames(QUAC_all_vcom_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[2]]]

#name columns
colnames(QUAC_all_vcom_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                "Number_of_Ind")

#loop for very common alleles 
for(vc in 1:length(QUAC_all_cat[[2]])){
  
  ##very common alleles
  #homozygotes
  QUAC_all_vcom_df[vc,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[2]][vc]] == 2, 
                                na.rm = TRUE)
  #heterozygotes
  QUAC_all_vcom_df[vc,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[2]][vc]] == 1,
                                na.rm = TRUE)
  
  #total ind
  QUAC_all_vcom_df[vc,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[2]][vc]] > 0,
                                na.rm = TRUE)
}

##common alleles 
#create matrix
QUAC_all_com_df <- matrix(nrow = length(QUAC_all_cat[[3]]),
                          ncol = 3)
#name rows
rownames(QUAC_all_com_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[3]]]
#name columns
colnames(QUAC_all_com_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                "Number_of_Ind")

#loop for common alleles
for(c in 1:length(QUAC_all_cat[[3]])){
  
  ##common alleles
  #homos
  QUAC_all_com_df[c,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[3]][c]] == 2, 
                              na.rm = TRUE)
  #hets
  QUAC_all_com_df[c,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[3]][c]] == 1, 
                              na.rm = TRUE)
  #all ind
  QUAC_all_com_df[c,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[3]][c]] > 0, 
                              na.rm = TRUE)
  
}

##low frequency alleles 
#create matrix 
QUAC_all_lowfreq_df <- matrix(nrow = length(QUAC_all_cat[[4]]),
                              ncol = 3)
#rownames
rownames(QUAC_all_lowfreq_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[4]]]
#name columns
colnames(QUAC_all_lowfreq_df) <- c("Number_of_Homos", "Number_of_Hets", 
                               "Number_of_Ind")

#loop to calculate the number of individuals with each allele 
for(lf in 1:length(QUAC_all_cat[[4]])){
  
  ##low frequency 
  #homos
  QUAC_all_lowfreq_df[lf,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[4]][lf]] == 2, 
                                   na.rm = TRUE)
  #hets
  QUAC_all_lowfreq_df[lf,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[4]][lf]] == 1, 
                                   na.rm = TRUE)
  #all ind
  QUAC_all_lowfreq_df[lf,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[4]][lf]] > 0, 
                                   na.rm = TRUE)
}

##rare alleles 
#create matrix
QUAC_all_rare_df <- matrix(nrow = length(QUAC_all_cat[[5]]),
                           ncol = 3)
#name rows
rownames(QUAC_all_rare_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[5]]]
#name columns
colnames(QUAC_all_rare_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                   "Number_of_Ind")

#loop to calculate how many alleles each individuals
for(r in 1:length(QUAC_all_cat[[5]])){
  
    ##rare alleles
    #homos
    QUAC_all_rare_df[r,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]]][r] == 2, 
                                     na.rm = TRUE)
    
    #hets
    QUAC_all_rare_df[r,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]][r]] == 1, 
                                 na.rm = TRUE)
  
    #all ind 
    QUAC_all_rare_df[r,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]]][r] > 0, 
                               na.rm = TRUE)
}

#now calculate the numbers of alleles that are in number of individuals 



#create a data frame to store results 
QUAC_ind_all_con_df <- as.data.frame(matrix(nrow = length(dup_reps),
                                         ncol = length(all_cat_list)))

#loop to sum how many individuals represent a copy of each allele on average
for(dup in dup_reps){
  
  #global
  QUAC_ind_all_con_df[dup+1,1] <- signif(sum(QUAC_all_glob_df[,3][QUAC_all_glob_df[,3] > dup])/length(QUAC_all_cat[[1]]),2)
  
  #very common
  QUAC_ind_all_con_df[dup+1,2] <- signif(sum(QUAC_all_vcom_df[,3][QUAC_all_vcom_df[,3] > dup])/length(QUAC_all_cat[[2]]),3)
  
  #common
  QUAC_ind_all_con_df[dup+1,3] <- signif(sum(QUAC_all_com_df[,3][QUAC_all_com_df[,3] > dup])/length(QUAC_all_cat[[3]]),2)
  
  #low frequency
  QUAC_ind_all_con_df[dup+1,4] <- signif(sum(QUAC_all_com_df[,3][QUAC_all_com_df[,3] > dup])/length(QUAC_all_cat[[4]]),2)
  
  
  #for rare alleles
  QUAC_ind_all_con_df[dup+1,5] <- signif(sum(QUAC_all_rare_df[,3][QUAC_all_rare_df[,3] > dup])/length(QUAC_all_cat[[5]]),1)
  
}

#organize df
colnames(QUAC_ind_all_con_df) <- all_cat_list
rownames(QUAC_ind_all_con_df) <- paste0(rownames(QUAC_ind_all_con_df),
                                        " or more copies")
write.csv(QUAC_ind_all_con_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_woK_ind_rep.csv") 
