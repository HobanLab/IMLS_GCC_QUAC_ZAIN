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

#################################################################
#     Calculate Individual Numbers That Contain Each Allele     #
#################################################################

#rename loci in each category - global
for(m in 1:length(QUAC_all_cat[[1]])){
 
  QUAC_all_cat[[1]][m] <- loci_names[[m]] 
  
}

#rename loci in each category - very common
for(m in 1:length(QUAC_all_cat[[2]])){
  
  QUAC_all_cat[[2]][m] <- loci_names[[m]] 
  
}

#rename loci in each category - common
for(m in 1:length(QUAC_all_cat[[3]])){
  
  QUAC_all_cat[[3]][m] <- loci_names[[m]] 
  
}

#rename loci in each category - low frequency
for(m in 1:length(QUAC_all_cat[[4]])){
  
  QUAC_all_cat[[4]][m] <- loci_names[[m]] 
  
}

#rename loci in each category - low frequency
for(m in 1:length(QUAC_all_cat[[5]])){
  
  QUAC_all_cat[[5]][m] <- loci_names[[m]] 
  
}

#allele categories
all_cat_list <-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")

loci_names <- unique(gsub("\\..*", "", colnames(QUAC_garden_genind@tab)))


QUAC_all_rep_df <- as.data.frame(seppop(QUAC_garden_wild_genind)$Garden@tab)

#first, create a data frame of all of the alleles 
#QUAC_all_rep_df %>%
#    select(matches(loci_names[[1]]))

#remove regional alleles 
QUAC_all_cat <- QUAC_all_cat[1:5]

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

##very common alleles
#create matrix
QUAC_all_vcom_df <- matrix(nrow = length(QUAC_all_cat[[2]]),
                           ncol = 3)

#name rows with allele
rownames(QUAC_all_vcom_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[2]]]

#name columns
colnames(QUAC_all_vcom_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                "Number_of_Ind")

##common alleles 
#create matrix
QUAC_all_com_df <- matrix(nrow = length(QUAC_all_cat[[3]]),
                          ncol = 3)
#name rows
rownames(QUAC_all_com_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[3]]]
#name columns
colnames(QUAC_all_com_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                "Number_of_Ind")

##low frequency alleles 
#create matrix 
QUAC_all_lowfreq_df <- matrix(nrow = length(QUAC_all_cat[[4]]),
                              ncol = 3)
#rownames
rownames(QUAC_all_lowfreq_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[4]]]
#name columns
colnames(QUAC_all_lowfreq_df) <- c("Number_of_Homos", "Number_of_Hets", 
                               "Number_of_Ind")

##rare alleles 
#create matrix
QUAC_all_rare_df <- matrix(nrow = length(QUAC_all_cat[[5]]),
                           ncol = 3)
#name rows
rownames(QUAC_all_rare_df) <- colnames(QUAC_all_rep_df)[QUAC_all_cat[[5]]]
#name columns
colnames(QUAC_all_rare_df) <- c("Number_of_Homos", "Number_of_Hets", 
                                   "Number_of_Ind")


QUAC_glob_copy_df <- matrix(nrow = length())


QUAC_rep_df <- matrix(nrow = 3, ncol = 5)

rownames(QUAC_rep_df) <- c("Individuals Ex Situ Representing Allele",
                           "% Homozygotes", "% Heterozygotes")
colnames(QUAC_rep_df) <- c("Global", "Very Common", "Common",
                           "Low Frequency", "Rare")

QUAC_all_global_count_df <- matrix(nrow = length(QUAC_all_cat[[1]]),
                                   ncol = 1)


  
  ##global alleles
  #homozygotes
  #QUAC_all_glob_df[g,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] == 2, 
  #                              na.rm = TRUE)
  
  #heterozygotes
  #QUAC_all_glob_df[g,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] == 1,
  #                             na.rm = TRUE)
  
  #total ind
  #QUAC_all_glob_df[g,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[1]][g]] > 0,
  #                            na.rm = TRUE)

QUAC_all_global_count_df <- matrix(nrow = length(QUAC_all_cat[[1]]),
                                   ncol = length(dup_reps))

#loop for global alleles - counting heterozygotes and homozygotes
for(g in 1:length(QUAC_all_cat[[1]])){
  
  for(d in dup_reps){
  
  QUAC_all_global_count_df[g,d+1] <- mean(length(which(QUAC_all_rep_df[,QUAC_all_cat[[1]]][g] > d)))

  }
}

QUAC_rep_df[1,1] <- round(sum(QUAC_all_glob_df[,3], na.rm = TRUE)/length(QUAC_all_cat[[1]]),0)
QUAC_rep_df[2,1] <- QUAC_rep_df[1,1]

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

QUAC_all_rare_df <- matrix(nrow = length(QUAC_all_cat[[5]]),
                           ncol = length(dup_reps))

QUAC_all_rare_count_df <- matrix(nrow = length(QUAC_all_cat[[5]]),
                                 ncol = 1)

for(r in 1:length(QUAC_all_cat[[5]])){
  
    ##rare alleles
    #homos
    #QUAC_all_rare_df[r,1] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]]][r] == 2, 
    #                                 na.rm = TRUE)
    #hets
    #QUAC_all_rare_df[r,2] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]][r]] == 1, 
    #                             na.rm = TRUE)
  
    #all ind 
    #QUAC_all_rare_df[r,3] <- sum(QUAC_all_rep_df[,QUAC_all_cat[[5]]][r] > 0, 
    #                           na.rm = TRUE)
    
  QUAC_all_rare_count_df[r,1] <- length(which(QUAC_all_rep_df[,QUAC_all_cat[[5]]][r] > 0))
  
}

##create data frames to store
#store number of individuals containing an allele 
QUAC_num_ind_rep_df <- matrix(nrow = length(dup_reps),
                              ncol = length(QUAC_all_cat))
colnames(QUAC_num_ind_rep_df) <- all_cat_list[1:5]


#store percent of individuals (out of total) containing alleles 
QUAC_all_ind_rep_per_df <- matrix(nrow = length(dup_reps),
                              ncol = length(QUAC_all_cat))
colnames(QUAC_all_ind_rep_per_df) <- all_cat_list[1:5]

#loop for 
for(d in dup_reps){
  
  ##store number of individuals that the alleles are represented in   
  #global
 
  sum(QUAC_all_rare_df[,1][QUAC_all_rare_df[,1] > 0])/277
  
}

##divide by the number of individuals to get the percent represented in 
#ex situ collections 
QUAC_all_ind_rep_per_df <- signif((QUAC_num_ind_rep_df/277)*100,3)

#export df matrix 
#QUAC_all_ind_rep_per_clean_df <- paste0(round(QUAC_all_ind_rep_per_df *100,1) , %)

QUAC_all_ind_rep_per_clean_df <- matrix(nrow = length(dup_reps),
                                        ncol = length(QUAC_all_cat))

#loop to clean up data frame and export 
for(d in dup_reps){
  for(c in 1:length(QUAC_all_cat)){
    
  QUAC_all_ind_rep_per_clean_df[d+1,c] <- paste0(QUAC_all_ind_rep_per_df[d+1,c], 
                                          "%", " (", QUAC_num_ind_rep_df[d+1,c],
                                          ")")
  }
}

QUAC_glob_df <- matrix(nrow = length(QUAC_all_cat[[1]]),
                       ncol = 10)

for(a in dup_reps){
  
  QUAC_glob_df[,a+1] <- QUAC_all_glob_df[,3] > a
  
}


 