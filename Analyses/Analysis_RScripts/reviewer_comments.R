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
#QUAC_garden_genind <- repool(seppop(QUAC_woK_genind)[1:17])

#rename pops
#levels(QUAC_garden_genind@pop) <- rep("Garden",17)

#create wild genind object
#QUAC_wild_genind <- repool(seppop(QUAC_woK_genind)[18:21])

#rename pops 
#levels(QUAC_wild_genind@pop) <- rep("Wild",4)

#recombine into garden/wild genind object
#QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

#convert to the wild genpop object
#QUAC_wild_genpop <- genind2genpop(seppop(QUAC_garden_wild_genind)[[2]])

#calculate how alleles are represented ex situ
QUAC_all_rep <- colSums(seppop(QUAC_garden_wild_genind)[[1]]@tab,na.rm=T)

#calculate the allele categories in the wild populations
QUAC_all_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, num_wild_ind, n_drop = 0, glob_only = TRUE)	

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

###################################
#     Trying with Sean's Code     #
###################################
#create garden genind
#num_garden_ind<-sum(table(QUAC_woK_genind@pop)[1:17])

#QUAC_garden_genind <- QUAC_woK_genind[1:num_garden_ind,]

#rename pops
#levels(QUAC_garden_genind@pop) <- rep("Garden",17)

#create wild genind object 
#num_wild_ind <- sum(table(QUAC_woK_genind@pop)[18:21])

#QUAC_wild_genind <- QUAC_woK_genind[(num_garden_ind+1):(num_garden_ind+num_wild_ind),]

#rename pops 
#levels(QUAC_wild_genind@pop) <- rep("Wild",4)

#recombine into garden/wild genind object
#QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

#convert to the wild genpop object
#QUAC_wild_genpop <- genind2genpop(seppop(QUAC_wild_genind)[[2]])

#calculate how alleles are represented ex situ
#QUAC_all_rep <- colSums(QUAC_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
#QUAC_all_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, num_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
#QUAC_all_cat <- QUAC_all_cat[1:5]

##create a list to store the individual numbers 
#list 

#num_rep_list <- list(list(), list(), list(), list(), list())

#for(cat in 1:length(QUAC_all_cat)){
  
#  num_alleles_in_cat <- length(QUAC_all_cat[[cat]])
  
#  for (a in 1:num_alleles_in_cat){
    
#      num_rep_list[[cat]][a] <- sum(QUAC_garden_genind@tab[,QUAC_all_cat[[cat]]][,a] > 0, na.rm=T)

#  }
#}

#create data frame to save results  
#QUAC_rep_df <- matrix(nrow = length(dup_reps),
#                      ncol = length(QUAC_all_cat))

#for(dup in dup_reps){
#  for(cat in 1:length(QUAC_all_cat)){
    
    #create data frame to store results 
#    QUAC_rep_df[dup+1,cat] <- sum(num_rep_list[[cat]]>dup)/length(QUAC_all_cat[[cat]])

    
#  }
#}

#QUAC_rep_df <- signif(QUAC_rep_df*100,3)
#colnames(QUAC_rep_df) <- all_cat_list
#rownames(QUAC_rep_df) <- paste0(c(1:10), " or more copies")

#write.csv(QUAC_rep_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rep_df.csv")

####ZAIN 
#load in ZAIN data file
#ZAIN_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen",
#                            ncode = 3)

#create garden genind
#ZAIN_garden_ind <- sum(table(ZAIN_genind@pop)[1:10])

#ZAIN_garden_genind <- ZAIN_genind[1:ZAIN_garden_ind,]

#create wild genind object 
#ZAIN_wild_ind <- sum(table(ZAIN_genind@pop)[c(11:19, 23:26, 28:32, 34:35)])

#ZAIN_wild_genind <- ZAIN_genind[(ZAIN_garden_ind+1):(ZAIN_garden_ind+ZAIN_wild_ind),]

#convert to the wild genpop object
#ZAIN_wild_genpop <- genind2genpop(ZAIN_wild_genind)

#calculate how alleles are represented ex situ
#ZAIN_all_rep <- colSums(ZAIN_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
#ZAIN_all_cat <- get.allele.cat(ZAIN_wild_genpop, 1, 1, ZAIN_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
#ZAIN_all_cat <- ZAIN_all_cat[1:5]

##create a list to store the individual numbers 
#list 
#ZAIN_num_rep_list <- list(list(), list(), list(), list(), list())

#for(cat in 1:length(ZAIN_all_cat)){
  
#  ZAIN_num_alleles_in_cat <- length(ZAIN_all_cat[[cat]])
  
#  for (a in 1:ZAIN_num_alleles_in_cat){
    
#    ZAIN_num_rep_list[[cat]][a] <- sum(ZAIN_garden_genind@tab[,ZAIN_all_cat[[cat]]][,a] > 0, na.rm=T)
    
#  }
#}

#create data frame to save results  
#ZAIN_rep_df <- matrix(nrow = length(dup_reps),
#                        ncol = length(ZAIN_all_cat))

#for(dup in dup_reps){
#  for(cat in 1:length(ZAIN_all_cat)){
    
    #create data frame to store results 
#    ZAIN_rep_df[dup+1,cat] <- sum(ZAIN_num_rep_list[[cat]]>dup)/length(ZAIN_all_cat[[cat]])
    
    
#  }
#}
#colnames(ZAIN_rep_df) <- all_cat_list
#rownames(ZAIN_rep_df) <- paste0(c(1:10)," or more copies")

#ZAIN_rep_df <- signif(ZAIN_rep_df*100,3)

###################################
#     Trying with Sean's Code     #
###################################
#create garden genind
num_garden_ind<-sum(table(QUAC_woK_genind@pop)[1:17])

QUAC_garden_genind <- QUAC_woK_genind[1:num_garden_ind,]

#rename pops
levels(QUAC_garden_genind@pop) <- rep("Garden",17)

#create wild genind object 
num_wild_ind <- sum(table(QUAC_woK_genind@pop)[18:21])

QUAC_wild_genind <- QUAC_woK_genind[(num_garden_ind+1):(num_garden_ind+num_wild_ind),]

#rename pops 
levels(QUAC_wild_genind@pop) <- rep("Wild",4)

#recombine into garden/wild genind object
#QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

#convert to the wild genpop object
QUAC_wild_genpop <- genind2genpop(QUAC_wild_genind)

#calculate how alleles are represented ex situ
QUAC_all_rep <- colSums(QUAC_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
QUAC_all_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, num_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
QUAC_all_cat <- QUAC_all_cat[1:5]

##create a list to store the individual numbers 
#list 

num_rep_list <- list(list(), list(), list(), list(), list())
num_rep_list_he <- list(list(), list(), list(), list(), list())
num_rep_list_ho <- list(list(), list(), list(), list(), list())

for(cat in 1:length(QUAC_all_cat)){
  
  num_alleles_in_cat <- length(QUAC_all_cat[[cat]])
  
  for (a in 1:num_alleles_in_cat){
    
    num_rep_list[[cat]][a] <- sum(QUAC_garden_genind@tab[,QUAC_all_cat[[cat]]][,a] > 0, na.rm=T)
    num_rep_list_he[[cat]][a] <- sum(QUAC_garden_genind@tab[,QUAC_all_cat[[cat]]][,a] == 1, na.rm=T)
    num_rep_list_ho[[cat]][a] <- sum(QUAC_garden_genind@tab[,QUAC_all_cat[[cat]]][,a] == 2, na.rm=T)
  }
}

#create data frame to save results  
QUAC_rep_df <- matrix(nrow = length(dup_reps),
                      ncol = length(QUAC_all_cat))

QUAC_rep_df_he <- matrix(nrow = length(dup_reps),
                         ncol = length(QUAC_all_cat))
QUAC_rep_df_ho <- matrix(nrow = length(dup_reps),
                         ncol = length(QUAC_all_cat))

for(dup in dup_reps){
  for(cat in 1:length(QUAC_all_cat)){
    
    #create data frame to store results 
    QUAC_rep_df[dup+1,cat] <- sum(num_rep_list[[cat]]>dup)/length(num_rep_list[[cat]])
    QUAC_rep_df_he[dup+1,cat] <- sum(num_rep_list_he[[cat]]>dup)/length(QUAC_all_cat[[cat]])
    QUAC_rep_df_ho[dup+1,cat] <- sum(num_rep_list_ho[[cat]]>dup)/length(QUAC_all_cat[[cat]])
    
  }
}


#write out matrix
QUAC_rep_df <- signif(QUAC_rep_df*100,3)
colnames(QUAC_rep_df) <- all_cat_list
rownames(QUAC_rep_df) <- paste0(c(1:10), " or more copies")

write.csv(QUAC_rep_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rep_df.csv")

#write out het matrix 
QUAC_rep_df_he <- signif(QUAC_rep_df_he*100,3)
colnames(QUAC_rep_df_he) <- all_cat_list
rownames(QUAC_rep_df_he) <- paste0(c(1:10), " or more copies")
write.csv(QUAC_rep_df_he, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rep_he_df.csv")

#write out homo matrix
QUAC_rep_df_ho <- signif(QUAC_rep_df_ho*100,3)
colnames(QUAC_rep_df_ho) <- all_cat_list
rownames(QUAC_rep_df_ho) <- paste0(c(1:10), " or more copies")
write.csv(QUAC_rep_df_ho, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rep_ho_df.csv")


####ZAIN 
#load in ZAIN data file
ZAIN_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen",
                            ncode = 3)

#create garden genind
ZAIN_garden_ind <- sum(table(ZAIN_genind@pop)[1:10])

ZAIN_garden_genind <- ZAIN_genind[1:ZAIN_garden_ind,]

#create wild genind object 
ZAIN_wild_ind <- sum(table(ZAIN_genind@pop)[c(11:19, 23:26, 28:32, 34:35)])

ZAIN_wild_genind <- ZAIN_genind[(ZAIN_garden_ind+1):(ZAIN_garden_ind+ZAIN_wild_ind),]

#convert to the wild genpop object
ZAIN_wild_genpop <- genind2genpop(ZAIN_wild_genind)

#calculate how alleles are represented ex situ
ZAIN_all_rep <- colSums(ZAIN_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
ZAIN_all_cat <- get.allele.cat(ZAIN_wild_genpop, 1, 1, ZAIN_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
ZAIN_all_cat <- ZAIN_all_cat[1:5]

##create a list to store the individual numbers 
#list 
ZAIN_num_rep_list <- list(list(), list(), list(), list(), list())
ZAIN_num_rep_list_he <- list(list(), list(), list(), list(), list())
ZAIN_num_rep_list_ho <- list(list(), list(), list(), list(), list())

for(cat in 1:length(ZAIN_all_cat)){
  
  ZAIN_num_alleles_in_cat <- length(ZAIN_all_cat[[cat]])
  
  for (a in 1:ZAIN_num_alleles_in_cat){
    
    ZAIN_num_rep_list[[cat]][a] <- sum(ZAIN_garden_genind@tab[,ZAIN_all_cat[[cat]]][,a] > 0, na.rm=T)
    ZAIN_num_rep_list_he[[cat]][a] <- sum(ZAIN_garden_genind@tab[,ZAIN_all_cat[[cat]]][,a] == 1, na.rm=T)
    ZAIN_num_rep_list_ho[[cat]][a] <- sum(ZAIN_garden_genind@tab[,ZAIN_all_cat[[cat]]][,a] == 2, na.rm=T)
    
  }
}

#create data frame to save results  
ZAIN_rep_df <- matrix(nrow = length(dup_reps),
                      ncol = length(ZAIN_all_cat))
ZAIN_rep_df_he <- matrix(nrow = length(dup_reps),
                         ncol = length(ZAIN_all_cat))
ZAIN_rep_df_ho <- matrix(nrow = length(dup_reps),
                         ncol = length(ZAIN_all_cat))

for(dup in dup_reps){
  for(cat in 1:length(ZAIN_all_cat)){
    
    #create data frame to store results 
    ZAIN_rep_df[dup+1,cat] <- sum(ZAIN_num_rep_list[[cat]]>dup)/length(ZAIN_all_cat[[cat]])
    
    
  }
}
colnames(ZAIN_rep_df) <- all_cat_list
rownames(ZAIN_rep_df) <- paste0(c(1:10)," or more copies")

ZAIN_rep_df <- signif(ZAIN_rep_df*100,3)

write.csv(ZAIN_rep_df, "../Analyses/Results/Garden_Wild_Comparison/ZAIN_rep_df.csv")

