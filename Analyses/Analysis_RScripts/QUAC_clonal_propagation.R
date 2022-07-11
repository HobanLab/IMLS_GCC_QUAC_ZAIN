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

#load in selected ind data frame 
QUAC_clo_prop_ind <- read.csv("G:/Shared drives/Emily_Schumacher/Project_Folders/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_prop_ind.csv")

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
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_allresampling_garden_wild.pdf", width = 8, height = 8)
plot(garden_all_mean[,1], col = "seagreen1", pch = 15, xlim = c(0,300), xlab = "Number of individuals",
     ylab = "Percent Allele Capture", ylim = c(0,100))
points(wild_all_mean[,1], col = "darkgreen", pch = 15)
#abline(h = 95, col = "darkslategray4", lty=2, lwd = 3)
legend('bottomright', legend = c("Garden", "Wild"), col = c("seagreen1", "darkgreen"),
       pch = 15)
#legend('bottomright', legend = c("Garden", "Wild", "95% diversity capture"), lwd = 2, cex = 1.2, 
#       pch = c(15, 15, NA), col = c("seagreen1", "darkgreen", "darkslategray4"), 
#       lty = c(NA, NA, 1), bty = 'n')
dev.off()

#plot zoom
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_allresampling_garden_wild_zoom.pdf", width = 8, height = 8)
plot(garden_all_mean[,1], col = "seagreen1", pch = 15, xlim = c(0,300), xlab = "Number of individuals",
     ylab = "Percent Allele Capture", ylim = c(60,100))
points(wild_all_mean[,1], col = "darkgreen", pch = 15)
abline(h = 95, col = "darkslategray4", lty=2, lwd = 3)

legend('bottom', legend = c("Garden", "Wild", "95% diversity capture"), lwd = 2, cex = 1.2, 
       pch = c(15, 15, NA), col = c("seagreen1", "darkgreen", "darkslategray4"), 
       lty = c(NA, NA, 1), bty = 'n')
dev.off()

###################################################################
#     Genetic Diversity Analyses on Clonal Prop Selected Ind      #
###################################################################
#load QUAC genind object
QUAC_genind <- read.genepop("Adegenet_Files/Garden_Wild/QUAC_woK_garden_wild_clean.gen", ncode = 3)
#load data frame
QUAC_df <- read.csv("Data_Frames/QUAC_woK_allpop_clean_df.csv")
rownames(QUAC_genind@tab) <- QUAC_df$Ind
levels(QUAC_genind@pop) <- unique(QUAC_df$Garden_Wild)

QUAC_garden_genind <- seppop(QUAC_genind)[[1]]

QUAC_clo_ind_genind <- QUAC_garden_genind[rownames(QUAC_garden_genind@tab) %in% 
                                     QUAC_clo_prop_ind$Hoban_Database_ID,]

levels(QUAC_clo_ind_genind@pop) <- "clo_selected_ind"

#create a random sample of 79 individuals 
QUAC_random_ind <- sample(rownames(seppop(QUAC_genind)[[1]]@tab), 77)

#limit genind by the random individuals 
QUAC_random_ind_genind <- QUAC_genind[rownames(QUAC_genind@tab) %in% 
                                        QUAC_random_ind,]

QUAC_random_ind_wild_genind <- repool(QUAC_random_ind_genind, seppop(QUAC_genind)[[2]])

levels(QUAC_random_ind_genind@pop) <- "random_ind"

#repool both genind objects 
QUAC_clo_random_genind <- repool(QUAC_clo_ind_genind, QUAC_random_ind_genind)

##genetic diversity analyses
#allelic richness 
QUAC_clo_random_gendiv_df <- cbind(unlist(summary(QUAC_clo_random_genind)[2]), 
                                   unlist(summary(QUAC_clo_random_genind)[4]),
                                   colMeans(allelic.richness(QUAC_clo_random_genind)$Ar) ,
                                   poppr(QUAC_clo_random_genind)[1:2,10])

#name rows and columns 
rownames(QUAC_clo_random_gendiv_df) <- c("Prop_Ind", "Random_Ind")
colnames(QUAC_clo_random_gendiv_df) <- c("N", "Na", "All_Rich", "HExp")
write.csv(QUAC_clo_random_gendiv_df, "../Analyses/Results/Sum_Stats/Prop_Ind_Gendiv.csv")

#########################################################
#     Summary of Selected Clonal Prop Ind Geography     #
#########################################################
#summarize individual table 
QUAC_str_pops <- c("Magazine1","Magazine2","Magazine3",
                   "Sugar_Loaf", "Porter","Pryor")
QUAC_clo_pop_df <- matrix(nrow = length(QUAC_clo_prop_ind[,1]), ncol = length(QUAC_str_pops))

#create a loop to check for each population 
for(ind in 1:length(QUAC_clo_prop_ind[,1])){
  
  QUAC_clo_pop_df[ind,1] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Magazine1")
  QUAC_clo_pop_df[ind,2] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Magazine2")
  QUAC_clo_pop_df[ind,3] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Magazine3")
  QUAC_clo_pop_df[ind,4] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Sugar_Loaf")
  QUAC_clo_pop_df[ind,5] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Porter")
  QUAC_clo_pop_df[ind,6] <- str_contains(QUAC_clo_prop_ind$Str_Pop[ind], "Pryor")
  
}
colnames(QUAC_clo_pop_df) <- QUAC_str_pops

#summarize this data frame 
QUAC_clo_prop_pop_sum_df <- matrix(nrow = 1, ncol = length(QUAC_str_pops))

for(pop in 1:length(QUAC_str_pops)){
  
  QUAC_clo_prop_pop_sum_df[,pop] <- length(QUAC_clo_pop_df[QUAC_clo_pop_df[,pop] == TRUE,][,1])
  
}
colnames(QUAC_clo_prop_pop_sum_df) <- QUAC_str_pops

write.csv(QUAC_clo_prop_pop_sum_df, "../Analyses/Results/Sum_Stats/QUAC_clo_prop_pop_sum_df.csv")

##count how many individuals come from the US
#list truths
QUAC_US_ind <- sum(sapply(QUAC_clo_prop_ind$Location, str_contains, "USA"))


