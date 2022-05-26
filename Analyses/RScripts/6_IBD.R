##This script calculates pairwise between wild Q. acerifolia and Z. 
#integrifolia individuals and their physical distance. 
#First we calculate the Fst between wild QUAC and ZAIN populations 
#and then regress it with the distance between populations. 
#We use "cleaned" data files here, which means genetic 
#data files have been cleaned for clones and individuals with 
#too much missing data (25% missing data)

#########################
#        Libraries      #
#########################

library(adegenet)
library(hierfstat)
library(geosphere)

#########################
#   Load Data Files     #
#########################
#set working directory
setwd("../../Data_Files")

#list data files 
#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames/", pattern = "_clean_df.csv")

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#############################
#     Fst Calculations      #
#############################

for(sp in 1:length(scenario_list)){
  
  #read in genepop file as a genind object
  sp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  #read in data frame 
  sp_df <- read.csv(paste0("Data_Frames/",sp_df_list[[sp]]))
  
  #name rows as individuals
  rownames(sp_genind@tab) <- sp_df[,1]
  
  #name populations 
  levels(sp_genind@pop) <- sp_df[,2]
  
  #limit data frame object to wild individuals  
  sp_wild_df <- sp_df[sp_df[,3] == "Wild",]
  
  #limit genind object to wild individuals
  sp_wild_genind <- sp_genind[sp_wild_df[,1] %in% rownames(sp_genind@tab),]

  #run hierfstat on 
#  sp_hierfstat <- genind2hierfstat(sp_wild_genind)
}

#convert to hierfstat format object 
QUAC_hierfstat <- genind2hierfstat(QUAC_wild_gen)

#run pairwise fst code 
QUAC_fst_df <- pairwise.neifst(QUAC_hierfstat)

#calculate geographic distances between mean locations
QUAC_dist <- matrix(nrow = length(QUAC_wildpop_names), ncol = length(QUAC_wildpop_names))

for(first in 1:length(QUAC_wildpop_names)){
  for(second in 1:length(QUAC_wildpop_names)){
    QUAC_dist[first,second] <-  distm(QUAC_coords[first,], QUAC_coords[second,], fun = distGeo)/1000
  }
}

#replacce NAs with zeroes 
QUAC_dist[is.na(QUAC_dist)] <- 0
QUAC_fst_df[is.na(QUAC_fst_df)] <- 0

#create a linear regression
QUAC_fst_dist <- lm(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)])

#visualize the isolation by distance relationship with p-value
pdf("../QUAC_analyses/Results/Clustering/QUAC_Dist_Fst.pdf")
plot(QUAC_fst_df[lower.tri(QUAC_fst_df)]~QUAC_dist[lower.tri(QUAC_dist)], pch = 17, ylim = c(0,0.13), 
     xlim = c(0,200),
     xlab = c("Distance (km)"), ylab = c("Fst"))
abline(QUAC_fst_dist)
legend('bottomleft', legend = c("R2 = -0.12","p-value = 0.865"), bty = 'n')
dev.off()
