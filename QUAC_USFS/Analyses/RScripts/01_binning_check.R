#####################
#     Libraries     #
#####################

library(diveRsity)
library(adegenet)
library(poppr)

#########################
#   Load Data Files     #
##########################
#set working directory
setwd("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/QUAC_USFS/Data_Files")

#read in genepop file as a genind object 
USFS_QUAC_gen <- read.genepop("Genepop_Files/USFS_QUAC_genepop.gen", ncode = 2)

#read in USFS QUAC data file 
USFS_QUAC_df <- read.csv("CSV_Files/USFS_QUAC_Scores_df.csv")
USFS_QUAC_df <- USFS_QUAC_df[1:206,]
#pop name list
pop_list <- unique(USFS_QUAC_df$Pop)

#name populations 
levels(USFS_QUAC_gen@pop) <- pop_list 

###############################################
#     Binning Reorganizing Genind Objects     #
###############################################
##separate by clones and wild individuals
#Scion individuals genind object
QUAC_scions_gen <- repool(seppop(USFS_QUAC_gen)[1:4])
#rename levels
levels(QUAC_scions_gen@pop) <- rep("Scions", rep(length(levels(QUAC_scions_gen@pop))))

#wild individuals wild genind object
QUAC_wild_gen <- repool(seppop(USFS_QUAC_gen)[5:8])
#rename levels 
levels(QUAC_wild_gen@pop) <- rep("Wild", rep(length(levels(QUAC_wild_gen@pop))))

#repool into one genind object
QUAC_scion_wild_gen <- repool(QUAC_scions_gen, QUAC_wild_gen)

#create genepop object 
QUAC_scion_wild_genpop <- genind2genpop(QUAC_scion_wild_gen)


