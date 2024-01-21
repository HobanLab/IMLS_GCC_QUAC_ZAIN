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
setwd("../../Data_Files")

#read in genepop file as a genind object 
USFS_QUAC_gen <- read.genepop("Genepop_Files/USFS_QUAC_rebinned_genepop.gen", 
                              ncode = 2)

#read in USFS QUAC data file 
USFS_QUAC_df <- read.csv("CSV_Files/USFS_QUAC_rebinned_df.csv")
USFS_QUAC_df <- USFS_QUAC_df[1:206,]
#pop name list
pop_list <- unique(USFS_QUAC_df$Pop)

#name populations 
levels(USFS_QUAC_gen@pop) <- pop_list 

###########################################################
#     Run clone check and remove individuals with MD      #
###########################################################

##run clone check 
#convert genind object to a genclone object 
QUAC_genclone <- as.genclone(USFS_QUAC_gen)

#identify multi-locus genotypes (non-clones)
QUAC_MLG <- mlg.id(QUAC_genclone)

#create clone index 
QUAC_clone_index <- which(sapply(QUAC_MLG, function(x) length(x)>1))

#list clones
QUAC_clone_list <- list()

#create a list of the clone individuals 
for(c in 1:length(QUAC_clone_index)) QUAC_clone_list[[c]] <- QUAC_MLG[[QUAC_clone_index[[c]]]]

#remove any individuals with more than 25% missing data 
QUAC_nomd_gen <- missingno(USFS_QUAC_gen, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) 

#save data file without individual with too much missing data
genind2genalex(QUAC_nomd_gen, "CSV_Files/USFS_QUAC_nomd_genalex.csv")
