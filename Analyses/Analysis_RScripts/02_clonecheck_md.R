###This script details the analyses run on Quercus acerifolia and Zamia integrifolia 
##(referred to as QUAC and ZAIN in this script for brevity) to prep individuals
#for genetic analyses. We remove individuals that are clones and with more than
#25% missing data, creating "clean" genotype files (genepop and tables), which all
#subsequent genetic analyses are run on.

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../Data_Files")

#now read in genepop file as a genind for adegenet 
sp_genind <- list.files(path = "Adegenet_Files", pattern = "allpop.gen")

#load relatedness CSV with ind and pop info
sp_df <- list.files(path = "CSV_files", pattern = "allpop_df.csv")

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

############################################################
#    Remove Clones and Individuals with Missing Data       #
############################################################

##loop to remove clones and individuals with missing data for genetic analyses
#loops over all scenarios for both species - QUAC and ZAIN 
for(sp in 1:length(sp_genind)){
  
  #load in genepop file as a genind object 
  sp_temp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind[[sp]]), ncode = 3)
  
  #load in genetic data frame  
  sp_temp_df <- read.csv(paste0("CSV_Files/",sp_df[[sp]]))
  
  #name rows in the genind object 
  rownames(sp_temp_genind@tab) <- sp_temp_df[,1]
  
  #create pop names list
  pop_names <- unique(sp_temp_df[,2])
  
  #name pops in genind object 
  levels(sp_temp_genind@pop) <- pop_names
  
  ##run clone check 
  #convert genind object to a genclone object 
  sp_temp_genclone <- as.genclone(sp_temp_genind)
  
  #identify multi-locus genotypes (non-clones)
  sp_mlg <- mlg.id(sp_temp_genclone)
  
  #create clone index 
  sp_clone_index <- which(sapply(sp_mlg, function(x) length(x)>1))
  
  #list clones
  sp_clone_id <- list()
  
  #create a list of the clone individuals 
  for(clones in 1:length(sp_clone_index)) sp_clone_id[[clones]] <- sp_mlg[[sp_clone_index[[clones]]]]
  
  #then create genind objects without clones 
  sp_genind_nocl <- clonecorrect(sp_temp_genind)
  
  #remove any individuals with more than 25% missing data 
  sp_genind_nomd <- missingno(sp_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) 
  
  #reduce data frame by nomd and clones
  sp_clean_df <- sp_temp_df[sp_temp_df[,1] %in% rownames(sp_genind_nomd@tab),]
  
  ##write out files cleaned for missing data and clones 
  #write out genalex
  genind2genalex(sp_genind_nomd, file = paste0("CSV_Files/", gsub("\\..*","",sp_genind[[sp]]), "_clean_genalex.csv"), 
              overwrite = TRUE)
  #write out data frame
  write.csv(sp_clean_df, paste0("CSV_Files/", gsub("\\..*","",sp_genind[[sp]]), "_clean_df.csv"), row.names = FALSE)

  #create an individual summary data frame 
  sp_ind_df <- cbind(summary(sp_temp_genind)[[2]], summary(sp_genind_nomd)[[2]])
  #name columns 
  colnames(sp_ind_df) <- c("before_datacleaning", "after_datacleaning")
  
  #write out summary individual data frame 
  write.csv(sp_ind_df, paste0("../Analyses/Results/Sum_Stats/", scenario_list[[sp]], "_ind_n.csv"))
  
}

