##This script details the analyses run on Quercus acerifolia and Zamia integrifolia 
#(referred to as QUAC and ZAIN in this script for brevity) genotype files 
#to prepare them for genetic diversity and structure analyses. 
#When files are referred to as "clean" that means individuals 
#that are clones and individuals with too much missing data have been removed. 
#When files and objects are titled "red" that means they have been reduced
#for half-siblings (12.5% relatedness coefficient)

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../Data_Files")

#now read in genepop file as a genind for adegenet 
sp_genind <- list.files(path = "Adegenet_Files", pattern = "allpop.gen$")

#load relatedness data frame for relatedness analysis 
sp_df <- list.files(path = "Data_Frames", pattern = "allpop_df.csv$")

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#pop type list 
pop_type <- c("Garden", "Wild")

#load in relatedness function 
source("../Analyses/RScripts/relatedness_analyses.R")

############################################################
#    Remove Clones and Individuals with Missing Data       #
############################################################

##loop to remove clones and individuals with missing data for genetic analyses
#loops over all scenarios for both species - QUAC and ZAIN 
for(sp in 1:length(sp_genind)){
  
  #load in genepop file as a genind object 
  sp_temp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind[[sp]]), ncode = 3)
  
  #load in genetic data frame  
  sp_temp_df <- read.csv(paste0("Data_Frames/",sp_df[[sp]]))
  
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
  genind2genalex(sp_genind_nomd, file = paste0("Data_Frames/", gsub("\\..*","",sp_genind[[sp]]), "_clean_genalex.csv"), 
              overwrite = TRUE)
  #write out data frame
  write.csv(sp_clean_df, paste0("Data_Frames/", gsub("\\..*","",sp_genind[[sp]]), "_clean_df.csv"), row.names = FALSE)


}

#################################
#      Relatedness Analysis     #
#################################
##code to reduce data frames by half-sibling relatedness 
#create a data frame to store results of half sib 
halfsib_garden_df <- matrix(nrow = length(scenario_list), ncol = length(pop_type))
halfsib_wild_df <- matrix(nrow = length(scenario_list), ncol = length(pop_type))

#list clean data frames 
sp_clean_df_list <- list.files(path = "Data_Frames", pattern = "clean_df.csv")

#list clean genepop files 
sp_clean_genepop_list <- list.files(path = "Adegenet_Files", pattern = "allpop_clean.gen")

#loop over species data frames 
for(sp in 1:length(scenario_list)){
  
  #load in garden data frame 
  sp_clean_temp_df <- read.csv(paste0("Data_Frames/",sp_clean_df_list[[sp]]))
  
  #load in genepop files into genind 
  sp_clean_temp_gen <- read.genepop(paste0("Adegenet_Files/",sp_clean_genepop_list[[sp]]), ncode = 3)
  #name individuals in genind 
  rownames(sp_clean_temp_gen@tab) <- sp_clean_temp_df[,1]
  #name pops 
  levels(sp_clean_temp_gen@pop) <- unique(sp_clean_temp_df[,2])
  
  ##Run relateness reduction code 
  #Garden
  halfsib_garden_df[sp,] <- halfsib_loiselle_sum_garden_df(sp_clean_temp_df)
  rownames(halfsib_garden_df) <- scenario_list
  colnames(halfsib_garden_df) <- c("Per_Halfsibs", "Tot_Ind")
  
  #Wild
  halfsib_wild_df[sp,] <- halfsib_loiselle_sum_wild_df(sp_clean_temp_df)
  rownames(halfsib_wild_df) <- scenario_list
  colnames(halfsib_wild_df) <- c("Per_Halfsibs", "Tot_Ind")
  
}

#write out summary table 
halfsib_df <- cbind(halfsib_garden_df, halfsib_wild_df)
