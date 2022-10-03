#This script details the analyses performed to determine which relatedness 
#analysis performs the best on these data. 

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

#load in relatedness function 
#source("../Analyses/RScripts/relatedness_analyses.R")
#################################
#      Relatedness Analysis     #
#################################

##first, load in one data file and run relatedness analysis to determine 
#how different population structures changes levels of relatedness 
#detected. 

#load in QUAC data file 
QUAC_rel_df <- read.csv("Data_Frames/QUAC_allpop_clean_df.csv")

##run relatedness analysis on different edits of the data file 








#list clean data frames 
sp_clean_df_list <- list.files(path = "Data_Frames", pattern = "clean_df.csv")

#list clean genepop files 
sp_clean_genepop_list <- list.files(path = "Adegenet_Files", pattern = "allpop_clean.gen")

#list of data frames
sp_clean_temp_df <- list()

#list of the genind objects
sp_clean_temp_gen <- list()

#create a list of the relatedness tests 
relatedness_analyses_list <- c("loiselle", "wang", "ritland")

#save summary df 
halfsib_df <- as.data.frame(matrix(nrow = length(species_list), ncol = length(relatedness_analyses_list)*2))

##loop over species data frames 
for(sp in 1:length(species_list)){
  ##loop over type of relatedness analysis 
  for(relate in 1:length(relatedness_analyses_list)){
    #load in clean data frames to perform relatedness analysis on 
    sp_clean_temp_df[[sp]] <- read.csv(paste0("Data_Frames/",sp_clean_df_list[[sp]]))
    
    #load in genepop files into genind 
    sp_clean_temp_gen[[sp]] <- read.genepop(paste0("Adegenet_Files/",sp_clean_genepop_list[[sp]]), ncode = 3)
    #name individuals in genind 
    rownames(sp_clean_temp_gen[[sp]]@tab) <- sp_clean_temp_df[[sp]][,1]
    #name pops 
    levels(sp_clean_temp_gen[[sp]]@pop) <- unique(sp_clean_temp_df[[sp]][,2])
    
    ##Run relate red code 
    #Garden
    sp_garden_clean_temp_df <- sp_clean_temp_df[[sp]][,-2]
    
    #limit data frame by garden only 
    sp_garden_only_temp_df <- sp_garden_clean_temp_df[sp_garden_clean_temp_df[,2] == "Garden",]
    
    #run relatedness analysis 
    sp_garden_rel_df <- Demerelate(sp_garden_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #use function to limit individuals  
    sp_garden_halfsib_list <- half_sibling_reduction(sp_garden_rel_df)
    
    #reduce data frame by list of half-siblings 
    sp_garden_relate_red_df <- sp_garden_clean_temp_df[!sp_garden_clean_temp_df[,1] %in% sp_garden_halfsib_list,]
    
    #write garden df 
    sp_garden_relate_red_df <- sp_garden_relate_red_df[sp_garden_relate_red_df[,2] == "Garden",]
    
    #save in df 
    halfsib_df[sp, relate] <- paste0(signif(((length(sp_garden_only_temp_df[,1])-length(sp_garden_relate_red_df[,1]))/
                                               length(sp_garden_only_temp_df[,1])),3)*100, "% (", 
                                     length(sp_garden_relate_red_df[,1]),")")
    
    #write out 
    write.csv(sp_garden_relate_red_df, paste0("Data_Frames/Relate_Red/",  species_list[[sp]], 
                                              "_Garden_relate_red_", relatedness_analyses_list[[relate]], "_df.csv"),
              row.names = FALSE)
    
    ##Wild relatedness analysis 
    sp_wild_clean_temp_df <- sp_clean_temp_df[[sp]][sp_clean_temp_df[[sp]][,3] == "Wild",][,-3]
    
    #now run relatedness analysis 
    sp_wild_rel_df <- Demerelate(sp_wild_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #limit by analysis halfsibs 
    sp_wild_halfsib_list <- half_sibling_reduction(sp_wild_rel_df)
    
    #reduce data frame by list of half-siblings 
    sp_wild_relate_red_df <- sp_wild_clean_temp_df[!sp_wild_clean_temp_df[,1] %in% sp_wild_halfsib_list,]
    
    #number of ind removed 
    halfsib_df[sp, relate+3] <-paste0(signif(((length(sp_wild_clean_temp_df[,1])-length(sp_wild_relate_red_df[,1]))/
                                                length(sp_wild_clean_temp_df[,1])),3)*100, "% (", 
                                      length(sp_wild_relate_red_df[,1]), ")")
    
    #write csv of the reduced data frame 
    write.csv(sp_wild_relate_red_df, paste0("Data_Frames/Relate_Red/", species_list[[sp]], 
                                            "_Wild_relate_red_", relatedness_analyses_list[[relate]], "_df.csv"),
              row.names = FALSE)
    
  }
}

#write out summary table 
rownames(halfsib_df) <- species_list
colnames(halfsib_df) <- rep(relatedness_analyses_list, 2)
write.csv(halfsib_df, "../Analyses/Results/Sum_Stats/halfsib_df.csv")