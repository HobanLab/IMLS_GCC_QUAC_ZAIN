##This script details an assessment of the number of individuals
#sourced from the same maternal plant in botanic garden 
#populations.

#########################
#        Libraries      #
#########################

library(tidyverse)
library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../Data_Files")

#load in accession records 
sp_accession_list <- list.files(path = "Data_Frames", 
                                pattern = "accession")

sp_df_list <- c("QUAC_woK_allpop_clean_df.csv","ZAIN_rebinned_allpop_clean_df.csv")

sp_list <- c("QUAC", "ZAIN")

###################################
#     Maternal plant records      #
###################################
#loop over scenarios to calculate the relatedness levels based on half siblings
for(m in 1:length(sp_df_list)){
  #loop over maternal accession information 
  for(a in 1:length(sp_accession_list)){

  #load in data frame for relatedness analysis  
  sp_temp_df <- read.csv(paste0("Data_Frames/",sp_df_list[[m]]))
  
  #load in data frame with maternal accession information 
  sp_accessions <- read.csv(paste0("Data_Frames/",sp_accession_list[[a]]))

  #run relatedness code to generate relatedness data frame 
  relate_df <- halfsib_relate_df_loiselle(sp_temp_df)
  
  #limit accession data by cleaned data 
  sp_temp_garden_df <- sp_temp_df[sp_temp_df$Garden_Wild == "Garden",]
  
  if(m == 1 & a == 1){
    
    sp_temp_accession_df <- sp_accessions[sp_accessions[,1] %in% sp_temp_garden_df[,1],]
    
    #run maternal accession code 
    mat_df <- accession_count(sp_temp_accession_df)
    
  rel_mat_df <- cbind(relate_df, mat_df)
  
  write.csv(rel_mat_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rel_mat_df.csv")
  }
  if(m == 2 & a == 2){

    sp_temp_accession_df <- sp_accessions[sp_accessions[,1] %in% sp_temp_garden_df[,1],]
    
    #run maternal accession code 
    mat_df <- accession_count(sp_temp_accession_df)
    rel_mat_df <- cbind(relate_df, mat_df)
    
    write.csv(rel_mat_df, "../Analyses/Results/Garden_Wild_Comparison/ZAIN_rel_mat_df.csv")
    
  }

  }

}
