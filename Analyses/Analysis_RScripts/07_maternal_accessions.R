##This script details an assessment of the number of individuals
#sourced from the same maternal plant in botanic garden 
#populations.

#########################
#        Libraries      #
#########################

library(tidyverse)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../Data_Files")

#load in accession records 
sp_accession_list <- list.files(path = "CSV_Files", pattern = "accession")

#create list for species
sp_list <- c("QUAC", "ZAIN")

#load in function to calculate maternal accessions for plants 
source("../Analyses/Functions/maternal_accession.R")

###################################
#     Maternal plant records      #
###################################
#loop over maternal accession information 
for(a in 1:length(sp_accession_list)){
  
  #load in data frame with maternal accession information 
  sp_accessions <- read.csv(paste0("CSV_Files/",sp_accession_list[[1]]))

  #run maternal accession code 
  mat_df <- mat_accession(sp_accessions)
    
  #write out maternal accession data frame
  write.csv(mat_df, paste0("../Analyses/Results/Garden_Wild_Comparison/", sp_list[[a]], "_mat_acc_df.csv"))
}
