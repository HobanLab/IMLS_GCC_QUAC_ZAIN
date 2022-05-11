#########################
#        Libraries      #
#########################

library(dplyr)

#########################################
#           Relatedness Df Code         #
#########################################
mat_accession <- function(x){

  #add a column with cleaned accession names 
  sp_garden_accessions <- x %>% mutate(accession_simple = gsub("\\*.*","",x$Accession))

  #create a list of all the botanic gardens for each species 
  sp_bg_names <- unique(x[,2])

  #create data frame to store results 
  sp_bg_maternal <- matrix(nrow = length(sp_bg_names), ncol = 1)

  #code calculate maternal lines by garden 
  for(garden in 1:length(sp_bg_names)){
  
    #count unique maternal lines for each botanic garden
    sp_bg_maternal[garden,1] <- length(unique(sp_garden_accessions[sp_garden_accessions[,2] == paste0(sp_bg_names[[garden]]),]$accession_simple))
  
  }
  rownames(sp_bg_maternal) <- sp_bg_names
  colnames(sp_bg_maternal) <- "mat_acc"
  return(sp_bg_maternal)
}
