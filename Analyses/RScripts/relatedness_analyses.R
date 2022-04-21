#########################
#        Libraries      #
#########################

library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("C:/Users/eschumacher/Documents/GitHub/QUAC_insitu_exsitu/QUAC_data_files/QUAC_data_frames")

#load relatedness df 
QUAC_df <- read.csv("Garden_Wild/QUAC_clean_df.csv")

#limit by garden 
QUAC_garden_df <- QUAC_df[QUAC_df$Garden_Wild == "Garden",]
QUAC_garden_df <- QUAC_garden_df[,-1]

##Function to reduce data frames by half-sibling relatedness
halfsib_relatedness_reduction_loiselle <- function(x){
    
    #first need to run the relatedness analysis 
    relatedness_df <- Demerelate(x, object = T, value = "loiselle")
    
    #next, determine the names of halfsibs
    halfsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness) > 0.25))
    
    #now clean the front 
    halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
    
    #clean the back for the list of halfsibs
    halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
    
    #create list of halfsibs 
    halfsib_list <- unique(halfsibs_clean_back)
    
    #Now create a data frame reduced by these things 
    #relate_red_df <- x[!x[,1] %in% halfsib_list,]
    
    
}
