#########################
#        Libraries      #
#########################

library(Demerelate)

#########################################
#           Relatedness Df Code         #
#########################################
##Function to reduce data frames by half-sibling relatedness
halfsib_relate_df_loiselle <- function(x){
    
    x <- x[order(x[,2]),]

    #first need to run the relatedness analysis 
    relatedness_df <- Demerelate(x[,-3], object = T, value = "loiselle", NA.rm	= TRUE)
    
    #create a population name list for each data frame 
    pop_names <- unique(x[x[,3] == "Garden",][,2])
    
    #create a matrix 
    relate_pop_df <- matrix(nrow = length(pop_names), ncol = 2)
    
    #then create a loop to take the name of every population and 
    #assess the level of relatedness in each pop or botanic garden 
    for(pop in 1:length(pop_names)){
        
        #next, determine the names of halfsibs
        halfsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness[pop_names[[pop]]]) > 0.125))
        
        #now clean the front 
        halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
        
        #clean the back for the list of halfsibs
        halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
        
        #create list of halfsibs 
        halfsib_list <- unique(halfsibs_clean_back)
        
        #create data frame with 
        relate_pop_df[pop,1] <- length(x[x[,2] == pop_names[[pop]],][,1])
        relate_pop_df[pop,2] <- relate_pop_df[pop,1] - length(halfsib_list)
        
     }
    rownames(relate_pop_df) <- pop_names
    colnames(relate_pop_df) <- c("Tot_Ind", "Half_Sibs")
    return(relate_pop_df)
}


