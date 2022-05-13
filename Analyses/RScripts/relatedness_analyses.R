#########################
#        Libraries      #
#########################

library(Demerelate)

#########################################
#           Relatedness Df Code         #
#########################################
#function to output data frames reduced by half-sib relatedness - garden dfs
halfsib_loiselle_garden_sum_bypop_df <- function(x){
    
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
        
        #create data frame that indicates the # of unique maternal accessions -
        #no halfsibs
        relate_pop_df[pop,1] <- length(x[x[,2] == pop_names[[pop]],][,1])
        relate_pop_df[pop,2] <- relate_pop_df[pop,1] - length(halfsib_list)
        
     }
    rownames(relate_pop_df) <- pop_names
    colnames(relate_pop_df) <- c("Tot_Ind", "Half_Sibs")
    return(relate_pop_df)
}

#function to output the % of halfsibs in each pop 
halfsib_loiselle_sum_garden_df <- function(x){
    
    #create data frame for the final output 
    allgarden_sum_df <- matrix(nrow = 1, ncol = 2)
    
    #first need to run the relatedness analysis 
    relatedness_df <- Demerelate(x[,-2], object = T, value = "loiselle", NA.rm	= TRUE)

    #next, determine the names of halfsibs
    halfsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness$Garden) > 0.125))
        
    #now clean the front 
    halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
        
    #clean the back for the list of halfsibs
    halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
        
    #create list of half-sibs 
    halfsib_list <- unique(halfsibs_clean_back)
  
    #calculate percent of halfsibs 
    allgarden_sum_df <- cbind((length(halfsib_list)/length(x[x[,3] == "Garden",][,1])), length(x[x[,3] == "Garden",][,1]))
    #name columns 
    colnames(allgarden_sum_df) <- c("Per_Halfsibs", "Tot_Ind")
    return(allgarden_sum_df)
}

#function to output the % of halfsibs in each pop 
halfsib_loiselle_sum_wild_df <- function(x){
    
    #create data frame for the final output 
    allwild_sum_df <- matrix(nrow = 1, ncol = 2)
    
    relate_wild_df <- x[x[,3] == "Wild",]
    
    #first need to run the relatedness analysis 
    relatedness_df <- Demerelate(relate_wild_df[,-3], object = T, value = "loiselle", NA.rm	= TRUE)
    
    #create a population name list for each data frame 
    pop_names <- unique(x[x[,3] == "Wild",][,2])
    
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
        
        #create list of half-sibs 
        halfsib_list <- unique(halfsibs_clean_back)
        
        #create data frame by pop with the total individual numbers and the # of halfsibs 
        relate_pop_df[pop,1] <- length(x[x[,2] == pop_names[[pop]],][,1])
        relate_pop_df[pop,2] <- length(halfsib_list)
        
    }
    
    #calculate percent of half sibs  
    allwild_relate_sum <- colSums(relate_pop_df)
    #save output 
    allwild_sum_df <- cbind(allwild_relate_sum[2]/allwild_relate_sum[1], allwild_relate_sum[1])
    #name columns 
    colnames(allwild_sum_df) <- c("Per_Halfsibs", "Tot_Ind")
    return(allwild_sum_df)
}
