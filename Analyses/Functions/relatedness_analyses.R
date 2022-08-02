###This function is used in the GCC QUAC ZAIN project to identify the number of
##siblings in genetic data. This is accomplished using the Loiselle statistic
#which corrects for small population sizes. Half siblings are related by 12.5%
#or more, while full siblings are related by 25% or more using this statistic. 
#This is accomplished in different ways depending on the "population type."
#Botanic gardens use all of the wild populations as their "reference populations"
#while wild populations use all other wild populations as their reference population. 

#########################
#        Libraries      #
#########################

library(Demerelate)

#########################################
#           Relatedness Df Code         #
#########################################
##halfsibs
#function to output the % of halfsibs in all garden pops 
halfsib_loiselle_sum_garden <- function(x){
    
    #first need to run the relatedness analysis with just garden/wild designation
    relatedness_df <- Demerelate(x[,-2], object = T, value = "loiselle", NA.rm	= TRUE)

    #next, determine the names of halfsibs
    halfsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness$Garden) > 0.125))
        
    #clean front characters
    halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
        
    #clean the duplicate name 
    halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
        
    #create list of unique individuals greater than and equal to half-sib relatedness
    halfsib_list <- unique(halfsibs_clean_back)
  
    #calculate percent of half-sibs 
    halfsib_garden_sum <- length(halfsib_list)/length(x[x[,3] == "Garden",][,1])
    #name columns
    return(halfsib_garden_sum)
}

#function to output the % of halfsibs in wild populations 
halfsib_loiselle_sum_wild <- function(x){
    
    #limit data frame to just wild populations 
    relate_wild_df <- x[x[,3] == "Wild",]
    
    #run relatedness analysis with just wild population names - remove pop type
    relatedness_df <- Demerelate(relate_wild_df[,-3], object = T, value = "loiselle", NA.rm	= TRUE, )
    
    #create a population name list for each data frame 
    pop_names <- unique(x[x[,3] == "Wild",][,2])
    
    #create a matrix to store # of related individuals
    relate_pop_df <- matrix(nrow = length(pop_names), ncol = 1)
    
    #then create a loop to take the name of every population and assess the level of relatedness in each pop 
    for(pop in 1:length(pop_names)){
        
        #next, determine the names of halfsibs
        halfsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness[pop_names[[pop]]]) > 0.125))
        
        #now clean the front 
        halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
        
        #clean the back for the list of halfsibs
        halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
        
        #create list of half-sibs 
        halfsib_list <- unique(halfsibs_clean_back)
        
        #create a list with all of the numbers of related individuals 
        relate_pop_df[pop,1] <- length(halfsib_list)
        
    }
    
    #calculate percent of half sibs  
    halfsib_wild_relate_sum <- sum(relate_pop_df)
    #save output 
    halfsib_wild_sum <- halfsib_wild_relate_sum/length(relate_wild_df[,1])
    #name columns 
    return(halfsib_wild_sum)
}

##full-sibs 
#function to output the % of full sibs in all garden pops 
fullsib_loiselle_sum_garden <- function(x){
  
  
  #run relatedness analysis on individuals with only pop type - garden or wild 
  relatedness_df <- Demerelate(x[,-2], object = T, value = "loiselle", NA.rm	= TRUE)
  
  #next, determine the names of full sibs 
  fullsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness$Garden) > 0.25))
  
  #now clean the front 
  fullsibs_clean_front <- gsub("^.*\\.","", fullsibs_names)
  
  #clean the back for the list of sibs
  fullsibs_clean_back <- gsub("^.*\\_","", fullsibs_clean_front)
  
  #create list of sibs 
  fullsib_list <- unique(fullsibs_clean_back)
  
  #calculate percent of sibs 
  fullsib_garden_sum <- length(fullsib_list)/length(x[x[,3] == "Garden",][,1])
  #return percent of sibs  
  return(fullsib_garden_sum)
}

#function to output the % of full sibs in wild populations 
fullsib_loiselle_sum_wild <- function(x){
  
  #limit data frame to just wild individuals 
  relate_wild_df <- x[x[,3] == "Wild",]
  
  #run relatedness analysis with wild individuals 
  relatedness_df <- Demerelate(relate_wild_df[,-3], object = T, value = "loiselle", NA.rm	= TRUE)
  
  #create a population name list for each data frame 
  pop_names <- unique(relate_wild_df[,2])
  
  #create a matrix for relatedness 
  relate_pop_df <- matrix(nrow = length(pop_names), ncol = 1)
  
  #assess the level of relatedness in wild each pop 
  for(pop in 1:length(pop_names)){
    
    #next, determine the names of halfsibs
    fullsibs_names <- names(which(unlist(relatedness_df$Empirical_Relatedness[pop_names[[pop]]]) > 0.25))
    
    #now clean the front 
    fullsibs_clean_front <- gsub("^.*\\.","", fullsibs_names)
    
    #clean the back for the list of sibs
    fullsibs_clean_back <- gsub("^.*\\_","", fullsibs_clean_front)
    
    #create list of full sibs 
    fullsib_list <- unique(fullsibs_clean_back)
    
    #create data frame by pop with the total individual numbers and the # of full sibs 
    relate_pop_df[pop,1] <- length(fullsib_list)
    
  }
  
  #calculate percent of full sibs  
  fullsib_wild_relate_sum <- sum(relate_pop_df)
  #save output 
  fullsib_wild_sum <- fullsib_wild_relate_sum/length(relate_wild_df[,1])
  #return percent of full siblings 
  return(fullsib_wild_sum)
}

