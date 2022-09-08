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
halfsib_garden_loiselle_list <- function(x){
    
    #next, determine the names of halfsibs
    halfsibs_names <- names(which(unlist(x$Empirical_Relatedness[1]) > 0.125))
        
    #clean front characters
    halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
        
    #clean the duplicate name 
    halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
        
    #create list of unique individuals greater than and equal to halfsib relatedness
    halfsib_list <- unique(halfsibs_clean_back)
    
    #return the list of halfsibs
    return(halfsib_list)
}

##halfsibs
#function to output the % of halfsibs in all wild pops 
halfsib_wild_loiselle_list <- function(x){
  
  #next, determine the names of halfsibs
  halfsibs_names <- names(which(unlist(x$Empirical_Relatedness) > 0.125))
  
  #clean front characters
  halfsibs_clean_front <- gsub("^.*\\.","", halfsibs_names)
  
  #clean the duplicate name 
  halfsibs_clean_back <- gsub("^.*\\_","", halfsibs_clean_front)
  
  #create list of unique individuals greater than and equal to halfsib relatedness
  halfsib_list <- unique(halfsibs_clean_back)
  
  #return the list of halfsibs
  return(halfsib_list)
}

##full-sibs 
#function to output the % of full sibs in all garden pops 
fullsib_garden_loiselle_list <- function(x){
  
  #next, determine the names of full sibs 
  fullsibs_names <- names(which(unlist(x$Empirical_Relatedness[1]) > 0.25))
  
  #now clean the front 
  fullsibs_clean_front <- gsub("^.*\\.","", fullsibs_names)
  
  #clean the back for the list of sibs
  fullsibs_clean_back <- gsub("^.*\\_","", fullsibs_clean_front)
  
  #create list of sibs 
  fullsib_list <- unique(fullsibs_clean_back)
  
  #return percent of sibs  
  return(fullsib_list)
}

#function to output the % of full sibs in all wild pops 
fullsib_wild_loiselle_list <- function(x){
  
  #next, determine the names of full sibs 
  fullsibs_names <- names(which(unlist(x$Empirical_Relatedness) > 0.25))
  
  #now clean the front 
  fullsibs_clean_front <- gsub("^.*\\.","", fullsibs_names)
  
  #clean the back for the list of sibs
  fullsibs_clean_back <- gsub("^.*\\_","", fullsibs_clean_front)
  
  #create list of sibs 
  fullsib_list <- unique(fullsibs_clean_back)
  
  #return percent of sibs  
  return(fullsib_list)
}
