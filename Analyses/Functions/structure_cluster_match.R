###This function is designed to match the cluster within structure to each wild
##population. It requires a Q matrix that results from Structure and only wild 
#populations to assign a cluster value to wild populations.  

###########################
#     Structure match     #
###########################

#structure cluster match function 
str_clust_match <- function(x){

  #first, if there are rows with NA in the structure output, remove them 
  x <- na.omit(x)
  
  ##create columns that will save the clusters and match assignments 
  #create a column for the cluster result from structure
  x$Cluster <- NA
  #create a column cluster with a population attached column 
  x$cluster_pop <- NA
  #create a column for matching  
  x$cluster_match <- NA

  ##loop to assign structure clusters to their referenced wild population 
  #for each row of the data frame, compare the cluster 
  for(c in 1:nrow(x)){
  
    #first, determine the structure cluster for each wild population
    #each structure output file should have the first column for individual name
    #the second column is populations
    #following that, all the numbers are cluster assignments 
    #3 is the first cluster column and the last 3 columns are the added columns 
    col_name <- x[c,3:ncol(x)-3] > 0.6
    
    #if there is a "true" value, that means there is a major cluster for the pop
    #assign it the name of the cluster that it is 
    if(any(col_name) == TRUE) {
    
        x$Cluster[c] <- paste0("Cluster", which(col_name))
    
      #if there is no true value, that means it was not assigned to a major cluster
      #name these cases 'none'
      }else{
    
        x$Cluster[c] <- "None"
    
      }
    }
    
    #now create a list of all the populations 
    pop_list <- unique(QUAC_k4$Pop)
  
    cluster_pop <- list()

#create a column to link populations to each cluster 
for(cl in 1:length(pop_list)){
  
  cluster_pop[[cl]] <- names(which.max(table(QUAC_k4[QUAC_k4$Pop == pop_list[[cl]],]$Cluster)))
  
  
  QUAC_k4[QUAC_k4$Pop == pop_list[[cl]],]$cluster_pop <- rep(paste0(cluster_pop[[cl]], "_", pop_list[[cl]]), 
                                                             nrow(QUAC_k4[QUAC_k4$Pop == pop_list[[cl]],]))
  
  #calculate if the match is correct 
  QUAC_k4$cluster_match <- gsub("_.*", "", QUAC_k4$cluster_pop) == QUAC_k4$Cluster
  
}
}