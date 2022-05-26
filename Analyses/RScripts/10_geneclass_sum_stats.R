##This script checks the number of successes with the Geneclass assignment 
#test. 

###########################
#     Load Data Files     #
###########################
setwd("G:/Shared drives/Emily_Schumacher/Project_Folders/GCC_QUAC_ZAIN/Data_Files/Assignment_Testing")

#load in all geneclass output files 
sp_geneclass_output <- list.files(pattern = "output.csv")

###################################
#     Assess assignment test      #
###################################
#create final matrix to store results
assign_success_df <- matrix(nrow = length(sp_geneclass_output), ncol = 1)

#create a column with a match area 
for(df in 1:length(sp_geneclass_output)){
  
  #load in the output df 
  sp_geneclass_df <- read.csv(sp_geneclass_output[[df]])
  
  ##data frame cleaning steps 
  #first, remove rows that do not have a wild source pop 
  sp_geneclass_df <- na.omit(sp_geneclass_df)
  
    #write code to remove 'NWS' - not wild sampled - individuals from the 
    #ZAIN data files 
    if(df == 9|df == 10|df == 11|df == 12|df == 13|df == 14| df ==15|df == 16){
    
      sp_geneclass_df <- sp_geneclass_df[-grep("NWS", sp_geneclass_df$Wild_Source),]
    }
  
    #loop to identify matching source populations for each botanic garden
    #individual
    for(atest in 1:length(sp_geneclass_df[,1])){
     
      #create a column that saves the number of successes for each
      sp_geneclass_df$match_TF[atest] <- sp_geneclass_df$First_Pop[atest] ==
                                        sp_geneclass_df$Wild_Source[atest]
    
    }
  
  #create a data frame that saves a success and failure
  assign_success <- length(sp_geneclass_df[sp_geneclass_df$match_TF == "TRUE",][,1])
  
  #total individual character 
  assign_total <- length(sp_geneclass_df[,1])
  
  #final data file with percent of success of geneclass assignment
  assign_success_df[df,] <- signif((assign_success/assign_total)*100, 3)
}

#name rows with the test function
rownames(assign_success_df) <- strsplit(sp_geneclass_output, split='_output.csv', fixed=TRUE)
#write out data frame
write.csv(assign_success_df, "assign_success_df.csv")
