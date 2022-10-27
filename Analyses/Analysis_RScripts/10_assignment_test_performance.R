###This script is used to assess the success of tests like structure and Geneclass
##at assigning botanic garden individuals to their referenced source population. 
#The wild populations are used as the reference source populations and the 
#garden individuals are assigned to these wild populations, and then compared 
#to the source population information recorded for each individual.

#########################
#        Libraries      #
#########################

library(diveRsity)
library(adegenet)
library(poppr)
library(Demerelate)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files 
setwd("../../Analyses")

#create a list of assignment test results for geneclass
sp_geneclass_output <- list.files(path = "../Analyses/Results/Clustering/Geneclass", pattern = ".csv")

###################################
#     Assess assignment test      #
###################################
#create final matrix to store results
assign_success_df <- matrix(nrow = length(sp_geneclass_output), ncol = 3)

#create a column with a match area 
for(df in 1:length(sp_geneclass_output)){
  
  #load in the output df 
  sp_geneclass_df <- read.csv(paste0("Results/Clustering/Geneclass/",sp_geneclass_output[[df]]))
  
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
      sp_geneclass_df$match1_TF[atest] <- sp_geneclass_df$First_Pop[atest] ==
                                        sp_geneclass_df$Wild_Source[atest]
      
      #create a column that saves the number of success in the top 2 
      sp_geneclass_df$match2_TF[atest] <- sp_geneclass_df$First_Pop[atest] ==
                                          sp_geneclass_df$Wild_Source[atest] |
                                          sp_geneclass_df$Second_Pop[atest] ==
                                          sp_geneclass_df$Wild_Source[atest]
      
      #create a column to calculate amount of time assigned pop was in top 3
      sp_geneclass_df$match3_TF[atest] <- sp_geneclass_df$First_Pop[atest] ==
                                          sp_geneclass_df$Wild_Source[atest] |
                                          sp_geneclass_df$Second_Pop[atest] ==
                                          sp_geneclass_df$Wild_Source[atest] |
                                          sp_geneclass_df$Third_Pop[atest] ==
                                          sp_geneclass_df$Wild_Source[atest]
      
    }
  
  #save successes and failures 
  assign_success_df[df,1] <- signif((length(sp_geneclass_df[sp_geneclass_df$match1_TF == TRUE,][,1])/length(sp_geneclass_df$match1_TF))*100, 3)
  assign_success_df[df,2] <- signif((length(sp_geneclass_df[sp_geneclass_df$match2_TF == TRUE,][,1])/length(sp_geneclass_df$match2_TF))*100, 3)
  assign_success_df[df,3] <- signif((length(sp_geneclass_df[sp_geneclass_df$match3_TF == TRUE,][,1])/length(sp_geneclass_df$match3_TF))*100, 3)
  
}

#name rows with the test function
rownames(assign_success_df) <- strsplit(sp_geneclass_output, split='_output.csv', fixed=TRUE)
#write out data frame
write.csv(assign_success_df, "Results/Sum_Stats/geneclass_assign_success_df.csv")

###QUAC structure testing assignment 
#read in data file 
QUAC_str_assign_df <- read.csv("Results/Clustering/Structure/QUAC/QUAC_str_assignment.csv")

#clean data file for NAs
QUAC_str_assign_df <- na.omit(QUAC_str_assign_df)

#add a column for correct assignment
QUAC_str_assign_df$str_success_assign <- QUAC_str_assign_df$Str_Cluster == QUAC_str_assign_df$Wild_Source

#calculate the percent of correct assignment 
QUAC_str_assign_success <- length(QUAC_str_assign_df[QUAC_str_assign_df$str_success_assign == TRUE,][,1])/
                            length(QUAC_str_assign_df$str_success_assign)

