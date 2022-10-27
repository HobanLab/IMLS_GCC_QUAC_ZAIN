##This script calculates pairwise between wild Q. acerifolia and Z. 
#integrifolia individuals and their physical distance. 
#First we calculate the Fst between wild QUAC and ZAIN populations 
#and then regress it with the distance between populations. 
#We use "cleaned" data files here, which means genetic 
#data files have been cleaned for clones and individuals with 
#too much missing data (25% missing data)

#########################
#        Libraries      #
#########################

library(adegenet)
library(hierfstat)
library(geosphere)
library(parallel)
library(foreach)
library(doParallel)

#########################
#   Load Data Files     #
#########################
#set working directory
setwd("../../Data_Files")

#list data files 
#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames/", pattern = "_clean_df.csv")

#wild df coords 
sp_coord_df <- list.files(path = "Data_Frames/", pattern = "wild_coord_df")

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#############################
#     Fst Calculations      #
#############################
#data frame for the relationship 
sp_IBD_df <- matrix(nrow = length(scenario_list), ncol = 2)

#list to store genind objects 
sp_wild_genind_list <- list()

##loop to calculate pairwise Fst between wild pops for all scenarios 
for(sp in 1:length(scenario_list)){
  
  #read in genepop file as a genind object
  sp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  #read in data frame 
  sp_df <- read.csv(paste0("Data_Frames/",sp_df_list[[sp]]))
  
  #name rows with individual names
  rownames(sp_genind@tab) <- sp_df[,1]
  
  #create list with population names 
  sp_pop_names <- unique(sp_df[,2])
  
  #name populations 
  levels(sp_genind@pop) <- sp_pop_names
  
  #limit data frame object to wild individuals  
  sp_wild_df <- sp_df[sp_df[,3] == "Wild",]
  
  #wild population name list 
  sp_wild_pop_names <- unique(sp_wild_df[,2])
  
  #limit genind object to wild individuals
  sp_wild_genind_list[[sp]] <- sp_genind[rownames(sp_genind@tab) %in% sp_wild_df[,1],]
  
  #use parallelization because this step takes quite a while without
  #Calculate the number of cores
  cores <- detectCores() - 1
  #Initiate cluster
  cl <- makeCluster(cores)
  
  #convert genind object to 
  sp_hierfstat_list <- parLapply(cl, sp_wild_genind_list, genind2hierfstat)
  
  #run pairwise Fst code 
  sp_pwfst_df <- parLapply(cl, sp_hierfstat_list, pairwise.neifst)
  
  #stop clustering
  stopCluster(cl)
  
  #write out pairwise fst df 
  write.csv(sp_pwfst_df[[sp]], paste0("../Analyses/Results/Sum_Stats/", scenario_list[[sp]], "_pwfst_df.csv"))

}

#create list to store distance data frames
sp_dist_df_list <- list()

#loop to calculate pairwise distance between wild pops 
for(sp in 1:length(sp_coord_df)){
 
  #load in the data file for the coordinates
  sp_wild_coord_df <- read.csv(paste0("Data_Frames/", sp_coord_df[[sp]]))
    
  #wild population names 
  sp_wild_pop_names <- unique(sp_wild_coord_df[,2])
    
  #calculate mean longitude and latitude for each population
  sp_mean_lon <- matrix()
  sp_mean_lat <- matrix()
    
  #loop to calculate mean longitude by pop
  for(pop in sp_wild_pop_names){
    
    #run code to calculate mean longitude by pop 
    sp_mean_lon[pop] <- mean(sp_wild_coord_df[sp_wild_coord_df[,2] == pop,][,3])  
    
  }
    
   #loop to calculate the mean lat by population 
   for(pop in sp_wild_pop_names){
    
     #run code to calculate mean latitude by pop
    sp_mean_lat[pop] <- mean(sp_wild_coord_df[sp_wild_coord_df[,2] == pop,][,4])  
    
   }
    
  #now combine into one data frame for mean lon/mean lat
  sp_com_coord_df <- cbind(sp_mean_lon, sp_mean_lat)[-1,]
      
  #create data frame to store pairwise distance between pops  
  sp_dist_df <- matrix(nrow = length(sp_wild_pop_names), 
                         ncol = length(sp_wild_pop_names))
    
  #loop to calculate pairwise distance 
  for(r1 in 1:length(sp_wild_pop_names)){ #compare lon/lat pop 1
     for(r2 in 1:length(sp_wild_pop_names)){ #compare lon/lat pop 2
      
       #code to compare r1 - pop 1 and r2 - pop 2 to calculate distance between
       sp_dist_df[r1,r2] <-  distm(sp_com_coord_df[r1,], sp_com_coord_df[r2,], fun = distGeo)/1000
      
       }
    }

  #list that stores pairwise comparisons 
  sp_dist_df_list[[sp]] <- sp_dist_df

}

#df_list[[2]][is.na(df_list[[2]])] <- 0

#df_list[[2]][lower.tri(df_list[[2]])]
#summary(lm(df_list[[2]][lower.tri(df_list[[2]])]~df_list[[1]][lower.tri(df_list[[1]])]))

##compare pairwise distance and pairwise Fst between pops
#create a matrix to store final results 
sp_fst_dist_df <- matrix(nrow = length(scenario_list), ncol = 2)

#loop to compare pairwise fst and distance 
for(sp in 1:length(sp_coord_df)){
  
  #first, replace NA values with zeroes in the pwfst dfs
  sp_pwfst_df[[sp]][is.na(sp_pwfst_df[[sp]])] <- 0
  
  #create regression 
  reg <- lm(sp_pwfst_df[[sp]][lower.tri(sp_pwfst_df[[sp]])]~sp_dist_df_list[[sp]][lower.tri(sp_dist_df_list[[sp]])])
  
  #summarize
  #sp_reg_sum <- summary(reg)
  
}

  #create geographic distance matrix 
  #sp_dist_df <- matrix(nrow = length(sp_wild_pop_names), ncol = length(sp_wild_pop_names))
  
  #replace NAs with 0s in the PW Fst data frame 
sp_pwfst_df[[1]][is.na(sp_pwfst_df[[1]])] <- 0
  
  ##linear regression between fst and distance
  #reg <- lm(sp_pwfst_df[[sp]][lower.tri(sp_pwfst_df[[sp]])] ~ sp_dist_df[lower.tri(sp_dist_df)])
  
  #saving summary statistics - R2 and p-value
  #sp_IBD_df[sp,1] <- as.numeric(summary(reg)[9])
  #sp_IBD_df[sp,2] <- summary(reg)$coefficients[2,4]
  
  #rp_IBD <- vector('expression',2)
  #rp_IBD[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
  #                        list(MYVALUE = format(sp_IBD_df[sp,1],dig=3)))[2]
  #rp_IBD[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
  #                       list(MYOTHERVALUE = format(sp_IBD_df[sp,2], digits = 2)))[2]
  
  #write out data frame 
  #pdf(paste0("../Analyses/Results/Sum_Stats/", scenario_list[[sp]], ".pdf"), width = 10, height = 8)
  
  #plot(sp_pwfst_df~sp_dist_df, pch = 16, xlim = c(0, max(sp_dist_df)), ylim = c(0,max(sp_pwfst_df)),
  #     xlab = "Distance (km)", ylab = "PW Fst")
  #abline(reg)
  
  #legend('topleft', legend = rp_IBD, bty = 'n')
  
  #dev.off

rownames(sp_IBD_df) <- scenario_list
colnames(sp_IBD_df) <- c("R2", "p-value")
write.csv(sp_IBD_df,"../Analyses/Results/Sum_Stats/sp_IBD_df.csv")
