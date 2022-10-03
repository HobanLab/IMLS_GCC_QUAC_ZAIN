#This script details the analyses performed to determine which relatedness 
#analysis performs the best on these data. 

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
setwd("../../Data_Files")

#now read in genepop file as a genind for adegenet 
sp_genind <- list.files(path = "Adegenet_Files", pattern = "allpop.gen$")

#load relatedness data frame for relatedness analysis 
sp_df <- list.files(path = "Data_Frames", pattern = "allpop_df.csv$")

#function to output the % of sibs in all garden pops 
fullsib_loiselle_rel_fun <- function(x){
  
  #clean front characters
  fullsib_clean_front <- gsub("^.*\\.","", x)
  
  #clean the duplicate name 
  fullsib_clean_back <- gsub("^.*\\_","", fullsib_clean_front)
  
  #create list of unique individuals greater than and equal to full sib relatedness
  fullsib_list <- unique(fullsib_clean_back)
  
  #return the list of full sibs
  return(fullsib_list)
}

#################################
#      Relatedness Analysis     #
#################################

##first, load in one data file and run relatedness analysis to determine 
#how different population structures changes levels of relatedness 
#detected. 

#load in QUAC data file 
QUAC_rel_df <- read.csv("Data_Frames/QUAC_allpop_clean_df.csv")

##run relatedness analysis on different organizations of the data file
#list scenarios: 
##gwtog indicates the data file for relatedness has the wild and garden
#individuals combined in the data file
##gwsep has removed the garden or wild depending on the pop of interest
##labpopsep indicates that every population is separately named - all wild
#pops and gardens pops 
##labpopgr indicates all wild and garden individuals are grouped into
#metapopulations of just "wild" and "garden"
rel_scen <- c("gwtog_labpopsep_garden", "gwtog_labpopsep_wild", 
              "gwtog_labpopgr_garden", "gwtog_labpopgr_wild", 
              "gwsep_labpopsep_garden", "gwsep_labpopsep_wild",
              "gwsep_labpopgr_garden", "gwsep_labpopgr_wild")

#populations of gardens
garden_pops <- unique(QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Garden",]$Pop)
wild_pops <- unique(QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Wild",]$Pop)

#create matrix to store results 
rel_popstr_table <- matrix(nrow = length(rel_scen), 
                           ncol = 1)

#loop to calculate relatedness in different pop structures
for(r in 1:length(rel_scen)){
  
  #wild and garden together with pops labeled separately 
  if(r == 1|r == 2){
    
    #remove pop type column to run relatedness - wild and garden individuals grouped
    #into large metapopulations 
    QUAC_scen_rel_df <- QUAC_rel_df[,-3]
      
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[1,1] <- 277 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness[garden_pops]) > 0.25))))
    #wild relatedness
    rel_popstr_table[2,1] <- 172 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness[wild_pops]) > 0.25))))
  }
  
  #wild and garden together with pops grouped by garden and wild 
  if(r == 3|r == 4){
    
    #remove separate garden and wild populations column  
    QUAC_scen_rel_df <- QUAC_rel_df[,-2]
    
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[3,1] <- 277 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness$Garden) > 0.25))))
    #wild relatedness
    rel_popstr_table[4,1] <- 172 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness$Wild) > 0.25))))
    
  }
  
  #garden separate with pops separated 
  if(r == 5){
    
    QUAC_scen_rel_df <- QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Garden",-3]
    
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[5,1] <- 277 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness[garden_pops]) > 0.25))))
    
  }
  
  #wild separate with pops separated
  if(r == 6){
    
    QUAC_scen_rel_df <- QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Wild",-3]
    
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[6,1] <- 172 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness[wild_pops]) > 0.25))))
    
  }
  
  #garden populations pooled with just garden individuals 
  if(r == 7){
    
    #select only garden individuals and remove the sep pop labels
    QUAC_scen_rel_df <- QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Garden",-2]
    
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[7,1] <- 277 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness$Garden) > 0.25))))
    
    
  }
  
  #wild populations pooled with just wild individuals 
  if(r == 8){
    
    #select only garden individuals and remove the sep pop labels
    QUAC_scen_rel_df <- QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Wild", -2]
    
    #run relatedness analysis for garden and wild pop types 
    QUAC_relate_df <- Demerelate(QUAC_scen_rel_df, tab.dist = "NA", object = TRUE, value = "loiselle")
    
    #garden relatedness
    rel_popstr_table[8,1] <- 172 - length(fullsib_loiselle_rel_fun(names(which(unlist(QUAC_relate_df$Empirical_Relatedness$Wild) > 0.25))))
    
    
  }
  
}



#list clean data frames 
sp_clean_df_list <- list.files(path = "Data_Frames", pattern = "clean_df.csv")

#list clean genepop files 
sp_clean_genepop_list <- list.files(path = "Adegenet_Files", pattern = "allpop_clean.gen")

#list of data frames
sp_clean_temp_df <- list()

#list of the genind objects
sp_clean_temp_gen <- list()

#create a list of the relatedness tests 
relatedness_analyses_list <- c("loiselle", "wang", "ritland")

#save summary df 
halfsib_df <- as.data.frame(matrix(nrow = length(species_list), ncol = length(relatedness_analyses_list)*2))

##loop over species data frames 
for(sp in 1:length(species_list)){
  ##loop over type of relatedness analysis 
  for(relate in 1:length(relatedness_analyses_list)){
    #load in clean data frames to perform relatedness analysis on 
    sp_clean_temp_df[[sp]] <- read.csv(paste0("Data_Frames/",sp_clean_df_list[[sp]]))
    
    #load in genepop files into genind 
    sp_clean_temp_gen[[sp]] <- read.genepop(paste0("Adegenet_Files/",sp_clean_genepop_list[[sp]]), ncode = 3)
    #name individuals in genind 
    rownames(sp_clean_temp_gen[[sp]]@tab) <- sp_clean_temp_df[[sp]][,1]
    #name pops 
    levels(sp_clean_temp_gen[[sp]]@pop) <- unique(sp_clean_temp_df[[sp]][,2])
    
    ##Run relate red code 
    #Garden
    sp_garden_clean_temp_df <- sp_clean_temp_df[[sp]][,-2]
    
    #limit data frame by garden only 
    sp_garden_only_temp_df <- sp_garden_clean_temp_df[sp_garden_clean_temp_df[,2] == "Garden",]
    
    #run relatedness analysis 
    sp_garden_rel_df <- Demerelate(sp_garden_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #use function to limit individuals  
    sp_garden_halfsib_list <- half_sibling_reduction(sp_garden_rel_df)
    
    #reduce data frame by list of half-siblings 
    sp_garden_relate_red_df <- sp_garden_clean_temp_df[!sp_garden_clean_temp_df[,1] %in% sp_garden_halfsib_list,]
    
    #write garden df 
    sp_garden_relate_red_df <- sp_garden_relate_red_df[sp_garden_relate_red_df[,2] == "Garden",]
    
    #save in df 
    halfsib_df[sp, relate] <- paste0(signif(((length(sp_garden_only_temp_df[,1])-length(sp_garden_relate_red_df[,1]))/
                                               length(sp_garden_only_temp_df[,1])),3)*100, "% (", 
                                     length(sp_garden_relate_red_df[,1]),")")
    
    #write out 
    write.csv(sp_garden_relate_red_df, paste0("Data_Frames/Relate_Red/",  species_list[[sp]], 
                                              "_Garden_relate_red_", relatedness_analyses_list[[relate]], "_df.csv"),
              row.names = FALSE)
    
    ##Wild relatedness analysis 
    sp_wild_clean_temp_df <- sp_clean_temp_df[[sp]][sp_clean_temp_df[[sp]][,3] == "Wild",][,-3]
    
    #now run relatedness analysis 
    sp_wild_rel_df <- Demerelate(sp_wild_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #limit by analysis halfsibs 
    sp_wild_halfsib_list <- half_sibling_reduction(sp_wild_rel_df)
    
    #reduce data frame by list of half-siblings 
    sp_wild_relate_red_df <- sp_wild_clean_temp_df[!sp_wild_clean_temp_df[,1] %in% sp_wild_halfsib_list,]
    
    #number of ind removed 
    halfsib_df[sp, relate+3] <-paste0(signif(((length(sp_wild_clean_temp_df[,1])-length(sp_wild_relate_red_df[,1]))/
                                                length(sp_wild_clean_temp_df[,1])),3)*100, "% (", 
                                      length(sp_wild_relate_red_df[,1]), ")")
    
    #write csv of the reduced data frame 
    write.csv(sp_wild_relate_red_df, paste0("Data_Frames/Relate_Red/", species_list[[sp]], 
                                            "_Wild_relate_red_", relatedness_analyses_list[[relate]], "_df.csv"),
              row.names = FALSE)
    
  }
}

#write out summary table 
rownames(halfsib_df) <- species_list
colnames(halfsib_df) <- rep(relatedness_analyses_list, 2)
write.csv(halfsib_df, "../Analyses/Results/Sum_Stats/halfsib_df.csv")