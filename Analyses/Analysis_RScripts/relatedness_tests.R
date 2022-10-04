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

#load relatedness data frame for relatedness analysis 
sp_df <- list.files(path = "Data_Frames", pattern = "allpop_df.csv$")

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

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

###################################
#     Demerelate Ref Pop Test     #
###################################

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


##write out data table from the different pop structures analysis 
#name rows and columns 
rownames(rel_popstr_table) <- rel_scen
colnames(rel_popstr_table) <- c("non_halfsibs")
#write out table
write.csv(rel_popstr_table, "../Analyses/Results/Relatedness/rel_popstr_table.csv")

###################################
#     Relate Indicator Tests      #
###################################
##This area details the part of the script testing how if relatedness indicators
#identified different numbers of siblings within each pop type (garden or wild)

#list clean data frames 
sp_clean_df_list <- list.files(path = "Data_Frames", pattern = "clean_df.csv")

#relatedness lists for all analyses
sp_garden_relate_levels <- list()
sp_wild_relate_levels <- list()

#create a list of the relatedness tests 
relatedness_analyses_list <- c("loiselle", "wang", "ritland")

#save summary df 
relate_ind_fullsib_df <- matrix(nrow = length(scenario_list), ncol = 
                                    length(relatedness_analyses_list)*2)


###this analysis compares the number of full siblings determined by different
##relatedness indicators 
#this loop runs over the main data files 
for(sp in 1:length(scenario_list)){
  #this loop runs over different relatedness analyses
  for(relate in 1:length(relatedness_analyses_list)){
    #load in clean data frames to perform relatedness analysis on 
    sp_clean_temp_df <- read.csv(paste0("Data_Frames/",sp_clean_df_list[[sp]]))
    
    ##Run relate red code 
    #Garden
    sp_garden_clean_temp_df <- sp_clean_temp_df[,-2]
    
    #limit data frame by garden only 
    sp_garden_only_temp_df <- sp_garden_clean_temp_df[sp_garden_clean_temp_df[,2] == "Garden",]
    
    #run relatedness analysis 
    sp_garden_rel_df <- Demerelate(sp_garden_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #save relatedness levels for each scenario to graph  
    sp_relate_levels[[sp]] <- sp_garden_rel_df$Empirical_Relatedness$Garden 
    
    #species full siblings 
    sp_garden_fullsib_list <- fullsib_loiselle_rel_fun(names(which(unlist(sp_garden_rel_df$Empirical_Relatedness$Garden) > 0.25)))
    
    #save in df 
    relate_ind_fullsib_df[sp, relate] <- paste0(signif(((length(sp_garden_fullsib_list))/length(sp_garden_only_temp_df[,1])),3)*100, "% (", 
                                                length(sp_garden_only_temp_df[,1]),")")
    
    
    ##Wild relatedness analysis 
    sp_wild_clean_temp_df <- sp_clean_temp_df[sp_clean_temp_df[,3] == "Wild",][,-3]
    
    #now run relatedness analysis 
    sp_wild_rel_df <- Demerelate(sp_wild_clean_temp_df, object = T, value = relatedness_analyses_list[[relate]])
    
    #limit by analysis halfsibs 
    sp_wild_fullsib_list <- fullsib_loiselle_rel_fun(names(which(unlist(sp_wild_rel_df$Empirical_Relatedness) > 0.25)))
     
    #number of ind removed 
    relate_ind_fullsib_df[sp, relate+3] <- paste0(signif(length(sp_wild_fullsib_list)/
                                                        length(sp_wild_clean_temp_df[,1])*100, 3), "% ",
                                                  "(", length(sp_wild_clean_temp_df[,1]) - length(sp_wild_fullsib_list), ")")
    
  }
}

#write out summary table 
rownames(relate_ind_fullsib_df) <- scenario_list
colnames(relate_ind_fullsib_df) <- rep(relatedness_analyses_list, 2)
write.csv(relate_ind_fullsib_df, "../Analyses/Results/Relatedness/relate_ind_fullsib_df.csv")


##making histograms of relatedness
#load in QUAC with Kessler data frame
QUAC_rel_df <- read.csv("Data_Frames/QUAC_allpop_clean_df.csv")

#list to save results 

#loop to run relatedness analyses 
for(r in 1:length(relatedness_analyses_list)){
  
  ##Garden
  #create garden df
  sp_garden_temp_df <- QUAC_rel_df[,-2]
  
  #run relatedness analysis
  sp_garden_rel_df <- Demerelate(sp_garden_temp_df, object = T, value = relatedness_analyses_list[[r]])
  
  #generate histogram 
  pdf(paste0("../Analyses/Results/Relatedness/QUAC_garden_", relatedness_analyses_list[[r]],
             "_distribution.pdf"), width = 8, height = 8)
  hist(sp_garden_rel_df$Empirical_Relatedness$Garden, xlim = c(-1,1), ylim = c(0, 12000), 
       col = "darkseagreen1", xlab = paste0(relatedness_analyses_list[[r]]," Garden"),
       main = paste0("QUAC garden ",relatedness_analyses_list[[r]], " Distribution"))
  dev.off()
  
  ##Wild
  #limit to wild pops and remove pop type column 
  sp_wild_temp_df <- QUAC_rel_df[QUAC_rel_df$Garden_Wild == "Wild",-3]
  
  #run relatedness analysis
  sp_wild_rel_df <- Demerelate(sp_wild_temp_df, object = T, value = relatedness_analyses_list[[1]])
  
  #generate histogram 
  pdf(paste0("../Analyses/Results/Relatedness/QUAC_wild_", relatedness_analyses_list[[r]],
             "_distribution.pdf"), width = 8, height = 8)
  hist(as.numeric(unlist(sp_wild_rel_df$Empirical_Relatedness)), 
       xlim = c(-1,1), ylim = c(0, 1000), 
       col = "forestgreen", xlab = paste0(relatedness_analyses_list[[r]]," Wild"),
       main = paste0("QUAC Wild ", relatedness_analyses_list[[r]], " Distribution"))
  dev.off()
  
}
