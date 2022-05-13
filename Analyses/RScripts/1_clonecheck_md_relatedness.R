##This script details the analyses run on Quercus acerifolia and Zamia integrifolia 
#(referred to as QUAC and ZAIN in this script for brevity) genotype files 
#to prepare them for genetic diversity and structure analyses. 
#When files are referred to as "clean" that means individuals 
#that are clones and individuals with too much missing data have been removed. 
#When files and objects are titled "red" that means they have been reduced
#for half-siblings (12.5% relatedness coefficient)

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

#create scenario list 
scenario_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#load in relatedness function 
source("../Analyses/RScripts/relatedness_analyses.R")

############################################################
#    Remove Clones and Individuals with Missing Data       #
############################################################

##loop to remove clones and individuals with missing data for genetic analyses
#loops over all scenarios for both species - QUAC and ZAIN 
for(sp in 1:length(sp_genind)){
  
  #load in genepop file as a genind object 
  sp_temp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind[[sp]]), ncode = 3)
  
  #load in genetic data frame  
  sp_temp_df <- read.csv(paste0("Data_Frames/",sp_df[[sp]]))
  
  #name rows in the genind object 
  rownames(sp_temp_genind@tab) <- sp_temp_df[,1]
  
  #create pop names list
  pop_names <- unique(sp_temp_df[,2])
  
  #name pops in genind object 
  levels(sp_temp_genind@pop) <- pop_names
  
  ##run clone check 
  #convert genind object to a genclone object 
  sp_temp_genclone <- as.genclone(sp_temp_genind)
  
  #identify multi-locus genotypes (non-clones)
  sp_mlg <- mlg.id(sp_temp_genclone)
  
  #create clone index 
  sp_clone_index <- which(sapply(sp_mlg, function(x) length(x)>1))
  
  #list clones
  sp_clone_id <- list()
  
  #create a list of the clone individuals 
  for(clones in 1:length(sp_clone_index)) sp_clone_id[[clones]] <- sp_mlg[[sp_clone_index[[clones]]]]
  
  #then create genind objects without clones 
  sp_genind_nocl <- clonecorrect(sp_temp_genind)
  
  #remove any individuals with more than 25% missing data 
  sp_genind_nomd <- missingno(sp_genind_nocl, type = "geno", cutoff = 0.25, quiet = FALSE, freq = FALSE) 
  
  #reduce data frame by nomd and clones
  sp_clean_df <- sp_temp_df[sp_temp_df[,1] %in% rownames(sp_genind_nomd@tab),]
  
  ##write out files cleaned for missing data and clones 
  #write out genalex
  genind2genalex(sp_genind_nomd, file = paste0("Data_Frames/", gsub("\\..*","",sp_genind[[sp]]), "_clean_genalex.csv"), 
              overwrite = TRUE)
  #write out data frame
  write.csv(sp_clean_df, paste0("Data_Frames/", gsub("\\..*","",sp_genind[[sp]]), "_clean_df.csv"), row.names = FALSE)


}

#################################
#      Relatedness Analysis     #
#################################
##code to reduce data frames by half-sibling relatedness 
#list clean data frames 
sp_clean_df_list <- list.files(path = "Data_Frames", pattern = "clean_df.csv")

#list clean genepop files 
sp_clean_genepop_list <- list.files(path = "Adegenet_Files", pattern = "allpop_clean.gen")

#create scenario list 
species_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#save summary df 
halfsib_df <- as.data.frame(matrix(nrow = length(species_list), ncol = 2))

#load in data frame 
QUAC_relate_wk_df <- read.csv("Data_Frames/QUAC_allpop_clean_df.csv")

#loop over species data frames 
for(sp in 1:length(species_list)){
  
  #load in garden data frame 
  sp_clean_temp_df <- read.csv(paste0("Data_Frames/",sp_clean_df_list[[sp]]))
  
  #load in genepop files into genind 
  sp_clean_temp_gen <- read.genepop(paste0("Adegenet_Files/",sp_clean_genepop_list[[sp]]), ncode = 3)
  #name individuals in genind 
  rownames(sp_clean_temp_gen@tab) <- sp_clean_temp_df[,1]
  #name pops 
  levels(sp_clean_temp_gen@pop) <- unique(sp_clean_temp_df[,2])
  
  ##Run relateness reduction code 
  #Garden
  sp_garden_clean_temp_df <- sp_clean_temp_df[,-2]
  
  #limit data frame by garden only 
  sp_garden_only_temp_df <- sp_garden_clean_temp_df[sp_garden_clean_temp_df[,2] == "Garden",]
  
  #limit genind object by garden 
  sp_garden_temp_genind <- sp_clean_temp_gen[rownames(sp_clean_temp_gen@tab) %in% sp_garden_only_temp_df[,1],]
  
  #run relatedness analysis 
  sp_garden_halfsib_list <- halfsib_relatedness_reduction_loiselle(sp_garden_only_temp_df)
  
  #reduce data frame by list of half-siblings 
  sp_garden_relate_red_df <- sp_garden_clean_temp_df[!sp_garden_clean_temp_df[,1] %in% sp_garden_halfsib_list,]
 
  #write garden df 
  sp_garden_relate_red_df <- sp_garden_relate_red_df[sp_garden_relate_red_df[,2] == "Garden",]
  
  #save in df 
  halfsib_df[sp,1] <- paste0(signif(((length(sp_garden_only_temp_df[,1])-length(sp_garden_relate_red_df[,1]))/
                                            length(sp_garden_only_temp_df[,1])),3)*100, "% (", 
                                   length(sp_garden_relate_red_df[,1]),")")
  
  #export data frame 
  write.csv(sp_garden_relate_red_df, paste0("Data_Frames/Relate_Red/", species_list[[sp]], 
                                            "_garden_relate_red_df.csv"), row.names = FALSE)
  
  #limit genind file  
  sp_relate_red_genind <- sp_garden_temp_genind[!rownames(sp_garden_temp_genind@tab) %in% sp_garden_halfsib_list,]
  
  #export genalex data frame 
  genind2genalex(sp_relate_red_genind, file = paste0("Data_Frames/Relate_Red/", species_list[[sp]], "_garden_relate_red_genalex.csv"), 
                 overwrite = TRUE)
  
  ##Wild relatedness analysis 
  sp_wild_clean_temp_df <- sp_clean_temp_df[sp_clean_temp_df[,3] == "Wild",][,-3]
  
  #limit genind object by wild individuals 
  sp_wild_genind <- sp_clean_temp_gen[rownames(sp_clean_temp_gen@tab) %in% sp_wild_clean_temp_df[,1]]

  #limit by analysis halfsibs 
  sp_wild_halfsib_list <- halfsib_relatedness_reduction_loiselle(sp_wild_clean_temp_df)
  
  #reduce data frame by list of half-siblings 
  sp_wild_relate_red_df <- sp_wild_clean_temp_df[!sp_wild_clean_temp_df[,1] %in% sp_wild_halfsib_list,]
  
  #number of ind removed 
  halfsib_df[sp, 2] <-paste0(signif(((length(sp_wild_clean_temp_df[,1])-length(sp_wild_relate_red_df[,1]))/
                                              length(sp_wild_clean_temp_df[,1])),3)*100, "% (", 
                                              length(sp_wild_relate_red_df[,1]), ")")
  #export data files 
  write.csv(sp_wild_relate_red_df, paste0("Data_Frames/Relate_Red/", species_list[[sp]], "_wild_relate_red_df.csv"))

  #limit genind object by half-siblings
  sp_wild_relate_red_genind <- sp_wild_genind[!rownames(sp_wild_genind@tab) %in% sp_wild_halfsib_list,]
  
  #export genalex data frame 
  genind2genalex(sp_wild_relate_red_genind, file = paste0("Data_Frames/Relate_Red/", species_list[[sp]], "_wild_relate_red_genalex.csv"), 
                 overwrite = TRUE)
  
}

#write out summary table 
rownames(halfsib_df) <- species_list
colnames(halfsib_df) <- c("Garden", "Wild")
write.csv(halfsib_df, "../Analyses/Results/Sum_Stats/halfsib_sum_df.csv")
