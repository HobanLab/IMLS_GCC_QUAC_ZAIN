##This script details the analyses run on Quercus acerifolia (referred to as 
#QUAC for brevity here) genotype files to prepare them for genetic diversity and 
#structure analyses. When files are referred to as "clean" that means individuals 
#that are clones and individuals with too much missing data have been removed. 
#When files and objects are titled "red" that means they have been reduced
#for relatedness (25% or more related individuals are reduced to one individual
#per phenotype)

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
setwd("G:/Shared drives/Emily_Schumacher/GCC_QUAC_ZAIN/Data_Files")

##convert QUAC all pop arlequin file into genepop file if needed
#arp2gen("C:\\Users\\eschumacher\\Documents\\GitHub\\QUAC_insitu_exsitu\\QUAC_data_files\\QUAC_adegenet_files\\QUAC_allpop.arp")

#list species 
species_list <- c("QUAC", "ZAIN")

#now read in genepop file as a genind for adegenet 
sp_genind <- list.files(path = "Adegenet_Files", pattern = "allpop.gen$")

#load relatedness data frame for relatedness analysis 
sp_df <- list.files(path = "Data_Frames", pattern = "allpop_df.csv$")

#load in relatedness function 
source("../Analyses/RScripts/relatedness_analyses.R")

############################################################
#    Remove Clones and Individuals with Missing Data       #
############################################################
setwd("G:/Shared drives/Emily_Schumacher/GCC_QUAC_ZAIN/Data_Files")

##loop over both QUAC and ZAIN to generate results  
for(sp in 1:length(species_list)){
  
  #load in genepop files as a genind object 
  sp_temp_genind <- read.genepop(paste0("Adegenet_Files/",sp_genind[[sp]]), ncode = 3)
  
  #load in data frames 
  sp_temp_df <- read.csv(paste0("Data_Frames/",sp_df[[sp]]))
  
  #name rows in the genind object 
  rownames(sp_temp_genind@tab) <- sp_temp_df[,1]
  
  #create pop names list
  pop_names <- unique(sp_temp_df[,2])
  
  #name pops in genind object 
  levels(sp_temp_genind@pop) <- pop_names
  
  ##run clone check 
  #convert to a genclone object 
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
  #genind2genalex(sp_genind_nomd, file = paste0("Data_Frames/", species_list[[sp]], "_clean_genalex.csv"), 
  #            overwrite = TRUE)
  #write out data frame
  #write.csv(sp_clean_df, paste0("Data_Frames/", species_list[[sp]], "_clean_df.csv"), row.names = FALSE)


}

#################################
#      Relatedness Analysis     #
#################################
setwd("G:/Shared drives/Emily_Schumacher/IMLS_QUAC_ZAIN/Data_Files")

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
