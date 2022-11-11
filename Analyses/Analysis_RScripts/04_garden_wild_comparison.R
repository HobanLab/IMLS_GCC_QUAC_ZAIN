##This script details the analysis of actually comparing genetic 
#diversity levels between garden and wild Q. acerifolia and 
#Zamia integrifolia populations. 
#We first calculated diversity levels throughout all garden  
#and wild populations (indicated by allelic richness and 
#expected heterozygosity) and then ran a t-test to assess significance. 
#In this script we use "cleaned" genetic files because they have been
#cleaned for clones and individuals with too much missing data 
#(25% or more)

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)
library(tidyr)

#######################
#     Load files      #
#######################
#set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames", pattern = "_clean_df.csv")

#load in function to calculate allele frequency categories
source("../Analyses/Functions/Fa_sample_funcs.R")

#create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)

##write in if files need to be converted 
#sp_arp_list <- list.files(path = "Adegenet_Files/Garden_Wild", pattern = "clean.arp$")

#for(sp in 1:length(sp_arp_list)){

# arp2gen(paste0("Adegenet_Files/Garden_Wild/", sp_arp_list[[sp]]))

#}

#list out species
species_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#list scenarios 
scenario_list <- c("Garden_allSSR", "Wild_allSSR", "Garden_gSSR", "Wild_gSSR",
                   "Garden_EST", "Wild_EST")

#pop list 
pop_list <- list(c(1:17), c(1:17), c(1:10), c(1:10),
                 c(18:22), c(18:21), c(11:19, 23:26,28:32,34:35), 
                 c(11:19, 23:26,28:32,34:35))

#initial lists 
QUAC_EST_loci <- c("FIR031", "GOT009", "POR016", "FIR013", "FIR043", "GOTO40", 
                   "PIE039", "FIR53", "FIR048", "PIE125")
QUAC_gSSR_loci <- c("0C11", "1G13", "G07", "1F02","QpZAG9")

#################################################
#     Comparing wild and garden populations     #
#################################################
#comparing wild and garden individuals 
allrich_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#comparing wild and garden individuals 
hexp_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#sum stats df for garden wild comp - allelic richness
allrich_df <- matrix(nrow = length(scenario_list), ncol = length(species_list))

#sum stats df for garden wild comp - hexp
hexp_df <- matrix(nrow = length(scenario_list), ncol = length(species_list))

#pvalue data frame 
sp_allrich_hexp_pvalue <- matrix(nrow = length(sp_genind_list), ncol = 2)
colnames(sp_allrich_hexp_pvalue) <- c("All_Rich", "Hexp")
rownames(sp_allrich_hexp_pvalue) <- species_list

#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(sp_genind_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  ##organize into pops - garden
  #separate into garden genind object 
  sp_garden_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp]]])
  #rename pops to be garden only 
  levels(sp_garden_genind@pop) <- rep("Garden", length(pop_list[[sp]]))
  
  ##organize into pop types 
  #separate into wild genind object 
  sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp+4]]])
  #rename to wild only 
  levels(sp_wild_genind@pop) <- rep("Wild", length(pop_list[[sp+4]]))
  
  #repool to calculate diversity stats 
  sp_garden_wild_genind <- repool(sp_garden_genind, sp_wild_genind)
  
  #if statement for EST and gSSR comparison - QUAC only, sp == 1, sp == 2
  if(sp == 1|sp == 2){
    
    ###calculate diversity stats for all scenarios 
    ##just determine wild and garden diversity levels 
    #calculate allelic richness 
    allrich_df[1:2,sp] <- colMeans(allelic.richness(sp_garden_wild_genind)$Ar)
    #calculate hexp 
    hexp_df[1:2,sp] <- colMeans(as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp, 
                                                    summary(seppop(sp_garden_wild_genind)[[2]])$Hexp))) 
    
    #create statistical analysis data frame 
    sp_allrich_df <- gather(allelic.richness(sp_garden_wild_genind)$Ar)
    #calculate heterozygosity and create data frame 
    sp_hexp <- as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp, summary(seppop(sp_garden_wild_genind)[[2]])$Hexp))
    colnames(sp_hexp) <- c("Garden", "Wild")
    sp_hexp_df <- gather(sp_hexp)
    #name for category 
    sp_allrich_df$cat_type <- paste0(sp_allrich_df[,1],"_allrich")
    sp_hexp_df$cat_type <- paste0(sp_hexp_df[,1],"hexp")
    
    ##separate by loci combination 
    #gSSR genind object
    sp_gSSR_genind <- sp_garden_wild_genind[loc = QUAC_gSSR_loci]
    #calculate allelic richness for gSSRs
    allrich_df[3:4,sp] <- colMeans(allelic.richness(sp_gSSR_genind)$Ar)
    #calculate expected heterozygosity
    hexp_df[3:4,sp] <- colMeans(as.data.frame(cbind(summary(seppop(sp_gSSR_genind)[[1]])$Hexp, 
                                 summary(seppop(sp_gSSR_genind)[[2]])$Hexp)))
    
    #create statistical analysis data frame
    sp_gSSR_allrich <- gather(allelic.richness(sp_gSSR_genind)$Ar)
    #calculate hexp and create data frame 
    sp_gSSR_hexp <- as.data.frame(cbind(summary(seppop(sp_gSSR_genind)[[1]])$Hexp,  summary(seppop(sp_gSSR_genind)[[2]])$Hexp))
    colnames(sp_gSSR_hexp) <- c("Garden", "Wild")
    sp_gSSR_hexp <- gather(sp_gSSR_hexp)
    #create rownames with sections of names 
    sp_gSSR_allrich$cat_type <- paste0(sp_gSSR_allrich[,1], "_gSSR_allrich")
    sp_gSSR_hexp$cat_type <- paste0(sp_gSSR_allrich[,1], "_gSSR_hexp")
    
    ###calculate diversity stats for all scenarios 
    ##separate by loci combination 
    #EST SSR genind object 
    sp_EST_genind <- sp_garden_wild_genind[loc = QUAC_EST_loci]
    #calculate allelic richness and create data frame 
    sp_EST_allrich <- gather(allelic.richness(sp_EST_genind)$Ar)
    #calculate hexp and create data frame 
    sp_EST_hexp <- as.data.frame(cbind(summary(seppop(sp_EST_genind)[[1]])$Hexp,
                                       summary(seppop(sp_EST_genind)[[2]])$Hexp))
    colnames(sp_EST_hexp) <- c("Garden", "Wild")
    sp_EST_hexp <- gather(sp_EST_hexp)
    #create rownames with sections of names 
    sp_EST_allrich$cat_type <- paste0(sp_EST_allrich[,1], "_EST_allrich")
    sp_EST_hexp$cat_type <- paste0(sp_EST_hexp[,1], "_EST_hexp")
    
    #store allelic richness in data frame
    allrich_df[5:6,sp] <- colMeans(allelic.richness(sp_EST_genind)$Ar)
    #store expected heterozygosity in a data frame 
    hexp_df[5:6,sp] <- colMeans(as.data.frame(cbind(summary(seppop(sp_EST_genind)[[1]])$Hexp, 
                                                    summary(seppop(sp_EST_genind)[[2]])$Hexp)))
    
    ##combine all categories for statistical tests 
    #all rich
    sp_allrich_allcat_df <- rbind(sp_allrich_df, sp_gSSR_allrich, sp_EST_allrich)
    #hexp 
    sp_hexp_allcat_df <- rbind(sp_hexp_df, sp_gSSR_hexp, sp_EST_hexp)
    
    #summary data frames 
    sp_allrich_hexp_pvalue[sp,1] <- kruskal.test(sp_allrich_allcat_df[,2]~sp_allrich_allcat_df[,3])[3]$p.value
    sp_allrich_hexp_pvalue[sp,2] <- kruskal.test(sp_hexp_allcat_df[,2]~sp_hexp_allcat_df[,3])[3]$p.value
    
  }else{
    
    #calculate wild and garden allelic richness 
    allrich_df[1:2, sp] <- colMeans(allelic.richness(sp_garden_wild_genind)$Ar)
    #calculate hexp for garden and wild 
    hexp_df[1:2, sp] <- colMeans(as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp,  
                                                     summary(seppop(sp_garden_wild_genind)[[2]])$Hexp)))
    
    #for ZAIN, calculate allelic richness and hexp and compare
    sp_allrich_df <- gather(allelic.richness(sp_garden_wild_genind)$Ar)
    
    #run t-test 
    sp_allrich_hexp_pvalue[sp,1] <- kruskal.test(sp_allrich_df[,2]~sp_allrich_df[,1])[3]$p.value
    
    #run hexp code
    sp_hexp <- as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp,  
                                   summary(seppop(sp_garden_wild_genind)[[2]])$Hexp))
    colnames(sp_hexp) <- c("Garden", "Wild")
    #create data frame for hexp
    sp_hexp_df <- gather(sp_hexp)
    #save p-value for hexp 
    sp_allrich_hexp_pvalue[sp,2] <- kruskal.test(sp_hexp_df[,2]~sp_hexp_df[,1])[3]$p.value
    
    
  }
}

###write out summary tables for allelic richness and hexp comparisons
##label data frames
colnames(allrich_df) <- species_list
rownames(allrich_df) <- scenario_list
colnames(hexp_df) <- species_list
rownames(hexp_df) <- scenario_list

#write out data frames
write.csv(allrich_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_allrich_df.csv")
write.csv(hexp_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_hexp_df.csv")
write.csv(sp_allrich_hexp_pvalue, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_allrich_hexp_pvalue_df.csv")

########################################
#     Allelic representation code      #
########################################
#list out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

##create table for % alleles captured by frequency and how many duplicates were present  
#create list with duplicates 
dup_reps <- c(0:9)

#create a table to store % alleles captured by gardens pops where no alleles are dropped 
sp_allele_cap_table_ndrop0 <- matrix(nrow = length(dup_reps), ncol = length(list_allele_cat))

#create a table to store % alleles captured by garden pops where alleles are dropped if there are fewer than 2
sp_allele_cap_table_ndrop2 <- matrix(nrow = length(dup_reps), ncol = length(list_allele_cat))

#create arrays and lists to store results 
sp_allele_cat <- list()
#create allele existing df
sp_all_exist_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
#create df of wild alleles captured by gardens
sp_wild_cap_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
##data frame to record allele capture code
sp_allele_cap <-matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))

##run loop to generate allelic capture table 
#the outer loop is calculating how many copies of each allele in each category exists
#the inner loop is calculating the percent capture of each allele in each frequency category 
for(sp in 1:length(species_list)){  #loop over every scenario
  for(ndrop in c(0,2)){     #loop to include very rare or not 
    
    #ndrop or not   
    if(ndrop == 0) n_drop_file <- "_ndrop0"
    if(ndrop == 2) n_drop_file <- "_ndrop2"
      
    #load genepop files as genind objects 
    sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
      
    #load data frames 
    sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]])) 
      
    #organize genind object
    rownames(sp_genind_temp@tab) <- sp_df_temp[,1]
    
    ##organize into pops - garden
    #separate into garden genind object 
    sp_garden_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp]]])
    #rename pops to be garden only 
    levels(sp_garden_genind@pop) <- rep("Garden", length(pop_list[[sp]]))
    
    ##organize into pop types 
    #separate into wild genind object 
    sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp+4]]])
    #rename to wild only 
    levels(sp_wild_genind@pop) <- rep("Wild", length(pop_list[[sp+4]]))
    
    #repool garden and wild individuals 
    sp_garden_wild_genind <- repool(sp_garden_genind, sp_wild_genind)
    
    #calculate number of individuals per pop
    n_ind_p_pop <- as.numeric(table(sp_wild_genind@pop))
      
    #convert the wild genind object to a genpop object
    sp_wild_genpop <- genind2genpop(sp_wild_genind)
      
    #create documents for comparison 
    n_ind_W <- nrow(sp_wild_genpop@tab);  n_ind_G <- nrow(sp_garden_genind@tab); 
    sp_alleles_cap <- colSums(sp_garden_genind@tab,na.rm=T)
      
    #first calculate the frequency categories of alleles in the wild individuals   	
    sp_allele_cat <- get.allele.cat(sp_wild_genpop, 1, 1, as.numeric(n_ind_p_pop), n_drop = ndrop, glob_only = TRUE)	
      
    #exterior loop to look at alleles by frequency category
    #interior loop to alleles by "duplication" amount - how many copies of each allele 
    for(cat in 1:length(list_allele_cat)){
      for(dup in 1:length(dup_reps)){
        
      #calculating alleles that exist by allelic category
      sp_all_exist_df[dup, cat] <- round(sum(sp_alleles_cap[sp_allele_cat[[cat]]] > dup_reps[[dup]]))
    
      #now determine how many wild alleles were captured per category 
      sp_wild_cap_df[dup, cat] <- round(sum(sp_alleles_cap[sp_allele_cat[[cat]]] > dup_reps[[dup]])/length(sp_allele_cat[[cat]]),4)
      
      #code to store as one data frame 
      sp_allele_cap[dup, cat] <- paste0(signif((sp_wild_cap_df[dup,cat]*100),3), "% (", signif(sp_all_exist_df[dup,cat],3), ")")
     
       }
    }
  
  ##format tables
  #alleles existing
  rownames(sp_all_exist_df) <- paste0(c(1:10), " or more copies")
  colnames(sp_all_exist_df) <- list_allele_cat
  #percent capture of allele types by gardens
  rownames(sp_wild_cap_df) <- paste0(c(1:10), " or more copies")
  colnames(sp_wild_cap_df) <- list_allele_cat
  #comparison of percent of wild alleles captured in garden 
  rownames(sp_allele_cap) <- paste0(c(1:10), " or more copies")
  colnames(sp_allele_cap) <- list_allele_cat
  
  ##write out data frames
  write.csv(sp_all_exist_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], n_drop_file, ".csv"))
  write.csv(sp_wild_cap_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], n_drop_file, ".csv"))
  write.csv(sp_allele_cap, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], n_drop_file, ".csv"))
  
  }
}

#write session info out
sessionInfo()
