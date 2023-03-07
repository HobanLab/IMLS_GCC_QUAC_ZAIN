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
species_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned", "ZAIN_red_sample")

#list scenarios 
scenario_list <- c("Garden_allSSR", "Wild_allSSR", "Garden_gSSR", "Wild_gSSR",
                   "Garden_EST", "Wild_EST")

#population lists for separating by garden/wild
#the first five are garden pops
#the last five are wild pops for both species 
pop_list <- list(c(1:17), c(1:17), c(1:10), c(1:10), c(1:10),
                 c(18:22), c(18:21), c(11:35), c(11:35), c(11:35))

ZAIN_garden_list <- list(c(1:10), c(1:10), c(1:10))
ZAIN_wild_red_list <- list( c(11:19, 23:26, 28:32, 34:35), 
                            c(11:19, 23:26, 28:32, 34:35),
                            c(11:19, 23:26, 28:32, 34:35))
                     

#QUAC loci lists - EST vs. gSSRs  
QUAC_EST_loci <- c("FIR031", "GOT009", "POR016", "FIR013", "FIR043", "GOTO40", 
                   "PIE039", "FIR53", "FIR048", "PIE125")
QUAC_gSSR_loci <- c("0C11", "1G13", "G07", "1F02","QpZAG9")

#allele frequency category lists 
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare",
                   "reg_rare","loc_com_d1","loc_com_d2","loc_rare")


#################################################
#     Comparing wild and garden populations     #
#################################################
#comparing wild and garden individuals 
allrich_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#comparing wild and garden individuals 
hexp_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#sum stats df for garden wild comp - allelic richness
allrich_df <- matrix(nrow = length(scenario_list), ncol = length(species_list)+3)

#sum stats df for garden wild comp - hexp
hexp_df <- matrix(nrow = length(scenario_list), ncol = length(species_list)+3)

#pvalue data frame 
sp_allrich_hexp_pvalue <- matrix(nrow = length(sp_genind_list)+3, ncol = 2)
colnames(sp_allrich_hexp_pvalue) <- c("All_Rich", "Hexp")
rownames(sp_allrich_hexp_pvalue) <- c(species_list, "ZAIN_og_wo_smallpops", 
                                      "ZAIN_rebinned_wo_smallpops", "ZAIN_sample_wo_smallpops")

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
  sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp+5]]])
  #rename to wild only 
  levels(sp_wild_genind@pop) <- rep("Wild", length(pop_list[[sp+5]]))
  
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
    
    #loop to create genind objects and diversity for genind objects 
    #without small ZAIN pops 
    for(pop in 1:length(ZAIN_garden_list)){
      
      #create a genind object without small pops for ZAIN, garden
      sp_genind_red_garden_temp <- repool(seppop(sp_genind_temp)[ZAIN_garden_list[[pop]]])
    
      #rename pops to just garden 
      levels(sp_genind_red_garden_temp@pop) <- rep("Garden", length(levels(sp_genind_red_garden_temp@pop)))
    
      #create a genind object without small pops for ZAIN, wild
      sp_genind_red_wild_temp <- repool(seppop(sp_genind_temp)[ZAIN_wild_red_list[[pop]]])
    
      #rename pops to just wild 
      levels(sp_genind_red_wild_temp@pop) <- rep("Wild", length(levels(sp_genind_red_wild_temp@pop)))
    
      #now repool 
      sp_genind_red_garden_wild <- repool(sp_genind_red_garden_temp, sp_genind_red_wild_temp)
    
      #calculate wild and garden allelic richness 
      allrich_df[1:2, sp+3] <- colMeans(allelic.richness(sp_genind_red_garden_wild)$Ar)
      
      #calculate hexp for garden and wild 
      hexp_df[1:2, sp+3] <- colMeans(as.data.frame(cbind(summary(seppop(sp_genind_red_garden_wild)[[1]])$Hexp,  
                                                     summary(seppop(sp_genind_red_garden_wild)[[2]])$Hexp)))
      
      #run hexp code
      sp_red_hexp <- as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp,  
                                     summary(seppop(sp_garden_wild_genind)[[2]])$Hexp))
      colnames(sp_red_hexp) <- c("Garden", "Wild")
      #create data frame for hexp
      sp_hexp_red_df <- gather(sp_red_hexp)
      
      #calculate p-values
      #for ZAIN, calculate allelic richness and hexp and compare
      sp_allrich_red_df <- gather(allelic.richness(sp_genind_red_garden_wild)$Ar)
      
      #run mann whitney u test for allelic richness 
      sp_allrich_hexp_pvalue[sp+3,1] <- kruskal.test(sp_allrich_red_df[,2]~sp_allrich_red_df[,1])[3]$p.value
      
     # run mann whitney u test for hexp
      sp_allrich_hexp_pvalue[sp+3,2] <- kruskal.test(sp_hexp_red_df[,2]~sp_hexp_red_df[,1])[3]$p.value
      
    }
  }
}

###write out summary tables for allelic richness and hexp comparisons
##label data frames
colnames(allrich_df) <- c(species_list, "ZAIN_og_wo_smallpops", 
                          "ZAIN_rebinned_wo_smallpops", "ZAIN_sample_wo_smallpops")
rownames(allrich_df) <- scenario_list
colnames(hexp_df) <- c(species_list, "ZAIN_og_wo_smallpops", 
                       "ZAIN_rebinned_wo_smallpops", "ZAIN_sample_wo_smallpops")
rownames(hexp_df) <- scenario_list

#write out data frames
write.csv(allrich_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_allrich_df.csv")
write.csv(hexp_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_hexp_df.csv")
write.csv(sp_allrich_hexp_pvalue, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_sp_allrich_hexp_pvalue_df.csv")

########################################
#     Allelic representation code      #
########################################
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

#without ZAIN small pops - create allele existing df
sp_all_red_exist_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
#without ZAIN small pops - create df of wild alleles captured by gardens
sp_wild_red_cap_df <- matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))
#without ZAIN small pops - data frame to record allele capture code
sp_allele_red_cap <-matrix(nrow = (length(dup_reps)), ncol = length(list_allele_cat))

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
      
    ##organize genind object
    #add individual names to each row of the tab 
    rownames(sp_genind_temp@tab) <- sp_df_temp[,1]
    #add pop names to the genind object 
    levels(sp_genind_temp@pop) <- unique(sp_df_temp$Pop)
    
    ##organize into pops - garden
    #separate into garden genind object 
    sp_garden_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp]]])
    #rename pops to be garden only 
    levels(sp_garden_genind@pop) <- rep("Garden", length(levels(sp_garden_genind@pop)))
    
    ##organize into pop types 
    #separate into wild genind object 
    sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[sp+5]]])
    #rename 
    levels(sp_wild_genind@pop) <- rep("Wild", length(levels(sp_wild_genind@pop)))
    
    #repool genind objects 
    sp_garden_wild_genind <- repool(sp_garden_genind, sp_wild_genind)
    
    #convert the wild genind object to a genpop object
    sp_wild_genpop <- genind2genpop(seppop(sp_garden_wild_genind)[2]$Wild)
      
    #create documents for comparison 
    n_ind_W <- nrow(sp_wild_genind@tab);  n_ind_G <- nrow(sp_garden_genind@tab); 
    sp_alleles_cap <- colSums(seppop(sp_garden_wild_genind)[[1]]@tab,na.rm=T)
      
    #first calculate the frequency categories of alleles in the wild individuals   	
    sp_allele_cat <- get.allele.cat(sp_wild_genpop, 1, 1, n_ind_W, n_drop = ndrop, glob_only = TRUE)	
    
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
  
  #add loop to calculate diversity in ZAIN without small pops 
  if(sp == 3|sp == 4|sp == 5){
    
    #loop to remove small pops from ZAIN and run diversity representation code
    for(pop in 1:length(ZAIN_wild_red_list)){
    
      ##organize into pop types 
      #separate into wild genind object 
      sp_wild_red_genind <- repool(seppop(sp_genind_temp)[ZAIN_wild_red_list[[pop]]])
      #rename 
      levels(sp_wild_red_genind@pop) <- rep("Wild", length(levels(sp_wild_red_genind@pop)))
    
      #repool genind objects 
      sp_garden_wild_red_genind <- repool(sp_garden_genind, sp_wild_red_genind)
    
      #convert the wild genind object to a genpop object
      sp_wild_red_genpop <- genind2genpop(seppop(sp_garden_wild_red_genind)[2]$Wild)
    
      #create documents for comparison 
      n_ind_W <- nrow(sp_wild_red_genind@tab);  n_ind_G <- nrow(sp_garden_genind@tab); 
      sp_alleles_red_cap <- colSums(seppop(sp_garden_wild_red_genind)[[1]]@tab,na.rm=T)
    
      #first calculate the frequency categories of alleles in the wild individuals   	
      sp_allele_red_cat <- get.allele.cat(sp_wild_red_genpop, 1, 1, n_ind_W, n_drop = ndrop, glob_only = TRUE)	
    
      #exterior loop to look at alleles by frequency category
      #interior loop to alleles by "duplication" amount - how many copies of each allele 
      for(cat in 1:length(list_allele_cat)){
        for(dup in 1:length(dup_reps)){
        
          #calculating alleles that exist by allelic category
          sp_all_red_exist_df[dup, cat] <- round(sum(sp_alleles_red_cap[sp_allele_red_cat[[cat]]] > dup_reps[[dup]]))
        
          #now determine how many wild alleles were captured per category 
          sp_wild_red_cap_df[dup, cat] <- round(sum(sp_alleles_red_cap[sp_allele_red_cat[[cat]]] > dup_reps[[dup]])/length(sp_allele_red_cat[[cat]]),4)
        
          #code to store as one data frame 
          sp_allele_red_cap[dup, cat] <- paste0(signif((sp_wild_red_cap_df[dup,cat]*100),3), "% (", signif(sp_all_red_exist_df[dup,cat],3), ")")
        
          
      }
    }
    
    }
    
    #without ZAIN small pops - alleles existing
    rownames(sp_all_red_exist_df) <- paste0(c(1:10), " or more copies")
    colnames(sp_all_red_exist_df) <- list_allele_cat
    #without ZAIN small pops - representing alleles
    rownames(sp_wild_red_cap_df) <- paste0(c(1:10), " or more copies")
    colnames(sp_wild_red_cap_df) <- list_allele_cat
    #without ZAIN small pops - comparing wild allele representation ex situ
    rownames(sp_allele_red_cap) <- paste0(c(1:10), " or more copies")
    colnames(sp_allele_red_cap) <- list_allele_cat
    
    write.csv(sp_all_red_exist_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_all_exist", n_drop_file, "_wo_smallpops.csv"))
    write.csv(sp_wild_red_cap_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_wildcap", n_drop_file, "_wo_smallpops.csv"))
    write.csv(sp_allele_red_cap, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_all_cap", n_drop_file, "_wo_smallpops.csv"))
    
    
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
  write.csv(sp_all_exist_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_all_exist", n_drop_file, ".csv"))
  write.csv(sp_wild_cap_df, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_wildcap", n_drop_file, ".csv"))
  write.csv(sp_allele_cap, paste0("../Analyses/Results/Garden_Wild_Comparison/",species_list[[sp]], "_all_cap", n_drop_file, ".csv"))
   
  
  
  }
}

#write session info out
sessionInfo()
