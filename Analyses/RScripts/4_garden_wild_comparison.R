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
##set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/Garden_Wild/", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames", pattern = "_clean_df.csv")

#load in function to calculate allele frequency categories
source("../Analyses/RScripts/Fa_sample_funcs.R")

#create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)

##write in if files need to be converted 
#sp_arp_list <- list.files(path = "Adegenet_Files/Garden_Wild", pattern = "clean.arp$")

#for(sp in 1:length(sp_arp_list)){

# arp2gen(paste0("Adegenet_Files/Garden_Wild/", sp_arp_list[[sp]]))

#}

#list of scenarios 
species_list <- c("QUAC_wK", "QUAC_wk_ESTSSR", "QUAC_wk_gSSR",
                  "QUAC_woK", "QUAC_woK_ESTSSR", "QUAC_woK_gSSR",
                  "ZAIN_og", "ZAIN_rebinned")

#list out the scenarios 
scenarios_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#################################################
#     Comparing wild and garden populations     #
#################################################
#comparing wild and garden individuals 
allrich_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#comparing wild and garden individuals 
hexp_garden_wild_df <- as.data.frame(matrix(nrow = 3, ncol = length(species_list)))

#sum stats df for garden wild comp 
allrich_hexp_df <- matrix(nrow = 6, ncol = length(species_list))

#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(species_list)){

  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/Garden_Wild/",sp_genind_list[[sp]]), ncode = 3)
  
  #load data frames 
  sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]])) 
  
  #organize genind object
  rownames(sp_genind_temp@tab) <- sp_df_temp[,1]
  levels(sp_genind_temp@pop) <- c("Garden", "Wild")
    
  #combining into a df 
  allrich_df <- gather(allelic.richness(sp_genind_temp)$Ar)
  
  #run t-test 
  allrich_pvalue <- as.numeric(kruskal.test(allrich_df[,2]~allrich_df[,1])[3])
 
  #name data frame
  allrich_hexp_df[1:2,sp] <- as.numeric(colMeans(as.data.frame(allelic.richness(sp_genind_temp)$Ar)))
  allrich_hexp_df[3,sp] <- allrich_pvalue
  
  #run hexp code
  hexp_df <- cbind(as.numeric(summary(seppop(sp_genind_temp)[[1]])$Hexp), as.numeric(summary(seppop(sp_genind_temp)[[2]])$Hexp))
  #name columns 
  colnames(hexp_df) <- c("Garden", "Wild")
  
  #transform the data frame for analyses 
  hexp_temp_df <- gather(as.data.frame(hexp_df))
  
  #save p-value 
  hexp_pvalue <- as.numeric(kruskal.test(hexp_temp_df[,2]~hexp_temp_df[,1])[3])
  
  #save in df 
  allrich_hexp_df[(1:2)+3,sp] <- as.numeric(colMeans(hexp_df))
  allrich_hexp_df[6,sp] <- hexp_pvalue
 
  #name colnames and rownames 
  colnames(allrich_hexp_df) <- species_list
  rownames(allrich_hexp_df) <- c("allrich_garden", "allrich_wild", "allrich_pvalue",
                                 "hexp_garden", "hexp_wild", "hexp_pvalue")
  
}


#write out df 
write.csv(allrich_hexp_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_garden_wild_df.csv")


#################################
#     Allelic capture code      #
#################################
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
for(sp in 1:length(scenarios_list)){  #loop over every scenario
  for(ndrop in c(0,2)){     #loop to include very rare or not 
    
    #ndrop or not   
    if(ndrop == 0) n_drop_file <- "_ndrop0"
    if(ndrop == 2) n_drop_file <- "_ndrop2"
      
    #load genepop files as genind objects 
    sp_genind_temp <- read.genepop(paste0("Adegenet_Files/Garden_Wild/",sp_genind_list[[sp]]), ncode = 3)
      
    #load data frames 
    sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]])) 
      
    #organize genind object
    rownames(sp_genind_temp@tab) <- sp_df_temp[,1]
    levels(sp_genind_temp@pop) <- unique(sp_df_temp[,3])
      
    #seppop genind - we only want to use the wild pops to calculate the alleles existing for capture
    sp_seppop <- seppop(sp_genind_temp)
      
    #calculate number of individuals per pop
    n_ind_p_pop <- as.numeric(table(sp_seppop[[2]]@pop))
      
    #convert the wild genind object to a genpop object
    sp_wild_genpop <- genind2genpop(sp_seppop[[2]])
      
    #create documents for comparison 
    n_ind_W<-table(sp_genind_temp@pop)[2];  n_ind_G<-table(sp_genind_temp@pop)[1]; 
    sp_alleles_cap <- colSums(sp_seppop[[1]]@tab,na.rm=T)
      
    #first calculate the frequency categories of alleles in the wild individuals   	
    sp_allele_cat <- get.allele.cat(sp_wild_genpop, 1, 1, as.numeric(n_ind_p_pop), n_drop = ndrop, glob_only = TRUE)	
      
    #create a data frame with all of the alleles existing by category
    for(cat in 1:length(list_allele_cat)){
      for(dup in 1:length(dup_reps)){
        
      #calculating 
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
  
  #write out duplicate code graph 
  #pdf(paste0("QUAC_wild_cap_dup_barplot",n_drop_file, '.pdf'), width = 8, height = 8)
  #barplot(t(QUAC_wild_cap_df[,c(1, 3, 5)])*100, las = 2, beside = TRUE,
  #        col = c("red", "coral", "deeppink4"), legend.text = c("Global Alleles", "Common Alleles", "Rare Alleles"))
  #dev.off()
  }
}

#####################
#     Plotting      #
#####################
#create mean data frames - allrich 
QUAC_allrich_mean_df <- as.data.frame(rbind(mean(QUAC_allrich[[1]][,1]), mean(QUAC_allrich[[2]][,1])))
QUAC_allrich_mean_df$pop_type <- NA
QUAC_allrich_mean_df$pop_type <- c("Garden", "Wild")

#calculate standard errors
allrich_garden_se <- sd(QUAC_allrich[[1]][,1])/sqrt(length(QUAC_allrich[[1]][,1]))
allrich_wild_se <- sd(QUAC_allrich[[2]][,1])/sqrt(length(QUAC_allrich[[2]][,1]))

#create mean data frame - hexp
QUAC_hexp_mean_df <- as.data.frame(rbind(mean(QUAC_hexp[[1]][,1]), mean(QUAC_hexp[[2]][,1])))
QUAC_hexp_mean_df$pop_type <- NA
QUAC_hexp_mean_df$pop_type <- c("Garden", "Wild")

#calculate standard errors
hexp_garden_se <- sd(QUAC_hexp[[1]][,1])/sqrt(length(QUAC_hexp[[1]][,1]))
hexp_wild_se <- sd(QUAC_hexp[[2]][,1])/sqrt(length(QUAC_hexp[[2]][,1]))


#allrich comparison boxplot
pdf("allrich_garden_wild_barplot.pdf", width = 8, height = 10)
barplot(QUAC_allrich_mean_df[,1], beside = TRUE, 
        ylim = c(0,15), col = c("darkgreen", "darkseagreen1"),
        names = c("Garden", "Wild"), 
        main = "Allelic Richness Compared Between Garden and Wild Populations", 
        xlab = "Population Type", ylab = "Allelic Richness")
arrows(x0 = 0.7, y0 = QUAC_allrich_mean_df[1,1] - allrich_garden_se, 
       x1 = 0.7, y1 = QUAC_allrich_mean_df[1,1] + allrich_garden_se,
       code=3, angle=90, length=0.1)

arrows(x0 = 1.9, y0 = QUAC_allrich_mean_df[2,1] - allrich_wild_se, 
       x1 = 1.9, y1 = QUAC_allrich_mean_df[2,1] + allrich_wild_se,
       code=3, angle=90, length=0.1)

abline(h = 0)
dev.off()

#hexp barplot

#hexp comparison boxplot
pdf("hexp_garden_wild_barplot.pdf", width = 8, height = 10)
barplot(QUAC_hexp_mean_df[,1], beside = TRUE, 
        ylim = c(0,1), col = c("darkgreen", "darkseagreen1"),
        names = c("Garden", "Wild"), 
        main = "Expected Heterozygosity Compared Between Garden and Wild Populations", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
arrows(x0 = 0.7, y0 = QUAC_hexp_mean_df[1,1] - hexp_garden_se, 
       x1 = 0.7, y1 = QUAC_hexp_mean_df[1,1] + hexp_garden_se,
       code=3, angle=90, length=0.1)

arrows(x0 = 1.9, y0 = QUAC_hexp_mean_df[2,1] - hexp_wild_se, 
       x1 = 1.9, y1 = QUAC_hexp_mean_df[2,1] + hexp_wild_se,
       code=3, angle=90, length=0.1)

abline(h = 0)
dev.off()


#write session info out
sessionInfo()
