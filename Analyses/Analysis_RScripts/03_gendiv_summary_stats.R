###In this script, we calculated genetic diversity summary statistics 
##First, we tested linkage disequilibrium, null alleles, and Hardy Weinberg 
#equilibrium. Next, a data frame including expected heterozygosity, 
#allelic richness, number of alleles, mean longitude and latitude 
#for wild populations, and individual numbers. 
#The result tables are included in full in the supplemental text of this manuscript.
#We used files in this script that are referred to as "clean," which means
#clones and individuals with too much missing data have been removed. 

#########################
#        Libraries      #
#########################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(diveRsity)
library(parallel)
library(doParallel)

#########################
#   Load Data Files     #
#########################
#set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "Data_Frames/", pattern = "_clean_df.csv")

#list out all the final scenarios for analyzing 
scenario_list <- c("QUAC_wK", "QUAC_woK", 
                   "ZAIN_og_allpops", "ZAIN_rebinned_allpops", "ZAIN_sample")

############################################################
#  Null Alleles, HWE Deviation, Linkage Disequilibrium     #
############################################################
#loop to save genind objects
sp_genind_temp <- list()

#loop to load in genind objects and run HWE deviation test
for(sp in 1:length(scen_data_clean_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp[[sp]] <- read.genepop(paste0("Adegenet_Files/", sp_genind_list[[sp]]), ncode = 3)
  
  #load data frames 
  sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]]))
  
  #organize genind object
  rownames(sp_genind_temp[[sp]]@tab) <- sp_df_temp[,1]
  levels(sp_genind_temp[[sp]]) <- unique(sp_df_temp[,2])
  
  ##HWE devitations 
  #run HWE deviations by pop 
  sp_HWE_pop <- seppop(sp_genind_temp[[sp]]) %>% lapply(hw.test, B = 1000)
  
  #create df by pop for HWE devitations
  sp_HWE_pop_df <- sapply(sp_HWE_pop, "[", i = TRUE, j = 3)
  
  ##HWE deviation data frame prep
  #name columns with populations
  colnames(sp_HWE_pop_df) <- unique(sp_df_temp[,2])
  #round to the 3rd digit
  sp_HWE_allpop_df <- signif(sp_HWE_pop_df, 3)
  #write out HWE deviation data files
  write.csv(sp_HWE_allpop_df, paste0("../Analyses/Results/Sum_Stats/", scen_data_clean_list[[sp]], 
                                     "_HWE_dev_pop.csv"))
  
  #calculate % of null alleles/locus
  #use parallelization because this step takes quite a while without
  #Calculate the number of cores
  cores <- detectCores() - 1
  #Initiate cluster
  cl <- makeCluster(cores)
  
  #run null allele calculations over all genind objects
  sp_null_all <- parLapply(cl, sp_genind_temp, null.all)
  
  #create null allele frequency summary data frame
  sp_null_all_df <- signif(data.frame(sp_null_all[[sp]]$null.allele.freq$summary2),3)
  
  #stop clustering
  stopCluster(cl) 
  
  #write out to CSV
  write.csv(sp_null_all_df, paste0("../Analyses/Results/Sum_Stats/", scen_data_clean_list[[sp]] , 
                                   "_null_all_df.csv"))
  
  ##calculate linkage disequilibrium 
  #use parallelization because this step takes quite a while without
  #Calculate the number of cores
  cores <- detectCores() - 1
  #Initiate cluster
  cl <- makeCluster(cores)
  #calculate linkage disequilbrium
  sp_ld <- parLapply(cl, sp_genind_temp, pair.ia, sample = 1000)
  
  #convert to a data frame
  sp_ld_df <- data.frame(round(sp_ld[[sp]], digits = 2))
  
  #write out 
  write.csv(sp_ld_df, paste0("../Analyses/Results/Sum_Stats/", 
                             scen_data_clean_list[[sp]], "_LD.csv"))
  
  #stop clustering
  stopCluster(cl)
  
}

###########################################
#          Genetic Stats by Pop           #
###########################################

#loop to generate genetic summary statistics for populations 
for(sp in 1:length(sp_genind_list)){
  
    #load genepop files as genind objects 
    sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
    
    #load data frames 
    sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]]))
  
    #organize genind 
    levels(sp_genind_temp@pop) <- unique(sp_df_temp[,2])
    
    ##start genetic analyses
    #create genetic summary of the genind file 
    sp_sum <- summary(sp_genind_temp)
    #create poppr file 
    sp_poppr <- poppr(sp_genind_temp)
    #save mean for final output table 
    sp_hexp_mean <- sp_poppr[1:length(levels(sp_genind_temp@pop)),10]
    #allele numbers by pop 
    sp_nall <- sp_sum$pop.n.all
    #individual numbers
    sp_ind <- sp_poppr[1:length(levels(sp_genind_temp@pop)), 2:3]
    #save allelic richness for comparison
    sp_allrich_list <- allelic.richness(sp_genind_temp)$Ar
    sp_allrich_mean <- colMeans(allelic.richness(sp_genind_temp)$Ar)	
  
    #create data frame 
    sp_allpop_gendiv_sumstat_df <- signif(cbind(sp_ind, sp_nall, sp_allrich_mean, sp_hexp_mean),3)
    
    #name rows 
    rownames(sp_allpop_gendiv_sumstat_df) <- levels(sp_genind_temp@pop)
    colnames(sp_allpop_gendiv_sumstat_df) <- c("Ind","MLG", "NAll", "All_Rich", "Hexp")
    
    #write out data frame
    write.csv(sp_allpop_gendiv_sumstat_df, paste0("../Analyses/Results/Sum_Stats/", scenario_list[[sp]],  
                                                  "_gendiv_sumstats.csv"))
    
    #also run ZAIN individuals with reduced pop size  
    if(sp == 3|sp == 4|sp == 5){
      
      sp_genind_red_temp <- repool(seppop(sp_genind_temp)[c(1:19, 23:26, 28:32, 34:35)])
      
      sp_sum <- summary(sp_genind_red_temp)
      #create poppr file 
      sp_poppr <- poppr(sp_genind_red_temp)
      #save mean for final output table 
      sp_hexp_mean <- sp_poppr[1:length(levels(sp_genind_red_temp@pop)),10]
      #allele numbers by pop 
      sp_nall <- sp_sum$pop.n.all
      #individual numbers
      sp_ind <- sp_poppr[1:length(levels(sp_genind_red_temp@pop)), 2:3]
      #save allelic richness for comparison
      sp_allrich_list <- allelic.richness(sp_genind_red_temp)$Ar
      sp_allrich_mean <- colMeans(allelic.richness(sp_genind_red_temp)$Ar)	
      
      #create data frame 
      sp_redpop_gendiv_sumstat_df <- signif(cbind(sp_ind, sp_nall, sp_allrich_mean, sp_hexp_mean),3)
      
      #name rows 
      rownames(sp_redpop_gendiv_sumstat_df) <- levels(sp_genind_red_temp@pop)
      colnames(sp_redpop_gendiv_sumstat_df) <- c("Ind","MLG", "NAll", "All_Rich", "Hexp")
      
      #write out data frame
      write.csv(sp_redpop_gendiv_sumstat_df, paste0("../Analyses/Results/Sum_Stats/",
                                                    scenario_list[[sp]], "gendiv_sumstats_redpop.csv"))
      
    }
    
}
