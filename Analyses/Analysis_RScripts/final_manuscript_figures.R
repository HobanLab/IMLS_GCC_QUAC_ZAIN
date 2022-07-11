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
#setwd("../../Data_Files")
setwd("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files/Garden_Wild/", pattern = "_clean.gen")

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

#list of scenarios 
species_list <- c("QUAC_wK", "QUAC_wK_ESTSSR", "QUAC_wK_garden_wild_gSSR",
                  "QUAC_woK", "QUAC_woK_ESTSSR", "QUAC_woK_garden_wild_gSSR",
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

#create lists 
allrich_df_list <- list()
hexp_df_list <- list()

#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(species_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/Garden_Wild/",sp_genind_list[[sp]]), ncode = 3)
  
  #name populations 
  levels(sp_genind_temp@pop) <- c(paste0(species_list[[sp]], "_Garden"), paste0(species_list[[sp]], "_Wild"))
  
  #combining into a df 
  allrich_df_list[[sp]] <- gather(allelic.richness(sp_genind_temp)$Ar)
  
  
  
  #run t-test 
  #allrich_pvalue <- as.numeric(kruskal.test(allrich_df_list[[sp]][,2]~allrich_df_list[[sp]][,1])[3])
  
  #name data frame
  # allrich_hexp_df[1:2,sp] <- as.numeric(colMeans(as.data.frame(allelic.richness(sp_genind_temp)$Ar)))
  #allrich_hexp_df[3,sp] <- allrich_pvalue
  
  #run hexp code
  # hexp_df <- cbind(as.numeric(summary(seppop(sp_genind_temp)[[1]])$Hexp), as.numeric(summary(seppop(sp_genind_temp)[[2]])$Hexp))
  #name columns 
  # colnames(hexp_df) <- c("Garden", "Wild")
  
  #transform the data frame for analyses 
  # hexp_temp_df[[sp]] <- gather(as.data.frame(hexp_df))
  
  #save p-value 
  # hexp_pvalue <- as.numeric(kruskal.test(hexp_temp_df[,2]~hexp_temp_df[,1])[3])
  
  #save in df 
  # allrich_hexp_df[(1:2)+3,sp] <- as.numeric(colMeans(hexp_df))
  # allrich_hexp_df[6,sp] <- hexp_pvalue
  
  #name colnames and rownames 
  #colnames(allrich_hexp_df) <- species_list
  #rownames(allrich_hexp_df) <- c("allrich_garden", "allrich_wild", "allrich_pvalue",
  #                               "hexp_garden", "hexp_wild", "hexp_pvalue")
}

#create a data frame of all QUAC allelic richness results 
QUAC_allrich_woK_df <- rbind(allrich_df_list[[4]], allrich_df_list[[5]], allrich_df_list[[6]])



#create a boxplot of only with Kessler individuals table
boxplot(QUAC_allrich_df[61:120,2]~QUAC_allrich_df[61:120,1])

#run kruskal-wallis 
QUAC_aov_test <- aov(QUAC_allrich_df[,2]~QUAC_allrich_df[,1])

#create data frame of allrich for ZAIN 
ZAIN_allrich_df <- rbind(allrich_df_list[[7]], allrich_df_list[[8]])

ZAIN_kw_test <- kruskal.test(ZAIN_allrich_df[,2]~ZAIN_allrich_df[,1])

##create a boxplot of just rebinned and woK QUAC 
#first, create a data frame with only the cases we care about 
allrich_final_df <- rbind(QUAC_allrich_df[61:120,], ZAIN_allrich_df[23:44,])

allrich_final_df <- allrich_final_df %>%
                      mutate(Treatment = gsub("^.*_","",allrich_final_df$key))

#separate by species 
allrich_final_df <- allrich_final_df %>%
                      mutate(Species = gsub('_.*','',allrich_final_df$Pop_Type))

colnames(allrich_final_df) <- c("Pop_Type", "Allelic_Richness", "Treatment")


ggplot(allrich_final_df, aes(x=Pop_Type, y=Allelic_Richness, fill=Treatment)) + 
  geom_boxplot() + xlab("Population Type") + ylab("Allelic Richness") + ylim(0,20) + 
  theme_bw() + facet_wrap(~Species, scale="free") + 
  scale_x_discrete(labels=c("Garden","Wild", "gSSR_Garden", "EST_Garden","gSSR_Wild", "EST_Wild")) + 
  scale_fill_manual(values = c("darkseagreen1", "darkgreen"))

#write out df 
write.csv(allrich_hexp_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_garden_wild_df.csv")


##new try 
ZAIN_allrich_df <- ZAIN_allrich_df %>%
                    mutate(Pop_Type = gsub("^.*_","",ZAIN_allrich_df$key))

allrich_df <- rbind(QUAC_allrich_woK_df, ZAIN_allrich_df[23:44,])

#add a species column 
allrich_df <- allrich_df %>%
                  mutate(Species = gsub('_.*','',allrich_df$key))


#name columns of the data frame 
colnames(allrich_df) <- c("Pop", "All_Rich", "Pop_Type", "Species")



#rename 
for(a in 1:length(allrich_df[,1])){
  
  tf <- allrich_df[a,1] == "QUAC_woK_Garden"
  tf2 <- allrich_df[a,1] == "QUAC_woK_Wild"
  
  if(tf == TRUE){
  
  
    allrich_df[a,1] <- "QUAC_A_woK_Garden"
  
  }
  if(tf2 == TRUE){
    
    allrich_df[a,1] <- "QUAC_A_woK_Wild"
    
  }
  
}

#ok let's try this again 
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_barplot.pdf", width = 12, height = 8)
ggplot(allrich_df, aes(x = Pop, y = All_Rich, fill = Pop_Type)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("darkseagreen1", "darkseagreen"),
                                     labels = expression("Garden", "Wild")) +
  scale_color_manual(labels=c("Garden","Wild",
                              "EST Garden","EST Wild",
                              "gSSR Garden", "gSSR Wild",
                              "Garden", "Wild")) +
  ylim(0,20) +
  xlab("Population Type") + ylab("Allelic Richness") + 
  facet_wrap(~Species, scale="free") + theme_bw() 
  dev.off()
  

#This function is run to output the resampling graphs for 
resampling_df <- resampling[-1,c(2:6)]
    
    #write PDF with name
    pdf(paste0(sp,"_resample_plot.pdf"))
    #add points
    plot(resampling_df[,1], col = "red", pch = 20, xlab = "Number of Individuals", 
         ylab = "Percent Diversity Capture", xlim = c(0,length(rownames(resampling))), ylim = c(0,100), cex = 1.2,
         main = "Percent Diversity Capture (All Alleles Included)")
    points(resampling_df[,2], col = "firebrick", pch = 20, cex = 1.2)
    points(resampling_df[,3], col = "darkorange3", pch = 20, cex = 1.2)
    points(resampling_df[,4], col = "coral", pch = 20, cex = 1.2)
    points(resampling_df[,5], col = "deeppink4", pch = 20, cex = 1.2)
    
    legend('bottomright', legend = c("Global", "Very Common", "Common", "Uncommon","Rare"),
           col = c("red", "firebrick", "darkorange3", "coral", "deeppink4"), pch = "_")

    
dev.off()
    

  


#write session info out
sessionInfo()
