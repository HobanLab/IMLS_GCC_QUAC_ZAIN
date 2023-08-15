###This script contains most of the work to generate final manuscript 

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)
library(poppr)
library(hierfstat)
library(tidyr)
library(dplyr)

#######################
#     Load files      #
#######################
#set working directory to load in data files
setwd("../../Data_Files")

#genind objects 
sp_genind_list <- list.files(path = "Adegenet_Files", pattern = "_clean.gen")

#df files 
sp_df_list <- list.files(path = "CSV_Files", pattern = "_clean_df.csv")

#load in function to calculate allele frequency categories
source("../Analyses/Functions/Fa_sample_funcs.R")

#create functions to run code 
colMax <- function(data) sapply(data, max, na.rm = TRUE)

##write in if files need to be converted 
#sp_arp_list <- list.files(path = "Adegenet_Files/Garden_Wild", pattern = "clean.arp$")

#for(sp in 1:length(sp_arp_list)){

# arp2gen(paste0("Adegenet_Files/Garden_Wild/", sp_arp_list[[sp]]))

#}

#list out the scenarios 
scenarios_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned", 
                    "ZAIN_rebinned_sample")

#species list 
species_list <- c("QUAC", "ZAIN")

#QUAC loci lists  
QUAC_EST_loci <- c("FIR031", "GOT009", "POR016", "FIR013", "FIR043", "GOTO40", 
                   "PIE039", "FIR53", "FIR048", "PIE125")
QUAC_gSSR_loci <- c("0C11", "1G13", "G07", "1F02","QpZAG9")

#####################
#     Figure 2      #
#####################
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
for(sp in 1:length(scenarios_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  #load data frame 
  sp_df_temp <- read.csv(paste0("CSV_Files/", sp_df_list[[sp]]))
  
  #reorganize genind object for QUAC without Kessler for figure 
  if(sp == 2){
    
    ##reorganize individuals into garden and wild pools, but also by loci type 
    #create QUAC garden/wild genind object 
    #first, separate garden genind object 
    QUAC_garden_genind <- repool(seppop(sp_genind_temp)[1:17])
    levels(QUAC_garden_genind@pop) <- rep("Garden", 17)
    
    #create wild pop genind 
    QUAC_wild_genind <- repool(seppop(sp_genind_temp)[18:21])
    levels(QUAC_wild_genind@pop) <- rep("Wild", 4)
    
    #repool wild/garden genind object 
    QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)
    
    ###calculate diversity stats for all scenarios 
    #create statistical analysis data frame 
    QUAC_allrich_df <- gather(allelic.richness(QUAC_garden_wild_genind)$Ar)
    #calculate heterozygosity and create data frame 
    QUAC_hexp <- as.data.frame(cbind(summary(seppop(QUAC_garden_wild_genind)[[1]])$Hexp, 
                                      summary(seppop(QUAC_garden_wild_genind)[[2]])$Hexp))
    colnames(QUAC_hexp) <- c("Garden", "Wild")
    QUAC_hexp_df <- gather(QUAC_hexp)
    
    #name for category 
    QUAC_allrich_df$cat_type <- paste0("QUAC_",QUAC_allrich_df[,1],"allrich")
    QUAC_hexp_df$cat_type <- paste0("QUAC_",QUAC_hexp_df[,1],"hexp")
    
    ##separate by loci combination 
    #gSSR genind object
    QUAC_gSSR_genind <- QUAC_garden_wild_genind[loc = QUAC_gSSR_loci]
    
    #create statistical analysis data frame
    QUAC_gSSR_allrich_df <- gather(allelic.richness(QUAC_gSSR_genind)$Ar)
    #calculate hexp and create data frame 
    QUAC_gSSR_hexp <- as.data.frame(cbind(summary(seppop(QUAC_gSSR_genind)[[1]])$Hexp,  
                                          summary(seppop(QUAC_gSSR_genind)[[2]])$Hexp))
    colnames(QUAC_gSSR_hexp) <- c("Garden", "Wild")
    QUAC_gSSR_hexp_df <- gather(QUAC_gSSR_hexp)
    
    #create rownames with sections of names 
    QUAC_gSSR_allrich_df$cat_type <- paste0("QUAC_",QUAC_gSSR_allrich_df[,1], "gSSRallrich")
    QUAC_gSSR_hexp_df$cat_type <- paste0("QUAC_",QUAC_gSSR_hexp_df[,1], "gSSRhexp")
    
    ###calculate diversity stats for all scenarios 
    ##separate by loci combination 
    #EST SSR genind object 
    QUAC_EST_genind <- QUAC_garden_wild_genind[loc = QUAC_EST_loci]
    #calculate allelic richness and create data frame 
    QUAC_EST_allrich_df <- gather(allelic.richness(QUAC_EST_genind)$Ar)
    #calculate hexp and create data frame 
    QUAC_EST_hexp <- as.data.frame(cbind(summary(seppop(QUAC_EST_genind)[[1]])$Hexp,
                                        summary(seppop(QUAC_EST_genind)[[2]])$Hexp))
    colnames(QUAC_EST_hexp) <- c("Garden", "Wild")
    QUAC_EST_hexp_df <- gather(QUAC_EST_hexp)
    
    #create rownames with sections of names 
    QUAC_EST_allrich_df$cat_type <- paste0("QUAC_",QUAC_EST_allrich_df[,1], "ESTallrich")
    QUAC_EST_hexp_df$cat_type <- paste0("QUAC_",QUAC_EST_hexp_df[,1], "ESThexp")
    
    ##combine all categories for statistical tests 
    #all rich
    QUAC_allrich_allLOCI_df <- rbind(QUAC_allrich_df, QUAC_gSSR_allrich_df, QUAC_EST_allrich_df)
    #hexp 
    QUAC_hexp_allLOCI_df <- rbind(QUAC_hexp_df, QUAC_gSSR_hexp_df, QUAC_EST_hexp_df)
    
  }
  
  #next, using ZAIN rebinned, without small pops, for diversity stats 
  if(sp == 4){
  
  #first create a garden/wild genind object 
  ZAIN_garden_genind <- repool(seppop(sp_genind_temp)[1:10])
  levels(ZAIN_garden_genind@pop) <- rep("Garden", 10)
  
  #create wild genind object 
  ZAIN_wild_genind <- repool(seppop(sp_genind_temp)[c(11:19, 23:26, 28:32, 34:35)])
  levels(ZAIN_wild_genind@pop) <- rep("Wild", 20)
  
  #repool to create a genind object 
  ZAIN_garden_wild_genind <- repool(ZAIN_garden_genind, ZAIN_wild_genind)
  
  #calculate allelic richness and create data frame 
  ZAIN_allrich <- gather(allelic.richness(ZAIN_garden_wild_genind)$Ar)
  #calculate hexp and create data frame 
  ZAIN_hexp <- as.data.frame(cbind(summary(seppop(ZAIN_garden_wild_genind)[[1]])$Hexp,
                                     summary(seppop(ZAIN_garden_wild_genind)[[2]])$Hexp))
  colnames(ZAIN_hexp) <- c("Garden", "Wild")
  ZAIN_hexp <- gather(ZAIN_hexp)
  
  #create rownames with sections of names 
  ZAIN_allrich$cat_type <- paste0("ZAIN_",ZAIN_allrich[,1],"allrich")
  ZAIN_hexp$cat_type <- paste0("ZAIN_",ZAIN_hexp[,1], "hexp")
  
  }
  
}

#run kruskal-wallis 
QUAC_aov_test <- aov(QUAC_allrich_allLOCI_df[,2]~QUAC_allrich_allLOCI_df[,3])

#create data frame of allrich for ZAIN 
ZAIN_kw_test <- kruskal.test(ZAIN_allrich[,2]~ZAIN_allrich[,1])

##create a boxplot of just rebinned and woK QUAC
#combined allelic richness data frame 
sp_allrich_df <- rbind(QUAC_allrich_allLOCI_df, ZAIN_allrich)

allrich_final_df <- sp_allrich_df %>%
                      mutate(Species = gsub("_.*","",sp_allrich_df$cat_type))

colnames(allrich_final_df) <- c("Pop_Type", "Allelic_Richness", "Cat", "Species")

allrich_final_df$Cat2 <- factor(allrich_final_df$Cat,     # Reorder factor levels
                                c("QUAC_Gardenallrich", "QUAC_Wildallrich", 
                                  "QUAC_GardengSSRallrich", "QUAC_WildgSSRallrich",
                                  "QUAC_GardenESTallrich", "QUAC_WildESTallrich",
                                  "ZAIN_Gardenallrich", "ZAIN_Wildallrich"))

#write out the plot 
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_allrich.pdf", 
    width = 10, height = 6)
ggplot(allrich_final_df, aes(x=Cat2, y=Allelic_Richness, fill=Pop_Type)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() + xlab("Population Type") + ylab("Allelic Richness") + ylim(0,20) + 
  theme_bw() + facet_wrap(~Species, scale="free") + 
  scale_x_discrete(labels=c("Garden", "Wild", 
                            "Garden_gSSR", "Garden_EST",
                           "Wild_gSSR ", "Wild_EST", 
                           "Garden",
                           "Wild")) + 

  scale_fill_manual(values = c("darkseagreen1", "darkseagreen4"))
dev.off()

#hexp data frame 
sp_hexp_df <- rbind(QUAC_hexp_allLOCI_df, ZAIN_hexp)

#create the output data frame for the results boxplot 
hexp_final_df <- sp_hexp_df %>%
                    mutate(Species = gsub("_.*","", sp_hexp_df$cat_type))

#name the columns for the data frame 
colnames(hexp_final_df) <- c("Pop_Type", "Hexp", "Cat", "Species")

#reorder data frame for the output 
hexp_final_df$Cat2 <- factor(hexp_final_df$Cat,     # Reorder factor levels
                              c("QUAC_Gardenhexp", "QUAC_Wildhexp",
                                "QUAC_GardengSSRhexp", "QUAC_WildgSSRhexp",
                                "QUAC_GardenESThexp", "QUAC_WildESThexp",
                                "ZAIN_Gardenhexp", "ZAIN_Wildhexp"))

#write out the plot 
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_hexp.pdf", 
    width = 10, height = 6)
ggplot(hexp_final_df, aes(x=Cat2, y=Hexp, fill=Pop_Type)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() + xlab("Population Type") + ylab("Expected Heterozygosity") + ylim(0,1) + 
  theme_bw() + facet_wrap(~Species, scale="free") + 
  scale_x_discrete(labels=c("Garden", "Wild", 
                            "Garden_gSSR", "Garden_EST",
                            "Wild_gSSR ", "Wild_EST", 
                            "Garden",
                            "Wild")) +
 scale_fill_manual(values = c("darkseagreen1", "darkseagreen4"))
dev.off()

#write out df 
write.csv(allrich_hexp_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_ZAIN_garden_wild_df.csv")

###############################
#     Clustering Analyses     #
###############################
####Final PCA code for QUAC 
##load in genind objects 
QUAC_genind <- read.genepop("Adegenet_Files/QUAC_woK_allpop_clean.gen",
                            ncode = 3)

#load in data frame woK
QUAC_df <- read.csv("CSV_Files/QUAC_woK_allpop_clean_df.csv")

#reorg the datafile 
QUAC_garden_ind <- sum(table(QUAC_genind@pop)[1:17])
QUAC_garden_genind <- QUAC_genind[1:QUAC_garden_ind,]
levels(QUAC_garden_genind@pop) <- rep("Garden", 17)

#combine with the QUAC wild genind 
QUAC_wild_ind <- sum(table(QUAC_genind@pop)[18:21])
QUAC_wild_genind <- QUAC_genind[(QUAC_garden_ind+1):(QUAC_garden_ind+QUAC_wild_ind),]
levels(QUAC_wild_genind@pop) <- c("Porter","Magazine", "Pryor", "Sugar_Loaf")

#recombine data files 
QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

##PCA code 
QUAC_tab <- tab(QUAC_garden_wild_genind, freq=TRUE, NA.method="mean")

#run PCA
QUAC_PCA <- dudi.pca(QUAC_tab, scale = FALSE, nf = 3, scannf = FALSE)

QUAC_PCA_df <- as.data.frame(cbind(as.numeric(QUAC_PCA$li$Axis1),
                                   as.numeric(QUAC_PCA$li$Axis2)))

QUAC_PCA_df$Pop <- c(QUAC_df$Garden_Wild[1:277],
                     QUAC_df$Pop[278:441])
colnames(QUAC_PCA_df) <- c("Axis1", "Axis2", "Pop", "Pop_Type")

#create pop type column 
QUAC_PCA_df$pop_type <- QUAC_df$Garden_Wild

QUAC_PCA_df$Axis3 <- QUAC_PCA$li$Axis3

#calculate % variation explained by axis 
QUAC_pc1 <- signif(((QUAC_PCA$eig[1])/sum(QUAC_PCA$eig))*100, 3)
QUAC_pc2 <- signif(((QUAC_PCA$eig[2])/sum(QUAC_PCA$eig))*100, 3)
QUAC_pc3 <- signif(((QUAC_PCA$eig[3])/sum(QUAC_PCA$eig))*100, 3)

##ZAIN
pdf("../Analyses/Results/Clustering/QUAC_woK_PCA2_3.pdf", width = 10, height = 8)

ggplot(QUAC_PCA_df, aes(as.numeric(Axis2), as.numeric(Axis3), col = Pop_Type, 
                        shape = Pop)) + 
  geom_point(size = 4) +
  xlab(paste0("PC2 (", QUAC_pc2, "%)")) +
  ylab(paste0("PC3 (", QUAC_pc3, "%)")) + 
  theme_bw() +  
  scale_color_manual(values = c("mediumseagreen", "black")) +
  scale_shape_manual(values = c(16,17,18,3,4))
dev.off()


###Final PCA code for ZAIN 
#load genind object
ZAIN_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen",
                            ncode = 3)



#load data frame 
ZAIN_df <- read.csv("CSV_Files/ZAIN_rebinned_allpop_clean_df.csv")

rownames(ZAIN_genind@tab) <- ZAIN_df$Sample.Name

#load pop df 
ZAIN_pop_df <- read.csv("CSV_Files/ZAIN_pop_df.csv")

##reorg the datafile 
#garden
ZAIN_garden_ind <- sum(table(ZAIN_genind@pop)[1:10])
ZAIN_garden_genind <- ZAIN_genind[1:ZAIN_garden_ind,]
levels(ZAIN_garden_genind@pop) <- rep("Garden", 10)

#wild
ZAIN_wild_ind <- sum(table(ZAIN_genind@pop)[c(11:35)])
ZAIN_wild_genind <- ZAIN_genind[(ZAIN_garden_ind+1):(ZAIN_garden_ind+ZAIN_wild_ind),]

#combine datafile 
ZAIN_garden_wild_genind <- repool(ZAIN_garden_genind,
                                  ZAIN_wild_genind)

#create tab object for genind 
ZAIN_tab <- tab(ZAIN_garden_wild_genind, freq=TRUE, NA.method="mean")

#run PCA
ZAIN_PCA <- dudi.pca(ZAIN_tab, scale = FALSE, nf = 3, scannf = FALSE)

#create PCA data frame 
ZAIN_PCA_df <- as.data.frame(cbind(as.numeric(ZAIN_PCA$li$Axis1), 
                                   as.numeric(ZAIN_PCA$li$Axis2),
                                   as.numeric(ZAIN_PCA$li$Axis3)))
#organize df 
rownames(ZAIN_PCA_df) <- rownames(tab(ZAIN_garden_wild_genind))

#limit df by reduced pops 
ZAIN_red_pop_df <- ZAIN_pop_df[ZAIN_pop_df$Sample.Name %in% rownames(ZAIN_PCA_df),]


#specify wild vs. garden individual
ZAIN_PCA_df$Pop_Type <- ZAIN_red_pop_df$Pop_Type
  
#add variety 
ZAIN_PCA_df$Variety <- ZAIN_red_pop_df$Variety

#name columns of the  PCA data frame 
colnames(ZAIN_PCA_df) <- c("Axis1","Axis2","Axis3","Pop_Type", "Variety")

#calculate % variation explained by axis 
ZAIN_pc1 <- signif(((ZAIN_PCA$eig[1])/sum(ZAIN_PCA$eig))*100, 3)
ZAIN_pc2 <- signif(((ZAIN_PCA$eig[2])/sum(ZAIN_PCA$eig))*100, 3)
ZAIN_pc3 <- signif(((ZAIN_PCA$eig[3])/sum(ZAIN_PCA$eig))*100, 3)

##ZAIN
pdf("../Analyses/Results/Clustering/PCA/ZAIN_PCA2_3.pdf", width = 10, height = 8)

ggplot(ZAIN_PCA_df, aes(Axis1, Axis3, col = Pop_Type, 
                        shape = Variety)) + 
  geom_point(size = 4) +
  xlab(paste0("PC2 (", ZAIN_pc1, "%)")) +
  ylab(paste0("PC3 (", ZAIN_pc3, "%)")) + 
  theme_bw() +  
  scale_color_manual(values = c("mediumseagreen", "black")) +
  scale_shape_manual(values = c(16,18,3))

dev.off()

###############################################################
#     Plotting Barplots of Genetic Diversity Differences      #
###############################################################

#calculate standard errors
allrich_garden_se <- sd(allrich_hexp_df[[1]][,1])/sqrt(length(allrich_hexp_df[[1]][,1]))
allrich_wild_se <- sd(allrich_hexp_df[[2]][,1])/sqrt(length(allrich_hexp_df[[2]][,1]))

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

#########################################
#     Allelic Representation Figure     #
#########################################

#load in data file of the allelic reps 
allsp_wildcap <- read.csv("Garden_Wild_Comparison/QUAC_ZAIN_allelic_rep_transposed.csv")

#convert to a dataframe for ggplot
allsp_wildcap_df <- as.data.frame(allsp_wildcap)

#reorder data frame for graph 
allsp_wildcap_df$allele_type <- factor(allsp_wildcap_df$allele_type,     # Reorder factor levels
                                       c("global", "glob_com",
                                         "glob_lowfr", "glob_rare"))

allsp_wildcap_df$copy_num <- factor(allsp_wildcap_df$copy_num,
                                    c("one_or_more",
                                      "two_or_more",
                                      "five_or_more",
                                      "ten_or_more"))
#mutiply 
allsp_wildcap_df$percent <- (allsp_wildcap_df$per_wild_rep)*100

pdf("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Analyses/Results/Garden_Wild_Comparison/allsp_allelic_rep.pdf", width = 8, height = 6)
ggplot(allsp_wildcap_df, aes(fill = allele_type,
                             x = copy_num,
                             y = percent)) +
  geom_bar(position="dodge",stat="identity") +
  facet_wrap(~species, scale="free") + 
  scale_fill_manual(values = c("red", "darkorange3","coral","deeppink")) + 
  xlab("Allele Copies") +
  ylab("Percent Allelic Representation") +
  theme_bw()
dev.off()

