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
##set working directory to load in data files
#setwd("../../Data_Files")
setwd("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Data_Files")

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

#list out the scenarios 
scenarios_list <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned", 
                    "ZAIN_rebinned_sample")

#QUAC loci lists  
QUAC_EST_loci <- c("FIR031", "GOT009", "POR016", "FIR013", "FIR043", "GOTO40", 
                   "PIE039", "FIR53", "FIR048", "PIE125")
QUAC_gSSR_loci <- c("0C11", "1G13", "G07", "1F02","QpZAG9")

QUAC_genind <- read.genepop("Adegenet_Files/QUAC_woK_allpop_clean.gen", ncode = 3)

ZAIN_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen", ncode = 3)

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
for(sp in 1:length(scenarios_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[sp]]), ncode = 3)
  
  #load data frame 
  sp_df_temp <- read.csv(paste0("Data_Frames/", sp_df_list[[sp]]))
  
  #reorganize genind object for QUAC without Kessler for figure 
  if(sp == 2){
    
    ##reorganize individuals into garden and wild pools, but also by loci type 
    #create QUAC garden/wild genind object 
    #first, separate garden genind object 
    QUAC_garden_genind <- repool(seppop(QUAC_genind)[1:17])
    levels(QUAC_garden_genind@pop) <- rep("Garden", 17)
    
    #create wild pop genind 
    QUAC_wild_genind <- repool(seppop(QUAC_genind)[18:21])
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
  ZAIN_garden_genind <- repool(seppop(ZAIN_genind)[1:10])
  levels(ZAIN_garden_genind@pop) <- rep("Garden", 10)
  
  #create wild genind object 
  ZAIN_wild_genind <- repool(seppop(ZAIN_genind)[c(11:19, 23:26, 28:32, 34:35)])
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
QUAC_resampling_df <- read.csv("../Analyses/Results/Garden_Wild_Comparison/QUAC_woK_resampling_df0.csv")
ZAIN_resampling_df <- read.csv("../Analyses/Results/Garden_Wild_Comparison/ZAIN_rebinned_resampling_df0.csv")    

#write PDF with name
pdf("../Analyses/Results/Garden_Wild_Comparison/QUAC_resample_plot_woK_ndrop0.pdf", width = 10, height = 8)
    
    #add points
    plot(QUAC_resampling_df[,2], col = "red", pch = 20, xlab = "Number of Individuals", 
         ylab = "Percent Diversity Capture", xlim = c(0,length(rownames(QUAC_resampling_df))), ylim = c(0,100), cex = 1.2,
         main = "Percent Diversity Capture (All Alleles Included)")
    points(QUAC_resampling_df[,4], col = "darkorange3", pch = 20, cex = 1.2)
    points(QUAC_resampling_df[,5], col = "coral", pch = 20, cex = 1.2)
    points(QUAC_resampling_df[,6], col = "deeppink4", pch = 20, cex = 1.2)
    
    #add line for 95% capture
    abline(h = 95, col = "darkslategray4", lty=2, lwd = 3)
    
    #add legend 
    legend('bottomright', legend = c("Global", "Common", "Low Frequency","Rare", "95% diversity capture"),
           lwd = 2, cex = 1.2, 
           col = c("red", "darkorange3", "coral", "deeppink4", "darkslategray4"), 
           lty = 1)

    
dev.off()
    
pdf("../Analyses/Results/Garden_Wild_Comparison/ZAIN_resample_plot_rebinned_ndrop0.pdf", width = 10, height = 8)
#add points
plot(ZAIN_resampling_df[,2], col = "red", pch = 20, xlab = "Number of Individuals", 
     ylab = "Percent Diversity Capture", xlim = c(0,length(rownames(ZAIN_resampling_df))), ylim = c(0,100), cex = 1.2,
     main = "Percent Diversity Capture (All Alleles Included)")
points(ZAIN_resampling_df[,4], col = "darkorange3", pch = 20, cex = 1.2)
points(ZAIN_resampling_df[,5], col = "coral", pch = 20, cex = 1.2)
points(ZAIN_resampling_df[,6], col = "deeppink4", pch = 20, cex = 1.2)

#add line for 95% capture
abline(h = 95, col = "darkslategray4", lty=2, lwd = 3)

#legend 
legend('bottomright', legend = c("Global", "Common", "Low Frequency","Rare", "95% diversity capture"),
       lwd = 2, cex = 1.2, 
       col = c("red", "darkorange3", "coral", "deeppink4", "darkslategray4"), 
       lty = 1)

dev.off()

#setwd
setwd("../../Data_Files")

##load in all the genind objects
sp_genind_list <- list.files(path = "Adegenet_Files", pattern = "_clean.gen")

#pop list 
pop_list <- list(c(1:17), c(1:17), c(1:10), c(1:10),
                 c(18:22), c(18:21), c(11:35), c(11:35))

#initial lists 
QUAC_EST_loci <- c("FIR031", "GOT009", "POR016", "FIR013", "FIR043", "GOTO40", 
                   "PIE039", "FIR53", "FIR048", "PIE125")
QUAC_gSSR_loci <- c("0C11", "1G13", "G07", "1F02","QpZAG9")

#pvalue data frame 
sp_allrich_hexp_pvalue <- matrix(nrow = length(sp_genind_list), ncol = 2)
colnames(sp_allrich_hexp_pvalue) <- c("All_Rich", "Hexp")
rownames(sp_allrich_hexp_pvalue) <- c("QUAC_wK", "QUAC_woK", "ZAIN_og", "ZAIN_rebinned")

#####################Plot QUAC second try 

#loop to compare diversity capture in wild and botanic garden populations
for(sp in 1:length(sp_genind_list)){
  
  #load genepop files as genind objects 
  sp_genind_temp <- read.genepop(paste0("Adegenet_Files/",sp_genind_list[[1]]), ncode = 3)
  
  ##organize into pops - garden
  #separate into garden genind object 
  sp_garden_genind <- repool(seppop(sp_genind_temp)[pop_list[[1]]])
  #rename pops to be garden only 
  levels(sp_garden_genind@pop) <- rep("Garden", length(pop_list[[1]]))
  
  ##organize into pop types 
  #separate into wild genind object 
  sp_wild_genind <- repool(seppop(sp_genind_temp)[pop_list[[1+4]]])
  #rename to wild only 
  levels(sp_wild_genind@pop) <- rep("Wild", length(pop_list[[1+4]]))
  
  #repool to calculate diversity stats 
  sp_garden_wild_genind <- repool(sp_garden_genind, sp_wild_genind)
  
  #if statement for EST and gSSR comparison - QUAC only, sp == 1, sp == 2
  if(sp == 1|sp == 2){
    
    ###calculate diversity stats for all scenarios 
    ##just determine wild and garden diversity levels 
    #calculate allelic richness and create data frame 
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
    #calculate allelic richness and create data frame 
    sp_gSSR_allrich <- gather(allelic.richness(sp_gSSR_genind)$Ar)
    #calculate hexp and create data frame 
    sp_gSSR_hexp <- as.data.frame(cbind( summary(seppop(sp_gSSR_genind)[[1]])$Hexp,  summary(seppop(sp_gSSR_genind)[[2]])$Hexp))
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
    sp_EST_hexp <- as.data.frame(cbind(summary(seppop(sp_EST_genind)[[1]])$Hexp,  summary(seppop(sp_EST_genind)[[2]])$Hexp))
    colnames(sp_EST_hexp) <- c("Garden", "Wild")
    sp_EST_hexp <- gather(sp_EST_hexp)
    #create rownames with sections of names 
    sp_EST_allrich$cat_type <- paste0(sp_EST_allrich[,1], "_EST_allrich")
    sp_EST_hexp$cat_type <- paste0(sp_EST_hexp[,1], "_EST_hexp")
    
    ##combine all categories for statistical tests 
    #all rich
    sp_allrich_allcat_df <- rbind(sp_allrich_df, sp_gSSR_allrich, sp_EST_allrich)
    #hexp 
    sp_hexp_allcat_df <- rbind(sp_hexp_df, sp_gSSR_hexp, sp_EST_hexp)
    
    #summary data frames 
    sp_allrich_hexp_pvalue[sp,1] <- kruskal.test(sp_allrich_allcat_df[,2]~sp_allrich_allcat_df[,3])[3]$p.value
    sp_allrich_hexp_pvalue[sp,2] <- kruskal.test(sp_hexp_allcat_df[,2]~sp_hexp_allcat_df[,3])[3]$p.value
    
  }else{
    
    #combining into a df 
    allrich_df <- gather(allelic.richness(sp_garden_wild_genind)$Ar)
    
    #run t-test 
    sp_allrich_hexp_pvalue[sp,1] <- kruskal.test(allrich_df[,2]~allrich_df[,1])[3]$p.value
    
    #run hexp code
    sp_hexp <- as.data.frame(cbind(summary(seppop(sp_garden_wild_genind)[[1]])$Hexp,  summary(seppop(sp_garden_wild_genind)[[2]])$Hexp))
    colnames(sp_hexp) <- c("Garden", "Wild")
    #create data frame for hexp
    sp_hexp_df <- gather(sp_hexp)
    #save p-value for hexp 
    sp_allrich_hexp_pvalue[sp,2] <- kruskal.test(sp_hexp_df[,2]~sp_hexp_df[,1])[3]$p.value
    
  }
}

##Plot QUAC allelic richness results 
#set factor levels to compare between  
sp_allrich_allcat_df$cat_type <- factor(sp_allrich_allcat_df$cat_type, levels=c("Garden_allrich", "Wild_allrich", 
                                          "Garden_gSSR_allrich", "Wild_gSSR_allrich",
                                          "Garden_EST_allrich", "Wild_EST_allrich"))

colnames(sp_allrich_allcat_df) <- c("Pop_Type", "All_Rich", "Category")


ggplot(sp_allrich_allcat_df, aes(x=Category, y=All_Rich, fill = Pop_Type))  + geom_boxplot() +
  xlab("Population Type") + ylab("Allelic Richness") + ylim(0,20)  +
  theme_bw() + 
  scale_x_discrete(labels=c("Garden", "Wild", "Garden gSSR", "Wild gSSR", "Garden EST", "Wild EST")) + 
  scale_fill_manual(values = c("darkseagreen1", "darkseagreen4"))


#####################
#     Plotting      #
#####################
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

