###This script visualizes the differences in microsatellite scoring for ZAIN 
##when scoring was performed by different labs. Wild ZAIN microsatellite scoring
#was performed by the XXXXX lab and at a different time than the scoring for
#garden individuals, which were scored by the XXXXX in 2021. Therefore,
#we attempted to determine if there were consistent scoring differences that 
#were due to different times and scorers. We then rebinned scores to create 
#consistency between the data files and created barplots to determine if 
#rebinning analysis lead to greater consistency among individuals. 

#####################
#     Libraries     #
#####################

library(adegenet)
library(diveRsity)

#######################
#     Load Files      #
#######################
setwd("../../Data_Files")

#convert arp2gen - comment in if needed
#arp2gen("ZAIN_adegenet_files/ZAIN_garden_wild_rebinned.arp")

#load in genepop file as genind object
ZAIN_og_genind <- read.genepop("Adegenet_Files/ZAIN_og_allpop.gen", ncode = 3)

#load in file with 2 pops - scored in 2021 (garden) and scored in 2011 (wild)
ZAIN_og_df <- read.csv("Data_Frames/ZAIN_og_allpop_df.csv")

#convert to genpop
ZAIN_og_genpop <- genind2genpop(ZAIN_og_genind)

#load rebinned genepop file as a genind object
ZAIN_rebinned_genind <- read.genepop("Adegenet_Files/Garden_Wild/ZAIN_rebinned_garden_wild.gen", ncode = 3)

#convert to genepop
ZAIN_rebinned_genpop <- genind2genpop(ZAIN_rebinned_genind)

#create a list of loci 
loci <- colnames(ZAIN_og_df)

#clean list to just locus name
loci <- unique(gsub("\\..*","",loci)[4:25])

#######################################
#     Initial Scoring Comparison      #
#######################################
##Code to compare microsatellite scores performed by different labs at different
#time points. 
#Garden individuals were scored in 2021
#Wild individuals were scored in 2010 - 2018

setwd("../Analyses/Results/Scoring_Comparison")
pdf("ZAIN_scoring_comparison_barplots.pdf",width=20,height=9)

#loop to compare scoring between 2021 Scoring and 2011 Scoring
for(a in loci){
  
  #create data frame of each locus
  ZAIN_scoring <- ZAIN_og_genpop[,which(grepl(a,colnames(ZAIN_og_genpop@tab)))]@tab
  
  for(p in 1:2) ZAIN_scoring[p,] <- ZAIN_scoring[p,]/sum(ZAIN_scoring[p,])
  #reorder data frame
  ZAIN_scoring <- ZAIN_scoring[,sort(colnames(ZAIN_scoring))]
  
  #now plot 
  ZAIN_barplot <- barplot(ZAIN_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                           legend.text =  c("2021","2011"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(ZAIN_scoring)))
  
}

dev.off()


############################################################
#     Scoring Comparison Following Rebinning Analysis      #
############################################################
##Code to plot ZAIN genind objects following rebinning analysis was performed
#create pdf 
pdf("ZAIN_scoring_comparison_barplots_post_rebinning.pdf",width=20,height=9)

#loop to compare scoring between 2021 Scoring and 2011 Scoring
for(a in loci){
  
  #create data frame of each locus 
  ZAIN_scoring <- ZAIN_rebinned_genpop[,which(grepl(a,colnames(ZAIN_rebinned_genpop@tab)))]@tab
  
  #loop to standardize score count by percent for each group - wild and garden individuals
  for(p in 1:2) ZAIN_scoring[p,] <- ZAIN_scoring[p,]/sum(ZAIN_scoring[p,])
  
  #reorder data frame for alleles to be plotted in numerical order
  ZAIN_scoring <- ZAIN_scoring[,sort(colnames(ZAIN_scoring))]
  
  #plot barplots by locus
  ZAIN_barplot <- barplot(ZAIN_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                          legend.text =  c("2021","2011"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(ZAIN_scoring)))
  
}

dev.off()

sessionInfo()
