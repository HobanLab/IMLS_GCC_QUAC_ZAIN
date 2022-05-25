##This script is used to determine if there are differences between scoring based
#on who scored them in Zamia integrifolia 

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

#load in data frame with 2 pops - scored by us vs. scored by the Griffith lab
ZAIN_garden_wild_gen <- read.genepop("ZAIN_adegenet_files/ZAIN_garden_wild.gen", ncode = 3)

#convert to genpop
ZAIN_garden_wild_pop <- genind2genpop(ZAIN_garden_wild_gen)

#load in garden/wild data frame 
ZAIN_garden_wild_df <- read.csv("ZAIN_data_frames/ZAIN_allpop_df.csv")

#load rebinned genepop file as a genind object
ZAIN_garden_wild_rebinned_gen <- read.genepop("ZAIN_adegenet_files/ZAIN_garden_wild_rebinned.gen", ncode = 3)

#convert to genepop
ZAIN_garden_wild_rebinned_pop <- genind2genpop(ZAIN_garden_wild_rebinned_gen)

#create a list of loci 
loci <- colnames(ZAIN_garden_wild_df)

#clean list 
loci <- unique(gsub("\\..*","",loci)[4:25])

###############################
#     Scoring Comparison      #
###############################

##loop to compare scoring between the Hoban Lab Scoring and the Griffith Lab
setwd("../Analyses/Results/Scoring_Comparison")
pdf("ZAIN_scoring_comparison_barplots.pdf",width=20,height=9)

for(a in loci){
  
  ZAIN_scoring <- ZAIN_garden_wild_pop[,which(grepl(a,colnames(ZAIN_garden_wild_pop@tab)))]@tab
  
  for(p in 1:2) ZAIN_scoring[p,] <- ZAIN_scoring[p,]/sum(ZAIN_scoring[p,])
  #reorder data frame
  ZAIN_scoring <- ZAIN_scoring[,sort(colnames(ZAIN_scoring))]
  
  #now plot 
  ZAIN_barplot <- barplot(ZAIN_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                           legend.text =  c("Hoban_Lab","Griffith_Lab"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(ZAIN_scoring)))
  
}

dev.off()


############################################################
#     Scoring Comparison Following Rebinning Analysis      #
############################################################
##loop to compare scoring between the Hoban Lab Scoring and the Griffith Lab

pdf("ZAIN_scoring_comparison_barplots_post_rebinning.pdf",width=20,height=9)
for(a in loci){
  
  ZAIN_scoring <- ZAIN_garden_wild_rebinned_pop[,which(grepl(a,colnames(ZAIN_garden_wild_rebinned_pop@tab)))]@tab
  
  for(p in 1:2) ZAIN_scoring[p,] <- ZAIN_scoring[p,]/sum(ZAIN_scoring[p,])
  #reorder data frame
  ZAIN_scoring <- ZAIN_scoring[,sort(colnames(ZAIN_scoring))]
  
  #now plot 
  ZAIN_barplot <- barplot(ZAIN_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                          legend.text =  c("Hoban_Lab","Griffith_Lab"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(ZAIN_scoring)))
  
}

dev.off()


sessionInfo()
