#####################
#     Libraries     #
#####################

library(diveRsity)
library(adegenet)
library(poppr)

#########################
#   Load Data Files     #
##########################
#set working directory
setwd("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/QUAC_USFS/Data_Files")

#read in genepop file as a genind object 
USFS_QUAC_gen <- read.genepop("Genepop_Files/USFS_QUAC_genepop.gen", ncode = 2)

#read in USFS QUAC data file 
USFS_QUAC_df <- read.csv("CSV_Files/USFS_QUAC_Scores_df.csv")
USFS_QUAC_df <- USFS_QUAC_df[1:206,]
#pop name list
pop_list <- unique(USFS_QUAC_df$Pop)

#name populations 
levels(USFS_QUAC_gen@pop) <- pop_list 

#create loci list 
loci <- colnames(USFS_QUAC_df)[3:length(colnames(USFS_QUAC_df))]

#clean list to just locus name
clean_loci <- unique(gsub("\\..*","",loci))
clean_loci[4] <- "1F02"

###############################################
#     Binning Reorganizing Genind Objects     #
###############################################
##separate by clones and wild individuals
#Scion individuals genind object
QUAC_scions_gen <- repool(seppop(USFS_QUAC_gen)[1:4])
#rename levels
levels(QUAC_scions_gen@pop) <- rep("Scions", rep(length(levels(QUAC_scions_gen@pop))))

#wild individuals wild genind object
QUAC_wild_gen <- repool(seppop(USFS_QUAC_gen)[5:8])
#rename levels 
levels(QUAC_wild_gen@pop) <- rep("Wild", rep(length(levels(QUAC_wild_gen@pop))))

#repool into one genind object
QUAC_scion_wild_gen <- repool(QUAC_scions_gen, QUAC_wild_gen)

#create genepop object 
QUAC_scion_wild_genpop <- genind2genpop(QUAC_scion_wild_gen)

#######################################
#     Initial Scoring Comparison      #
#######################################
###Scoring barplot to look to see if there are scoring differences 
##Wild individuals were were scored by Emily Schumacher in 2020
#Scions were scored by Emily Schumacher in 2023

#set working directory
setwd("../Analyses/Results/Scoring_Comparison")
pdf("QUAC_scoring_comparison_barplots.pdf",width=20,height=9)

#2023 scoring barplots 
for(a in clean_loci){
  
  #create data frame of each locus
  QUAC_scoring <- QUAC_scion_wild_genpop[,which(grepl(a,colnames(QUAC_scion_wild_genpop@tab)))]@tab
  
  #turn into percentages
  for(p in 1:2) QUAC_scoring[p,] <- QUAC_scoring[p,]/sum(QUAC_scoring[p,])
  
  #reorder data frame
  QUAC_scoring <- QUAC_scoring[,sort(colnames(QUAC_scoring))]
  
  #now plot 
  QUAC_barplot <- barplot(QUAC_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                          legend.text =  c("2023","2020"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(QUAC_scoring)))
  
}

dev.off()

##############3
#     Rebinned
##############
setwd("../../../Data_Files/Genepop_Files")
##load in the rebinned genind object 
QUAC_rebinned_genind <- read.genepop("USFS_QUAC_rebinned_genepop.gen", ncode = 2)

#repool object into wild and garden
#scion
QUAC_rb_scion <- repool(seppop(QUAC_rebinned_genind)[1:4])
#rename levels
levels(QUAC_rb_scion@pop) <- rep("Scion", length(levels(QUAC_rb_scion@pop)))

#wild
QUAC_rb_wild <- repool(seppop(QUAC_rebinned_genind)[5:8])
#rename level
levels(QUAC_rb_wild@pop) <- rep("Wild", length(levels(QUAC_rb_wild@pop)))

#repool scion and garden 
QUAC_rb_scion_wild_genind <- repool(QUAC_rb_scion, QUAC_rb_wild)

#convert to genind object
QUAC_rebinned_genepop <- genind2genpop(QUAC_rb_scion_wild_genind)

##Code to plot ZAIN genind objects following rebinning analysis was performed
#create pdf 
setwd("../../Analyses/Results/Scoring_Comparison")
pdf("QUAC_scoring_comparison_barplots_post_rebinning.pdf",width=20,height=9)

#loop to compare scoring between 2021 Scoring and 2011 Scoring
for(a in clean_loci){
  
  #create data frame of each locus 
  QUAC_scoring <- QUAC_rebinned_genepop[,which(grepl(a,colnames(QUAC_rebinned_genepop@tab)))]@tab
  
  #loop to standardize score count by percent for each group - wild and garden individuals
  for(p in 1:2) QUAC_scoring[p,] <- QUAC_scoring[p,]/sum(QUAC_scoring[p,])
  
  #reorder data frame for alleles to be plotted in numerical order
  QUAC_scoring <- QUAC_scoring[,sort(colnames(QUAC_scoring))]
  
  #plot barplots by locus
  QUAC_rebinned_barplot <- barplot(QUAC_scoring, las = 2, beside = TRUE, col = c("darkgreen", "darkseagreen1"),
                          legend.text =  c("2023","2020"), ylim = c(0,1), main = paste0(a),
                          names = gsub("^.*\\.","",colnames(QUAC_scoring)))
  
}
dev.off()
