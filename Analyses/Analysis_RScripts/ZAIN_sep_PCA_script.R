
#
setwd("C:/Users/eschumacher/Documents/GitHub/GCC_QUAC_ZAIN/Data_Files")

#load genind object
ZAIN_rebinned_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen", ncode = 3)

#load data frame 
ZAIN_rebinned_df <- read.csv("Data_Frames/ZAIN_rebinned_allpop_clean_df.csv")

#reorg 
rownames(ZAIN_rebinned_genind@tab) <- ZAIN_rebinned_df$Sample.Name
levels(ZAIN_rebinned_genind@pop) <- unique(ZAIN_rebinned_df$Pop)
wild_pops <- ZAIN_rebinned_df[ZAIN_rebinned_df$Garden_Wild == "Wild",]$Pop

#save just wild pops 
ZAIN_wild_rebinned_genind <- repool(seppop(ZAIN_rebinned_genind)[11:35])

#create tab object for genind 
ZAIN_wild_rebinned_tab <- tab(ZAIN_wild_rebinned_genind, freq=TRUE, NA.method="mean")

#run PCA
ZAIN_PCA <- dudi.pca(ZAIN_wild_rebinned_tab, scale = FALSE, nf = 2, scannf = FALSE)

#create PCA data frame 
ZAIN_PCA_df <- as.data.frame(cbind(as.numeric(ZAIN_PCA$li$Axis1), 
                                      as.numeric(ZAIN_PCA$li$Axis2), wild_pops))

colnames(ZAIN_PCA_df) <- c("Axis1","Axis2","Pop")

#calculate % variation explained by axis 
ZAIN_pc1 <- signif(((ZAIN_PCA$eig[1])/sum(ZAIN_PCA$eig))*100, 3)
ZAIN_pc2 <- signif(((ZAIN_PCA$eig[2])/sum(ZAIN_PCA$eig))*100, 3)

##replace colors in the PCA document 
for(1:nrow(ZAIN_PCA_df)){
  
  #repalce names 
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFB",][,3] <- "Bronson"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFBP1",][,3] <- "Bulow1"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFBP2",][,3] <- "Bulow2"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFCC",][,3] <- "Cayo Costa"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFCNS",][,3] <- "Canaveral"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZCF",][,3] <- "Chapman Field"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFCR",][,3] <- "Crystal River"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFENP",][,3] <- "Everglades"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFFDSP",][,3] <- "Faver-Dykes"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFI",][,3] <- "Ichetucknee"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFK",][,3] <- "Koreshan"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFLB",][,3] <- "Broward"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFLTI",][,3] <- "Little Talbot Island"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFONF1",][,3] <- "Ocala1"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFONF2",][,3] <- "Ocala2"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFO",][,3] <- "Oscar Scherer"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFBP",][,3] <- "Palm Beach"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFSW1",][,3] <- "Suwannee1"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFSW2",][,3] <- "Suwannee2"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFTS",][,3] <- "Tide Swamp"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFTSP1",][,3] <- "Tomoka1"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFTSP2",][,3] <- "Tomoka2"
  ZAIN_PCA_df[ZAIN_PCA_df$Pop == "ZFWB",][,3] <- "Werner-Boyce"
  
}


umbrosa_pops <- c("Bulow1", "Bulow2","Canaveral",
                  "Faver-Dykes", "Little Talbot Island",
                  "Ocala1", "Ocala2", "Tomoka1", "Tomoka2")

for(p in 1:length(umbrosa_pops)){
  
    ZAIN_PCA_df[ZAIN_PCA_df$Pop == umbrosa_pops[[p]],]$variety <- "umbrosa"
  
}

pdf("../Analyses/Results/Clustering/PCA/ZAIN_wildpop_PCA.pdf", width = 8,
    height = 6)
ggplot(ZAIN_PCA_df, aes(as.numeric(Axis1), as.numeric(Axis2), col = variety)) + 
  geom_vline(xintercept=0,lwd=1,colour="black") +
  geom_hline(yintercept=0,lwd=1,colour="black") +
  geom_point(size = 3) + 
  xlab(paste0("PC1 (", ZAIN_pc1, "%)")) +
  ylab(paste0("PC2 (", ZAIN_pc2, "%)")) + 
  theme_bw() +  
  scale_color_manual(values=c("dodgerblue","cadetblue3"))
dev.off()
