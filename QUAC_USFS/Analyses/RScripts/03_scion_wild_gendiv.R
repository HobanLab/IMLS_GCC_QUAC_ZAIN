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

#load genind object 
QUAC_gen <- read.genepop("Genepop_Files/USFS_QUAC_gen_clean.gen", ncode = 2)

#load csv file as a data frame 
QUAC_df <- read.csv("CSV_Files/USFS_QUAC_rebinned_df.csv")

#remove individual with too much missing data
QUAC_df <- QUAC_df[QUAC_df$ID %in% rownames(QUAC_gen@tab),]

#create pop list name
pop_list <- unique(QUAC_df$Pop)

#name levels in genind object 
levels(QUAC_gen@pop) <- pop_list

###############################
#     Run gendiv analyses     #
###############################

##reorganize data files
#scions
QUAC_scion_gen <- repool(seppop(QUAC_gen)[1:4])
#levels(QUAC_scion_gen@pop) <- rep("Scions", 4)

#wild 
QUAC_wild_gen <- repool(seppop(QUAC_gen)[5:8])
#levels(QUAC_wild_gen@pop) <- rep("Wild", 4)

#combine
QUAC_scion_wild_gen <- repool(QUAC_scion_gen, QUAC_wild_gen)

##create df for allrich
QUAC_allrich_gather <- rbind(gather(allelic.richness(QUAC_scion_gen)$Ar, na.rm = TRUE),
                            gather(allelic.richness(QUAC_wild_gen)$Ar, na.rm = TRUE))

#allelic richness data frame
QUAC_allrich_gather$Pop_Type <- NA
QUAC_allrich_gather$Pop_Type[1:57] <- "Scion"
QUAC_allrich_gather$Pop_Type[58:117] <- "Wild"
QUAC_allrich_gather$key <- gsub("_.*", "",QUAC_allrich_gather$key)
QUAC_allrich_gather$Analysis <- NA
QUAC_allrich_gather$Analysis <- "All_Rich"
colnames(QUAC_allrich_gather) <- c("Pop", "Stat", "Pop_Type", "Analysis")

col_order <- c("Stat", "Pop", "Pop_Type", "Analysis")

QUAC_allrich_gather <- QUAC_allrich_gather[,col_order]

#create boxplot - write out 
png("../Analyses/Results/QUAC_allrich.png", 
    width     = 5,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4)
ggplot(QUAC_allrich_gather, aes(x = Pop, y = Stat, fill = Pop_Type)) + 
  geom_boxplot() + xlab("Population") + ylab("Allelic Richness") + 
  scale_fill_manual(values = c("darkseagreen2", "darkgreen")) + 
  theme_bw()
dev.off()

#hexp list
hexp_wild_list <- list()
hexp_scion_list <- list()

#create heterozygosity barplot
for(f in 1:length(seppop(QUAC_wild_gen))){
  
  hexp_scion_list[[f]] <- as.data.frame(cbind(summary(seppop(QUAC_scion_gen)[[f]])$Hexp, pop_list[[f+4]],
                                "Scion", "Hexp"))
  #name columns 
  colnames(hexp_scion_list[[f]]) <- c("Stat", "Pop", "Pop_Type", "Analysis")
  #name columns
  hexp_wild_list[[f]] <- as.data.frame(cbind(summary(seppop(QUAC_wild_gen)[[f]])$Hexp, 
                                             pop_list[[f+4]], "Wild", "Hexp"))
  #name columns
  colnames(hexp_wild_list[[f]]) <- c("Stat", "Pop", "Pop_Type", "Analysis")
  
}

#combine all lists
QUAC_hexp_df <- rbind(hexp_scion_list[[1]], hexp_scion_list[[2]],
                      hexp_scion_list[[3]], hexp_scion_list[[4]],
                      hexp_wild_list[[1]], hexp_wild_list[[2]],
                      hexp_wild_list[[3]], hexp_wild_list[[4]])

#convert hexp measurements to numeric
QUAC_hexp_df$Stat <- as.numeric(QUAC_hexp_df$Stat)

#write out barplot
png("../Analyses/Results/QUAC_hexp.png", 
    width     = 5,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4)
ggplot(QUAC_hexp_df, aes(x = Pop, y = Stat, fill = Pop_Type)) + 
  geom_boxplot() + xlab("Population") + ylab("Expected Heterozygosity") + 
  scale_fill_manual(values = c("darkseagreen2", "darkgreen")) + 
  theme_bw()
dev.off()
