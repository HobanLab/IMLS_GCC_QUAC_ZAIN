###This code is used to generate structure diagrams for wild and botanic garden
##individuals for both QUAC and ZAIN. 

#####################
#     Libraries     #
#####################

#######################
#     Load files      # 
#######################
setwd("../../Data_Files")

#create a data frame list 
QUAC_wk_list <- list.files(path = "Structure_Files/QUAC/QUAC_wild_wK", pattern = ".csv")
QUAC_woK_list <- list.files(path = "Structure_Files/QUAC/QUAC_wild_woK", pattern = ".csv")
QUAC_garden_wild_wk_list <- list.files(path = "Structure_Files/QUAC/QUAC_allpop_wK", pattern = ".csv")
QUAC_garden_wild_woK_list <- list.files(path = "Structure_Files/QUAC/QUAC_allpop_woK", pattern = ".csv")
ZAIN_garden_wild_og_list <- list.files(path = "Structure_Files/ZAIN/ZAIN_og_allpop_STR", pattern = ".csv")

#####################################
#     Create structure diagrams     #
#####################################
##QUAC wild pops with Kessler  
pop_names <- c("Porter", "Magazine", "Pryor", "Sugar Loaf","Kessler")
label_pos <- c(25,55,92,140,173)

for(k in 1:length(QUAC_wk_list)){
  
  QUAC_k <- read.csv(paste0("Structure_Files/QUAC_wild_wK/",QUAC_wk_list[[k]]))
  
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  colors <- viridis::viridis(k+1)
  
  pdf(paste0("../Analyses/Results/Clustering/QUAC_wild_wK_STR_k",k+1,".pdf"), width = 20, height = 10)
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))

  for(i in 1:length(QUAC_k_ready)){
    
    if(i==1){
    
    #Initial barplot
    barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
            border = NA)
   off.value <- QUAC_k_ready[,i]
  }else{
    # Subsequent barplots with offset
   barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
           border = NA)
    off.value <- off.value + QUAC_k_ready[,i]
    
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
    
    
      #}
    text(x=label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=pop_names, cex=1.2)  
    
    }
  }
  
   dev.off()
}

##QUAC structure diagram without Kessler 
QUAC_pop_names_woK <- c("Porter", "Magazine", "Pryor", "Sugar Loaf")
label_pos <- c(25,55,92,140)

for(k in 1:length(QUAC_garden_wild_wk_list)){
  
  QUAC_k <- read.csv(paste0("Structure_Files/QUAC_wild_woK/",QUAC_woK_list[[k]]))
  
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  colors <- viridis::viridis(k+1)
  
  pdf(paste0("../Analyses/Results/Clustering/QUAC_wild_woK_STR_k",k+1,".pdf"), width = 20, height = 10)
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  for(i in 1:length(QUAC_k_ready)){
    
    if(i==1){
      
      #Initial barplot
      barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
              border = NA)
      off.value <- QUAC_k_ready[,i]
    }else{
      # Subsequent barplots with offset
      barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
              border = NA)
      off.value <- off.value + QUAC_k_ready[,i]
     
      
      #}
      
    }
  }
  
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  text(x=label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=QUAC_pop_names_woK, cex=1.2)  
  dev.off()
}

##QUAC structure diagram with Kessler and garden individuals
#QUAC_garden_wild_wK <- c("AP", "Arnold", "Bartlett",
#                         "CBG", "DBG", "FGT",
 #                        "HBG","Meis","MBG", "Moore",
#                         "MA", "Peckerwood",
#                         "TMA","TTA", "USNA","UWBG",
#                         "ZBGP", "Porter", "Magazine", 
#                         "Pryor", "Sugar Loaf","Kessler")
#label_pos <- c(2, 15, 29, )

##QUAC structure diagram with Kessler population individuals 
#loop to write out all k values 
for(k in 1:length(QUAC_garden_wild_wk_list)){
  
  QUAC_k <- read.csv(paste0("Structure_Files/QUAC_allpop_wK/",QUAC_garden_wild_wk_list[[k]]))
  
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  colors <- viridis::viridis(k+1)
  
  pdf(paste0("../Analyses/Results/Clustering/QUAC_garden_wild_wK_STR_k",k+1,".pdf"), width = 20, height = 10)
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  for(i in 1:length(QUAC_k_ready)){
    
    if(i==1){
      
      #Initial barplot
      barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
              border = NA)
      off.value <- QUAC_k_ready[,i]
    }else{
      # Subsequent barplots with offset
      barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
              border = NA)
      off.value <- off.value + QUAC_k_ready[,i]
      
      
      #}
      
    }
  }
  
 # axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  #text(x=label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=QUAC_pop_names_woK, cex=1.2)  
  dev.off()
}

##QUAC structure diagram for all pops (garden and wild) without Kessler pop 
#loop to write out all k values 
for(k in 1:length(QUAC_garden_wild_woK_list)){
  
  QUAC_k <- read.csv(paste0("Structure_Files/QUAC/QUAC_allpop_woK/",QUAC_garden_wild_woK_list[[k]]))
  
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  colors <- viridis::viridis(k+1)
  
  pdf(paste0("../Analyses/Results/Clustering/Structure/QUAC/QUAC_garden_wild_woK_STR_k",k+1,".pdf"), width = 20, height = 10)
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  for(i in 1:length(QUAC_k_ready)){
    
    if(i==1){
      
      #Initial barplot
      barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
              border = NA)
      off.value <- QUAC_k_ready[,i]
    }else{
      # Subsequent barplots with offset
      barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
              border = NA)
      off.value <- off.value + QUAC_k_ready[,i]
      
      
      #}
      
    }
  }
  
  #axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  #text(x=label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=QUAC_pop_names_woK, cex=1.2)  
  dev.off()
}

##ZAIN wild and botanic garden structure diagrams 
for(k in 1:length(ZAIN_garden_wild_og_list)){
  
  ZAIN_k <- read.csv(paste0("Structure_Files/ZAIN/ZAIN_og_allpop_STR/",ZAIN_garden_wild_og_list[[k]]))
  
  ZAIN_k_ready <- ZAIN_k[,-c(1:2)]
  
  colors <- viridis::viridis(k+1)
  
  pdf(paste0("../Analyses/Results/Clustering/Structure/ZAIN/ZAIN_garden_wild_og_k",k+1,".pdf"), width = 20, height = 10)
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  for(i in 1:length(ZAIN_k_ready)){
    
    if(i==1){
      
      #Initial barplot
      barplot(ZAIN_k_ready[,i], xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
              border = NA)
      off.value <- ZAIN_k_ready[,i]
    }else{
      # Subsequent barplots with offset
      barplot(ZAIN_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
              border = NA)
      off.value <- off.value + ZAIN_k_ready[,i]
      
      
      #}
      
    }
  }
  
  #axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  #text(x=label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=ZAIN_pop_names_woK, cex=1.2)  
  dev.off()
}
