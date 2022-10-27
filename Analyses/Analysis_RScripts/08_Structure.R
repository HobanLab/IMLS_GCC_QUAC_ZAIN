###This code is used to generate structure diagrams for wild and botanic garden
##individuals for both QUAC and ZAIN. 

#####################
#     Libraries     #
#####################

library(randomcoloR)
library(RColorBrewer)

#######################
#     Load files      # 
#######################
setwd("../Analyses/Results/Clustering/Structure")

##create lists of Q values for each k clustering value 
#QUAC k-value lists with Kessler pop
QUAC_wk_garden_wild_klist <- list.files(path = "QUAC/QUAC_allpop_wK", pattern = ".csv")
#QUAC k-value lists without Kessler pop
QUAC_woK_garden_wild_klist <- list.files(path = "QUAC/QUAC_allpop_woK", pattern = ".csv")
#QUAC k-value lists, wild only with Kessler
QUAC_wk_wild_klist <- list.files(path = "QUAC/QUAC_wild_wK", pattern = ".csv")
#QUAC k-value lists, wild only without Kessler
QUAC_woK_wild_klist <- list.files(path = "QUAC/QUAC_wild_woK", pattern = ".csv")

###lists for ZAIN - only best supported k Q values are used for ZAIN 
##due to the large numbe rof populations 
#ZAIN k-value lists, original scores garden and wild 
ZAIN_og_garden_wild_klist <- list.files(path = "ZAIN/ZAIN_og_allpop_STR", pattern = ".csv")
#ZAIN k-value lists, rebinned garden and wild Q values  
ZAIN_rebinned_garden_wild_list <- list.files(path = "ZAIN/ZAIN_rebinned_allpop_str", pattern = ".csv")
#ZAIN k-value lists, wild only 
ZAIN_rebinned_wild_klist <- list.files(path = "ZAIN/ZAIN_rebinned_wild", pattern = ".csv")

#####################################
#     Create structure diagrams     #
#####################################

####QUAC structure, wild populations with garden pops  
###Structure diagram with kessler mountain pop 

#create lists with label names and position
QUAC_wK_pop_names <- c("Garden","Porter", "Magazine", "Pryor", "Sugar Loaf","Kessler")
QUAC_wK_label_pos <- c(150,316, 352, 391, 431, 468)

#loop to create structure diagrams for each k clustering value 
for(k in 1:length(QUAC_wk_garden_wild_klist)){
  
  #load in list with each k value
  QUAC_k <- read.csv(paste0("QUAC/QUAC_allpop_wK/",QUAC_wk_garden_wild_klist[[k]]))
  
  #remove columns with population and individual names  
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  #generate colors for each cluster
  colors <- viridis::viridis(k+1)
  
  ##create PDF for the structure diagram
  pdf(paste0("QUAC/QUAC_garden_wild_wK_STR_k",k+1,".pdf"), width = 20, height = 10)
  #indicate graphing margins 
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))

  #loop for each column of the structure diagram 
  for(i in 1:length(QUAC_k_ready)){
    
    #first column 
    if(i==1){
      
    #create barplot
    barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
            border = NA)
    #offset each column after the first column  
    off.value <- QUAC_k_ready[,i]
    
    }else{
    
    #add all other columns in the barplot
    barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
           border = NA)
    #offset each column after the first column  
    off.value <- off.value + QUAC_k_ready[,i]
    
    #set up axis dimensions
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
    
    #now add population labels to the x-axis
    text(x=QUAC_wK_label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=QUAC_wK_pop_names, cex=1.2)  
    
    }
  }
  
   dev.off()
}

##Structure diagram for QUAC k values without Kessler Mountain pop
#create lists with label names and position
QUAC_woK_pop_names <- c("Garden","Porter", "Magazine", "Pryor", "Sugar Loaf")
QUAC_woK_label_pos <- c(150,316, 352, 391, 431)

#loop to create structure diagrams for each k clustering value 
for(k in 1:length(QUAC_woK_garden_wild_klist)){
  
  #load in list with each k value
  QUAC_k <- read.csv(paste0("QUAC/QUAC_allpop_woK/",QUAC_woK_garden_wild_klist[[k]]))
  
  #remove columns with individual and population names 
  QUAC_k_ready <- QUAC_k[,-c(1:2)]
  
  #generate colors for each cluster
  colors <- viridis::viridis(k+1)
  
  ##create PDF for QUAC woK structure diagrams 
  pdf(paste0("QUAC/QUAC_garden_wild_woK_STR_k",k+1,".pdf"), width = 20, height = 10)
  
  #indicate graphing margins 
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  #loop to generate structure diagram
  for(i in 1:length(QUAC_k_ready)){
    
    #create first column 
    if(i==1){
      
      #Initial barplot
      barplot(QUAC_k_ready[,i], xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, beside=F, 
              col=colors[i], axisnames=T, space=0, yaxt= "n", main=paste0("K = ", k+1),
              border = NA)
      
      #create first column offset 
      off.value <- QUAC_k_ready[,i]
      
    }else{
      
      #create columns other than first column
      barplot(QUAC_k_ready[,i], offset=off.value, add=T, beside=F, xlim=c(0,length(QUAC_k_ready[,1])), horiz=F, col=colors[i], yaxt= "n", space=0,
              border = NA)
      #offset value from the first column 
      off.value <- off.value + QUAC_k_ready[,i]
    }
  }

  #set up axis dimensions
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  
  #label axes
  text(x=QUAC_woK_label_pos, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=QUAC_woK_pop_names, cex=1.2)  
  
  dev.off()
}

#####ZAIN og structure diagram, garden and wild  

#create lists with population labels
ZAIN_og_pop_names <- unique(read.csv(paste0("ZAIN/ZAIN_og_allpop_STR/", ZAIN_og_garden_wild_klist[[1]]))[,2])
#create list for label positions
ZAIN_label_positions <- c(190,313,322,326,331,341,355,367,376,387,406,
                          442,483,512,563,603,640,686,718,737,752,767,
                          788,815,859,884,899,900, 912,946,984,1031,1066,
                          1083,1107,1146)

#loop to create struture diagrams for all delta K clustering values 
for(k in 1:length(ZAIN_og_garden_wild_klist)){
  
  #load in k Q value list
  ZAIN_k <- read.csv(paste0("ZAIN/ZAIN_og_allpop_STR/",ZAIN_og_garden_wild_klist[[k]]))
  
  #remove columns with ind and pop names 
  ZAIN_k_ready <- ZAIN_k[,-c(1:2)]
  
  #generate colors for clusters
  colors <- randomColor(length(ZAIN_k_ready))

  ##create PDF for structure diagrams 
  pdf(paste0("ZAIN/ZAIN_og_garden_wild_STR_k",length(ZAIN_k_ready),".pdf"), width = 20, height = 10)
  
  #indicate graphing margins
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  #loop to create structure diagrams
  for(i in 1:length(ZAIN_k_ready)){
    
    #create first column 
    if(i==1){
      
      #generate barplot with the first column 
      barplot(ZAIN_k_ready[,i], xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, 
              beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", 
              main=paste0("K = ", length(ZAIN_k_ready)),border = NA)
      #create offset for other columns 
      off.value <- ZAIN_k_ready[,i]
    }else{
      #create all other columns 
      barplot(ZAIN_k_ready[,i], offset=off.value, add=T, beside=F, 
              xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, col=colors[i], 
              yaxt= "n", space=0,border = NA)
      off.value <- off.value + ZAIN_k_ready[,i]
      
    }
  }
  
  #label axes with populations
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), 
       cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  text(x=ZAIN_label_positions, y=-0.031, srt=35, adj=1, xpd=TRUE, 
       labels=ZAIN_og_pop_names, cex=1.2)  
  dev.off()
}

####ZAIN structure for rebinned data 
#create lists with population labels
ZAIN_rebinned_pop_names <- unique(read.csv(paste0("ZAIN/ZAIN_rebinned_allpop_str/", ZAIN_rebinned_garden_wild_list[[1]]))[,2])
#remove NA from this list 
ZAIN_rebinned_pop_names <- ZAIN_rebinned_pop_names[c(1,3:36)]
#create list for label positions
ZAIN_rebinned_label_positions <- c(190,313,322,326,331,341,355,367,376,387,406,
                          442,483,522,553,573,603,650,687,721,760,802,843,
                          888,912,937,951,975,1003,1016,1030,1065,1105,1125,1148)

#loop to create struture diagrams for all delta K clustering values 
for(k in 1:length(ZAIN_rebinned_garden_wild_list)){
  
  #load in k Q value list
  ZAIN_k <- read.csv(paste0("ZAIN/ZAIN_rebinned_allpop_str/",ZAIN_rebinned_garden_wild_list[[k]]))
  
  #remove columns with ind and pop names 
  ZAIN_k_ready <- ZAIN_k[,-c(1:2)]
  
  #generate colors for clusters
  colors <- distinctColorPalette(length(ZAIN_k_ready))
  
  ##create PDF for structure diagrams 
  pdf(paste0("ZAIN/ZAIN_rebinned_garden_wild_STR_k",length(ZAIN_k_ready),".pdf"), width = 20, height = 10)
  
  #indicate graphing margins
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  #loop to create structure diagrams
  for(i in 1:length(ZAIN_k_ready)){
    
    #create first column 
    if(i==1){
      
      #generate barplot with the first column 
      barplot(ZAIN_k_ready[,i], xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, 
              beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", 
              main=paste0("K = ", length(ZAIN_k_ready)),border = NA)
      #create offset for other columns 
      off.value <- ZAIN_k_ready[,i]
    }else{
      #create all other columns 
      barplot(ZAIN_k_ready[,i], offset=off.value, add=T, beside=F, 
              xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, col=colors[i], 
              yaxt= "n", space=0,border = NA)
      off.value <- off.value + ZAIN_k_ready[,i]
      
    }
  }
  
  #label axes with populations
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), 
       cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  text(x=ZAIN_rebinned_label_positions, y=-0.031, srt=35, adj=1, xpd=TRUE, 
       labels=ZAIN_rebinned_pop_names, cex=1.2)  
  dev.off()
}

##ZAIN wild pops only structure diagrams 
#create lists with population labels
ZAIN_rebinned_wild_pop_names <- unique(read.csv(paste0("ZAIN/ZAIN_rebinned_wild/", ZAIN_rebinned_wild_klist[[1]]))[,2])
#remove NA from this list 
ZAIN_rebinned_wild_pop_names <- ZAIN_rebinned_wild_pop_names[c(1,3:26)]
#create list for label positions
ZAIN_rebinned_wild_label_positions <- c(19,52,89,130,160,181,216,257,294,334,
                                        366,408,452,490,524,544,560,582,610,624,
                                        638,665,707,733,756)

#loop to create struture diagrams for all delta K clustering values 
for(k in 1:length(ZAIN_rebinned_wild_klist)){
  
  #load in k Q value list
  ZAIN_k <- read.csv(paste0("ZAIN/ZAIN_rebinned_wild/",ZAIN_rebinned_wild_klist[[k]]))
  
  #remove columns with ind and pop names 
  ZAIN_k_ready <- ZAIN_k[,-c(1:2)]
  
  #generate colors for clusters
  colors <- distinctColorPalette(length(ZAIN_k_ready))
  
  ##create PDF for structure diagrams 
  pdf(paste0("ZAIN/ZAIN_rebinned_wild_STR_k",length(ZAIN_k_ready),".pdf"), width = 20, height = 10)
  
  #indicate graphing margins
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  
  #loop to create structure diagrams
  for(i in 1:length(ZAIN_k_ready)){
    
    #create first column 
    if(i==1){
      
      #generate barplot with the first column 
      barplot(ZAIN_k_ready[,i], xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, 
              beside=F, col=colors[i], axisnames=T, space=0, yaxt= "n", 
              main=paste0("K = ", length(ZAIN_k_ready)),border = NA)
      #create offset for other columns 
      off.value <- ZAIN_k_ready[,i]
    }else{
      #create all other columns 
      barplot(ZAIN_k_ready[,i], offset=off.value, add=T, beside=F, 
              xlim=c(0,length(ZAIN_k_ready[,1])), horiz=F, col=colors[i], 
              yaxt= "n", space=0,border = NA)
      off.value <- off.value + ZAIN_k_ready[,i]
      
    }
  }
  
  #label axes with populations
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), 
       cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  text(x=ZAIN_rebinned_wild_label_positions, y=-0.031, srt=35, adj=1, xpd=TRUE, 
       labels=ZAIN_rebinned_wild_pop_names, cex=1.2)  
  dev.off()
}

