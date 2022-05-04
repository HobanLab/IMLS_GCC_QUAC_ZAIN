#####################
#     Libraries     #
#####################

#######################
#     Load files      # 
#######################
setwd("../../Data_Files")

#create a data frame list 
QUAC_wk_list <- list.files(path = "Structure_Files/QUAC_wild_wK", pattern = ".csv")
QUAC_pop_list <- list()

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


