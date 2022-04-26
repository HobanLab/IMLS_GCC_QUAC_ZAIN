#####################
#     Libraries     #
#####################

###################################
#     Plot structure diagrams     #
###################################
# Reading in K7_subset file, an updated version of the CLUMPAK K7.output file. 
# Comes from JuneSubset ipyrad/STRUCTURE run
# This contains samples in their presented order, and has capillaris and one Troy Peak sample removed from file
k7 <- read.table("K7_plot.csv", header=T)

# Make a vector of names
names <- c("cusickiana_SRP","cusickiana_Jarbidge","cusickiana_Owyhee","maguirei","domensis","nevadensis_GRBA","nevadensis_Troy", "parryi")
# Make a vector for name positions on graph
labelPos <- c(16, 32, 37, 50, 71, 82, 88, 94)

# Plotting function for K of any value
plot_k <- function(klist,labelPositions,...){
  # List of colors, which are combinations of RGB components, in hexadecimal  
  colors <- c("#2171B5","#D95F02","#7570B3","#E7298A","#66A61E","#8C510A","#666666")
  # Graphing parameters
  par(mar=c(7,2,10,1)+0.1, mgp = c(3,1,1))
  for(i in 1:length(klist)){
    if(i==1){
      # Initial barplot
      barplot(klist[,i], xlim=c(0,100), horiz=F, beside=F, col=colors[i], axisnames=T, space=0.2, yaxt= "n", main="K=7")
      off.value <- klist[,i]
    }else{
      # Subsequent barplots with offset
      barplot(klist[,i], offset=off.value, add=T, beside=F, xlim=c(0,100), horiz=F, col=colors[i], yaxt= "n")
      off.value <- off.value + klist[,i]
    }
  }
  # y axis
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.50", "0.75", "1.00"), cex.axis = 1, las = 2, pos = -0.2, xpd=T)
  # Add group labels
  text(x=labelPositions, y=-0.031, srt=35, adj=1, xpd=TRUE, labels=names, cex=1.2)
}
dev.off()
plot_k(k7[,4:10], labelPos)
# Plotting lines beneath groups
lineWidth <- 1.9; lineHeight <- rep(-0.015,2)

lines(x = c(0.3,28.3), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(29.4,34.4), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(35.1,39.1), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(40.0,60.7), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(61.5,78.8), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(79.5,84.7), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(85.5,89.5), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)
lines(x = c(90.5,98.0), y = lineHeight, lwd = lineWidth, col = "black", xpd = NA)

# %%% PLOTTING DIFFERENT CLUSTERS FOR K=6 %%%----
# Code for plotting function
plot_k6 <- function(klist,main,...){
  cols <- c('#A8FFFD','#B862D3', '#A39D9D','#FFFF00', '#69C261', '#FF59AC', '#26CDCD',  '#C1C6FF') # Combination of RGB components, in hexadecimal  
  #par(mfrow = c(length(klist),1), mar=c(0.1,0.1,0.1,0.1)+0.1, oma = c(15,0,0,0), mgp = c(1,1,0))
  for(i in 1:length(klist)){
    if(i==1){
      barplot(klist[,i], xlim=c(0,100), horiz=F, beside=F, col=cols[i], axisnames=T, space=0.2, yaxt= "n", main=main)
      off.value <- klist[,i]
    }else{
      barplot(klist[,i], offset=off.value, add=T, beside=F, xlim=c(0,100), horiz=F, col=cols[i], yaxt= "n")
      off.value <- off.value + klist[,i]
    }
  }
  # y axis
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), cex.axis = 1, las = 2, pos = -0.2)
}

# Major cluster
k6 <- read.table("ClumppIndFile.output", header=F)
# Remove useless columns
k6 <-k6[,-(1:5)]
plot_k6(k6, main="K=6 Major Cluster")

# Minor cluster 1
k6 <- read.table("ClumppIndFile.MinorCluster1", header=F)
# Remove useless columns
k6 <-k6[,-(1:5)]
plot_k6(k6, main="K=6 Minor Cluster 1")

# Minor cluster 2
k6 <- read.table("ClumppIndFile.MinorCluster2", header=F)
# Remove useless columns
k6 <-k6[,-(1:5)]
plot_k6(k6, main="K=6 Minor Cluster 2")
