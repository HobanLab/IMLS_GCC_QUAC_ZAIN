###################################
#     Trying with Sean's Code     #
###################################
#create garden genind
num_garden_ind<-sum(table(QUAC_woK_genind@pop)[1:17])

QUAC_garden_genind <- QUAC_woK_genind[1:num_garden_ind,]

#rename pops
levels(QUAC_garden_genind@pop) <- rep("Garden",17)

#create wild genind object 
num_wild_ind <- sum(table(QUAC_woK_genind@pop)[18:21])

QUAC_wild_genind <- QUAC_woK_genind[(num_garden_ind+1):(num_garden_ind+num_wild_ind),]

#rename pops 
levels(QUAC_wild_genind@pop) <- rep("Wild",4)

#recombine into garden/wild genind object
#QUAC_garden_wild_genind <- repool(QUAC_garden_genind, QUAC_wild_genind)

#convert to the wild genpop object
QUAC_wild_genpop <- genind2genpop(seppop(QUAC_wild_genind)[[2]])

#calculate how alleles are represented ex situ
QUAC_all_rep <- colSums(QUAC_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
QUAC_all_cat <- get.allele.cat(QUAC_wild_genpop, 1, 1, num_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
QUAC_all_cat <- QUAC_all_cat[1:5]

##create a list to store the individual numbers 
#list 

num_rep_list <- list(list(), list(), list(), list(), list())

for(cat in 1:length(QUAC_all_cat)){
  
  num_alleles_in_cat <- length(QUAC_all_cat[[cat]])
  
  for (a in 1:num_alleles_in_cat){
    
    num_rep_list[[cat]][a] <- sum(QUAC_garden_genind@tab[,QUAC_all_cat[[cat]]][,a] > 0, na.rm=T)
    
  }
}

#create data frame to save results  
QUAC_rep_df <- matrix(nrow = length(dup_reps),
                      ncol = length(QUAC_all_cat))

for(dup in dup_reps){
  for(cat in 1:length(QUAC_all_cat)){
    
    #create data frame to store results 
    QUAC_rep_df[dup+1,cat] <- sum(num_rep_list[[cat]]>dup)/length(QUAC_all_cat[[cat]])
    
    
  }
}

QUAC_rep_df <- signif(QUAC_rep_df*100,3)
colnames(QUAC_rep_df) <- all_cat_list
rownames(QUAC_rep_df) <- paste0(c(1:10), " or more copies")

write.csv(QUAC_rep_df, "../Analyses/Results/Garden_Wild_Comparison/QUAC_rep_df.csv")

####ZAIN 
#load in ZAIN data file
ZAIN_genind <- read.genepop("Adegenet_Files/ZAIN_rebinned_allpop_clean.gen",
                            ncode = 3)

#create garden genind
ZAIN_garden_ind <- sum(table(ZAIN_genind@pop)[1:10])

ZAIN_garden_genind <- ZAIN_genind[1:ZAIN_garden_ind,]

#create wild genind object 
ZAIN_wild_ind <- sum(table(ZAIN_genind@pop)[c(11:19, 23:26, 28:32, 34:35)])

ZAIN_wild_genind <- ZAIN_genind[(ZAIN_garden_ind+1):(ZAIN_garden_ind+ZAIN_wild_ind),]

#convert to the wild genpop object
ZAIN_wild_genpop <- genind2genpop(ZAIN_wild_genind)

#calculate how alleles are represented ex situ
ZAIN_all_rep <- colSums(ZAIN_garden_genind@tab,na.rm=T)

#calculate the allele categories in the wild populations
ZAIN_all_cat <- get.allele.cat(ZAIN_wild_genpop, 1, 1, ZAIN_wild_ind, n_drop = 0, glob_only = TRUE)	

#remove regional alleles 
ZAIN_all_cat <- ZAIN_all_cat[1:5]

##create a list to store the individual numbers 
#list 
ZAIN_num_rep_list <- list(list(), list(), list(), list(), list())

for(cat in 1:length(ZAIN_all_cat)){
  
  ZAIN_num_alleles_in_cat <- length(ZAIN_all_cat[[cat]])
  
  for (a in 1:ZAIN_num_alleles_in_cat){
    
    ZAIN_num_rep_list[[cat]][a] <- sum(ZAIN_garden_genind@tab[,ZAIN_all_cat[[cat]]][,a] > 0, na.rm=T)
    
  }
}

#create data frame to save results  
ZAIN_rep_df <- matrix(nrow = length(dup_reps),
                      ncol = length(ZAIN_all_cat))

for(dup in dup_reps){
  for(cat in 1:length(ZAIN_all_cat)){
    
    #create data frame to store results 
    ZAIN_rep_df[dup+1,cat] <- sum(ZAIN_num_rep_list[[cat]]>dup)/length(ZAIN_all_cat[[cat]])
    
    
  }
}
colnames(ZAIN_rep_df) <- all_cat_list
rownames(ZAIN_rep_df) <- paste0(c(1:10)," or more copies")

ZAIN_rep_df <- signif(ZAIN_rep_df*100,3)

write.csv(ZAIN_rep_df, "../Analyses/Results/Garden_Wild_Comparison/ZAIN_rep_df.csv")


