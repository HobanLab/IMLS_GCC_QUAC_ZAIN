#####################
#     Libraries     #
#####################

library(adegenet)

###########################
#     Load Data Files     #
###########################
 
#load in fa sample functions
source("../Analyses/Functions/Fa_sample_funcs.R")

#These numbers designate the populations for wild and garden, for Quercus and Zamia
#Used in code below to subset the genpop objects by wild and garden
#The ZAIN has some wild populations excluded- NOTE come back to this to discuss
garden_pop_numbers<-list(1:17,1:10)
wild_pop_numbers<-list(18:21,c(11:19, 23:26, 28:32, 34:35))
 

 #We  have a loop around all species, NOTE this is so code does not diverge when doing each species separate
 
 for (sp in 1:2){
 	#Right now just doing the first two files
	 gen_inp_filenames<-c("Adegenet_Files/QUAC_woK_allpop_clean.gen", "Adegenet_Files/ZAIN_rebinned_allpop_clean.gen", "Adegenet_Files/ZAIN_rebinned_sample_clean.gen")
	 outp_filenames<-c("QUAC_woK_indiv_rep_percents", "ZAIN_indiv_rep_percents", "ZAIN_red_indiv_rep_percents")
	
	sp_genind <- read.genepop(gen_inp_filenames[sp], ncode = 3)
		 
	#allele categories list
	all_cat_list <-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")
	 
	#vector of duplicate reps (for 1 to 10 individuals)
	 dup_reps <- c(0:9)
	
	#create garden genind. 
	#Note this cannot be done with seppop, because seppop will count only alleles in your subset, thus possibly "losing" wild alleles
	num_garden_ind<-sum(table(sp_genind@pop)[garden_pop_numbers[[sp]]])
	sp_garden_genind <- sp_genind[1:num_garden_ind,]
	#rename pops
	levels(sp_garden_genind@pop) <- rep("Garden",length(garden_pop_numbers[[sp]]))
	
	#create wild genind object 
	#NOTE- This is actually quite tricky because unlike garden populations, wild populations kept are scattered through the dataset for ZAIN
	#So we have to create two vectors- a vector of "starting individuals" and "ending individuals"
	#and use those vectors to populate a list which will go from starting to ending individual of each population, all glued together
	#then you have to unlist this list to make a vector
	wild_ind_list<-list()
	for (i in wild_pop_numbers[[sp]]){
		wild_ind_list[[i]]<-(cumsum(table(sp_genind@pop))-table(sp_genind@pop)+1)[i]:cumsum(table(sp_genind@pop))[i]
	}
	sp_wild_genind <- sp_genind[unlist(wild_ind_list),]
	#rename pops 
	levels(sp_wild_genind@pop) <- rep("Wild",length(wild_pop_numbers[[sp]]))
	
	num_wild_ind <- as.numeric(table(sp_wild_genind@pop))
	
	#convert to the wild genpop object
	sp_wild_genpop <- genind2genpop(sp_wild_genind)
	
	#calculate the allele categories in the wild populations
	sp_all_cat <- get.allele.cat(sp_wild_genpop, 1, 1, num_wild_ind, n_drop = 0, glob_only = TRUE)	
	
	#subset to allele of interest e.g. exlcuding regional alleles 
	sp_all_cat <- sp_all_cat[1:5]
	
	#################################
	#	How many individuals have each allele
	#################################
	
	#create a list to store the number of individuals representing each allele 
	#This list is length of 5, the five allele categories we are concerned with
	#the elements of the list are vectors... the vector is length equal to the number of alleles in each category
	# each element of the vector will be the number of individuals having that allele
	#for example, [[1]][1:3] might be 5,1, 10 which means five individuals, 0 individuals, and 10 individuals have those first three alleles, resepectively
	num_indiv_rep_list <- list(vector(), vector(), vector(), vector(), vector())
	num_indiv_rep_list_he <-  list(vector(), vector(), vector(), vector(), vector())
	num_indiv_rep_list_ho <-  list(vector(), vector(), vector(), vector(), vector())
	
	#This for loop goes through the 5 allele categories of interest
	#It then goes through all the alleles in a category
	#It then counts the number of individuals, the number of homozygotes, and the number of heterozygotes for that allele
	#recall that sp_garden_genind@tab is a matrix of nrows= number of individuals and ncols = number of alleles
	#so every cell of the matrix is an individual-allele combination, and the data in the cell are the number of copies of that allele in that individual
	#An individual can have the allele in 2 copies (homozygote), 1 copy (heterozygote), or not have the allele (0)
	for(cat in 1:length(sp_all_cat)){
	  
	  num_alleles_in_cat <- length(sp_all_cat[[cat]])
	  
	  for (a in 1:num_alleles_in_cat){
	    
	    num_indiv_rep_list[[cat]][a] <- sum(sp_garden_genind@tab[,sp_all_cat[[cat]]][,a] > 0, na.rm=T)	#either he or ho
	    num_indiv_rep_list_he[[cat]][a] <- sum(sp_garden_genind@tab[,sp_all_cat[[cat]]][,a] == 1, na.rm=T) 	#he
	    num_indiv_rep_list_ho[[cat]][a] <- sum(sp_garden_genind@tab[,sp_all_cat[[cat]]][,a] == 2, na.rm=T)	#ho
	  }
	}
	
	############################################################
	# Percent of alleles represented in greater than "dup" number of individuals
	############################################################
	
	#create data frame to save results  
	#In this case the results are the percent of alleles present in greater than "dup" number of individuals
	#The he and ho stand for individuals in the heterozygous and homozygous states
	percent_indiv_results <- matrix(nrow = length(dup_reps),
	                      ncol = length(sp_all_cat))
	
	percent_indiv_results_he <- matrix(nrow = length(dup_reps),
	                         ncol = length(sp_all_cat))
	percent_indiv_results_ho <- matrix(nrow = length(dup_reps),
	                         ncol = length(sp_all_cat))
	
	#This loop goes through the number of "dups" from 1 to 10 (with 1 meaning no "backup")
	#Then through the loop of allele categories (5)
	#Within that loop it determines if the number of individuals with that allele is greater than "dup", 
	#then divides the number of alleles meeting that criteria by the total number of alleles 
	#thus returning the proportion of alleles contained in more than "dup" individuals
	#The he and ho stand for individuals in the heterozygous and homozygous states
	for(dup in dup_reps){
	  for(cat in 1:length(sp_all_cat)){
	    
	    #create data frame to store results 
	    percent_indiv_results[dup+1,cat] <- sum(num_indiv_rep_list[[cat]]>dup)/length(sp_all_cat[[cat]])
	    percent_indiv_results_he[dup+1,cat] <- sum(num_indiv_rep_list_he[[cat]]>dup)/length(sp_all_cat[[cat]])
	    percent_indiv_results_ho[dup+1,cat] <- sum(num_indiv_rep_list_ho[[cat]]>dup)/length(sp_all_cat[[cat]])
	    
	  }
	}
	
	
	#round off
	percent_indiv_results <- signif(percent_indiv_results*100,3)
	percent_indiv_results_he <- signif(percent_indiv_results_he*100,3)
	percent_indiv_results_ho <- signif(percent_indiv_results_ho*100,3)
	
	#label columns and rows
	colnames(percent_indiv_results) <- all_cat_list
	colnames(percent_indiv_results_he) <- all_cat_list
	colnames(percent_indiv_results_ho) <- all_cat_list
	rownames(percent_indiv_results) <- paste0(c(1:10), " or more copies")
	rownames(percent_indiv_results_he) <- paste0(c(1:10), " or more copies")
	rownames(percent_indiv_results_ho) <- paste0(c(1:10), " or more copies")
	
	#save percents as output files
	write.csv(percent_indiv_results, paste("../Analyses/Results/Garden_Wild_Comparison/percent_indiv_results",outp_filenames[sp],".csv"))
	write.csv(percent_indiv_results_he, paste("../Analyses/Results/Garden_Wild_Comparison/percent_indiv_results",outp_filenames[sp],"_he",".csv"))
	write.csv(percent_indiv_results_ho, paste("../Analyses/Results/Garden_Wild_Comparison/percent_indiv_results",outp_filenames[sp],"_ho",".csv"))
}
