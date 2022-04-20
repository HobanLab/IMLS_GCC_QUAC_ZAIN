#######################
#---#DETERMINE WHAT ALLELES FALL IN WHAT CATEGORIES---#
#######################

get.allele.cat <- function(UK_genpop, region_makeup, N_REGIONS, n_ind_p_pop, n_drop=2,glob_only=F){

	n_pops<-length(n_ind_p_pop)
	#--Set up categories for GLOBAL ALLELES
	allele_freqs<-colSums(UK_genpop@tab)/(sum(n_ind_p_pop)*2)	
	min_freq<-(n_drop/(2*sum(n_ind_p_pop)))
	glob<-as.vector(which(allele_freqs>min_freq))
	glob_v_com<-as.vector(which(allele_freqs>0.10))
	glob_com<-as.vector(which(allele_freqs>0.05))
	glob_lowfr<-as.vector(which(allele_freqs<0.10&allele_freqs>0.01))
	glob_rare<-as.vector(which((allele_freqs<0.01)&(allele_freqs>min_freq)))

	if (!(glob_only)) {

	#--REGIONAL ALLELES--#	#NOT USED IN WHITE ASH MANUSCRIPT- USED IN OTHER PROJECTS
	reg_com<-vector(length=length(allele_freqs));	reg_rare<-vector(length=length(allele_freqs))
	reg_com_int<-NULL #These are the interesting ones- common in one region, rare in others
	allele_freq_by_reg<-matrix(nrow=N_REGIONS,ncol=length(allele_freqs))

	#this will go region by region and determine allele frequencies in that region 
	#then it determines for which regions the allele frequency is >0.20
	#then it gives each allele a "+1" when it is at >0.20 in a region
	#alleles for which this count==1 at the end, are only present at >0.20 for one region
	for (nr in 1:N_REGIONS){ 
		#get the number of individuals in the region (sum across populations in the region)
		n_ind_this_region<- sum(n_ind_p_pop[region_makeup[nr][[1]]])	#note that the "1" just gets us the list item, its ok
		#divide allele counts by number of individuals (times 2 because diploid)
		allele_freq_by_reg[nr,]<-colSums(UK_genpop[region_makeup[[nr]],]@tab)/(n_ind_this_region*2)
		#record if above threshold for common
		reg_com<- as.numeric(reg_com) + as.numeric(allele_freq_by_reg[nr,]>=0.20)
	#check: print(as.vector(which((allele_freq_by_reg[nr,]>=thresh_freq_H))))
		#record if below threshold for rare
		reg_rare<- as.numeric(reg_rare) + as.numeric(allele_freq_by_reg[nr,]<0.10)
	}  
	#this is the number of alleles that are regionally common
	#check: which(reg_com==1); which(reg_rare==N_REGIONS)
	#check: ((reg_com==1)|(reg_com==2)|(reg_com==3))&((reg_rare==5)|(reg_rare==4))
	reg_com_int<-which((reg_com==1)&(reg_rare>=(N_REGIONS-3)))
	#check: allele_freq_by_reg[,reg_com_int]


	#--LOCAL ALLELES--#
	#get the populations having alleles which occur at frequency greater than a given frequency in focal pop'n
	#and not more than another frequency in more than two other populations
	local_freqs<-UK_genpop@tab/as.vector(n_ind_p_pop*2)
	#identify columns (alleles) for which at least one population exceeds min threshold
	loc_com_d1<-as.vector(which((colMax(as.data.frame((local_freqs)))>=0.20)&(apply(local_freqs,2,sort)[(n_pops*.95),]<0.10),arr.ind=F))
	loc_com_d2<-as.vector(which((colMax(as.data.frame((local_freqs)))>=0.15)&(apply(local_freqs,2,sort)[(n_pops*.90),]<0.05),arr.ind=F))
	loc_rare<-as.vector(which((colMax(as.data.frame((local_freqs)))<0.10)&(colMax(as.data.frame((local_freqs)))>min_freq)))
	} else reg_com_int<-NA; loc_com_d1<-NA; loc_com_d2<-NA; loc_rare<-NA
	list(glob, glob_v_com, glob_com, glob_lowfr, glob_rare, reg_com_int, loc_com_d1, loc_com_d2, loc_rare)
}		


#####################
#--GENERAL SAMPLING FUNCTION--#
#####################

#we pass it a vector of populations and a vector of sample sizes, as well as the genind object
sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
	for (p in 1:length(vect_pop_ID))
	alleles[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),])
	alleles
}
