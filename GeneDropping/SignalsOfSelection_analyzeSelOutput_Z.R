#script to analyze gene dropping results for tests of selection
#Rose Driscoll and Nancy Chen
#Last updated: 19 January 2024


library(dplyr)
library(plyr)

## First, the sex-linked SNPs

snplist<-data.frame(SNP=1:250, Chr=rep("Z", times=249),stringsAsFactors=FALSE)

#net change in allele freq between 1999-2013
pval_1999to2013_sexlinked<-data.frame(snp=snplist$SNP,stringsAsFactors=FALSE)

#look at obs vs sim change in allele freq between adjacent years over time
pval_change_sexlinked<-data.frame(snp=rep(snplist$SNP,each=14),year=rep(c(1999:2012),length(snplist$SNP)),stringsAsFactors=FALSE)

for (i in snplist$SNP)
{
  # read in data
	obs<-read.table(file=paste('working_files/intermediate_files/seltest_Zsexlinked.',i,'.drop.data.txt',sep=''),header=TRUE)
	sim<-read.table(file=paste('working_files/intermediate_files/seltest_Zsexlinked.',i,'.drop.sim.txt',sep=''),header=TRUE)
	# get the simulation results for just 1999-2013 and divide by the total number of alleles each year
	simfreq<-mapply('/',sim[,11:25],obs[obs$allele==2 & obs$cohort_year>1998,'all_alleles_count'])

	# test for net selection between 1999-2013
	# calculate the difference between the observed allele frequency in 2013 and observed allele frequency in 1999
	obsDelta1999<-obs[obs$allele==2 & obs$cohort_year==2013,'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==1999,'frequency_of_allele']
	# calculate the difference between the simulated allele frequency in 2013 and observed allele frequency in 1999 (for each of 1,000,000 simulations)
	simDelta1999<-simfreq[,15]-simfreq[,1]
	# add the observed change and median simulated change to `pval_1999to2013_sexlinked` table
	pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'obsChange']<-obsDelta1999
	pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'simChange']<-median(simDelta1999)
	# calculate the difference between the median simulated change and the simulated change from each simulation, as well as the observed change
	simdiff1999<-simDelta1999 - median(simDelta1999)
	obsdiff1999<-obsDelta1999 - median(simDelta1999)

	# if the observed change was smaller than the median simulated change
	if (obsdiff1999 < 0)
	{
	  # we will say that this is the - direction
		pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'dir']<-'-'
		# calculate the p value: proportion of simulations where the simulated difference was smaller than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
		pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'pval']<-sum(simdiff1999<obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
	} 
	# if the observed change was larger than the median simulated change
	else 
	{
	  # we will say that this is the + direction
		pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'dir']<-'+'
		# calculate the p value: proportion of simulations where the simulated difference was larger than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
		pval_1999to2013_sexlinked[pval_1999to2013_sexlinked$snp==i,'pval']<-sum(simdiff1999>obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
	}

	#test for selection in adjacent years
	# for the 14 pairs of adjacent years we are looking at (1999-2000, 2000-2001, ..., 2012-2013)
	for(x in c(1:14)) 
	{
	  # get the first year in the pair
		yr<-1998 + x
		# calculate the difference between the observed allele frequency in the second year and observed allele frequency in the first year
		obsDelta<-obs[obs$allele==2 & obs$cohort_year==(yr+1),'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==yr,'frequency_of_allele']
		# calculate the difference between the simulated allele frequency in the second year and observed allele frequency in the first year (for each of 1,000,000 simulations)
		simDelta<-simfreq[,x+1]-simfreq[,x]
		# add the observed change and median simulated change to `pval_change_sexlinked` table
		pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'obsChange']<-obsDelta
		pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'simChange']<-median(simDelta)
		# calculate the difference between the median simulated change and the simulated change from each simulation, as well as the observed change
		simdiff<-simDelta - median(simDelta)
		obsdiff<-obsDelta - median(simDelta)

		# if the observed change was smaller than the median simulated change
		if (obsdiff < 0)
		{
		  # we will say that this is the - direction
		  pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'dir']<-'-'
		  # calculate the p value: proportion of simulations where the simulated difference was smaller than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
		  pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'pval']<-sum(simdiff<obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
		} 
		# if the observed change was larger than the median simulated change
		else 
		{
		  # we will say that this is the + direction
		  pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'dir']<-'+'
		  # calculate the p value: proportion of simulations where the simulated difference was larger than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
		  pval_change_sexlinked[pval_change_sexlinked$snp==i & pval_change_sexlinked$year==yr,'pval']<-sum(simdiff>obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
		}
	}

}

## Now the pseudoautosomal SNPs

snplist<-data.frame(SNP=1:19, Chr=rep("Z", times=19),stringsAsFactors=FALSE)

#net change in allele freq between 1999-2013
pval_1999to2013_pseudoautosomal<-data.frame(snp=snplist$SNP,stringsAsFactors=FALSE)

#look at obs vs sim change in allele freq between adjacent years over time
pval_change_pseudoautosomal<-data.frame(snp=rep(snplist$SNP,each=14),year=rep(c(1999:2012),length(snplist$SNP)),stringsAsFactors=FALSE)

for (i in snplist$SNP)
{
  # read in data
  obs<-read.table(file=paste('working_files/intermediate_files/seltest_Zpseudoautosomal.',i,'.drop.data.txt',sep=''),header=TRUE)
  sim<-read.table(file=paste('working_files/intermediate_files/seltest_Zpseudoautosomal.',i,'.drop.sim.txt',sep=''),header=TRUE)
  # get the simulation results for just 1999-2013 and divide by the total number of alleles each year
  simfreq<-mapply('/',sim[,11:25],obs[obs$allele==2 & obs$cohort_year>1998,'all_alleles_count'])
  
  # test for net selection between 1999-2013
  # calculate the difference between the observed allele frequency in 2013 and observed allele frequency in 1999
  obsDelta1999<-obs[obs$allele==2 & obs$cohort_year==2013,'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==1999,'frequency_of_allele']
  # calculate the difference between the simulated allele frequency in 2013 and observed allele frequency in 1999 (for each of 1,000,000 simulations)
  simDelta1999<-simfreq[,15]-simfreq[,1]
  # add the observed change and median simulated change to `pval_1999to2013_pseudoautosomal` table
  pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'obsChange']<-obsDelta1999
  pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'simChange']<-median(simDelta1999)
  # calculate the difference between the median simulated change and the simulated change from each simulation, as well as the observed change
  simdiff1999<-simDelta1999 - median(simDelta1999)
  obsdiff1999<-obsDelta1999 - median(simDelta1999)
  
  # if the observed change was smaller than the median simulated change
  if (obsdiff1999 < 0)
  {
    # we will say that this is the - direction
    pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'dir']<-'-'
    # calculate the p value: proportion of simulations where the simulated difference was smaller than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
    pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'pval']<-sum(simdiff1999<obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
  } 
  # if the observed change was larger than the median simulated change
  else 
  {
    # we will say that this is the + direction
    pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'dir']<-'+'
    # calculate the p value: proportion of simulations where the simulated difference was larger than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
    pval_1999to2013_pseudoautosomal[pval_1999to2013_pseudoautosomal$snp==i,'pval']<-sum(simdiff1999>obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
  }
  
  #test for selection in adjacent years
  # for the 14 pairs of adjacent years we are looking at (1999-2000, 2000-2001, ..., 2012-2013)
  for(x in c(1:14)) 
  {
    # get the first year in the pair
    yr<-1998 + x
    # calculate the difference between the observed allele frequency in the second year and observed allele frequency in the first year
    obsDelta<-obs[obs$allele==2 & obs$cohort_year==(yr+1),'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==yr,'frequency_of_allele']
    # calculate the difference between the simulated allele frequency in the second year and observed allele frequency in the first year (for each of 1,000,000 simulations)
    simDelta<-simfreq[,x+1]-simfreq[,x]
    # add the observed change and median simulated change to `pval_change_pseudoautosomal` table
    pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'obsChange']<-obsDelta
    pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'simChange']<-median(simDelta)
    # calculate the difference between the median simulated change and the simulated change from each simulation, as well as the observed change
    simdiff<-simDelta - median(simDelta)
    obsdiff<-obsDelta - median(simDelta)
    
    # if the observed change was smaller than the median simulated change
    if (obsdiff < 0)
    {
      # we will say that this is the - direction
      pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'dir']<-'-'
      # calculate the p value: proportion of simulations where the simulated difference was smaller than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
      pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'pval']<-sum(simdiff<obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
    } 
    # if the observed change was larger than the median simulated change
    else 
    {
      # we will say that this is the + direction
      pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'dir']<-'+'
      # calculate the p value: proportion of simulations where the simulated difference was larger than the observed difference plus half the proportion of simulations where the simulated and observed differences were the same
      pval_change_pseudoautosomal[pval_change_pseudoautosomal$snp==i & pval_change_pseudoautosomal$year==yr,'pval']<-sum(simdiff>obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
    }
  }
  
}

## Combine and save tables

# re-number pseudoautosomal snps to come after sex-linked snps
pval_1999to2013_pseudoautosomal$snp <- pval_1999to2013_pseudoautosomal$snp + 250
pval_change_pseudoautosomal$snp <- pval_change_pseudoautosomal$snp + 250

# combine tables
pval_change <- bind_rows(pval_change_sexlinked, pval_change_pseudoautosomal)
pval_1999to2013 <- bind_rows(pval_1999to2013_sexlinked, pval_1999to2013_pseudoautosomal)

# save to file
save(pval_change,pval_1999to2013,file="working_files/intermediate_files/SignalsOfSelection_Results.rdata")

