## geneDrop_runner_RD.py
## python script to trim pedigree, make input files for geneDrop_RD.c, and run geneDrop_RD.c (compiled as a.out)
## run: python geneDrop_runner_final.py [ped file] [cohort file] [prefix for output files] [number of drops] [filter ungenotyped individuals?] [output type] [XY, ZW, or autosome data?] [proportion of unsexed individuals to assign as male]

# import packages
import numpy as np
import pandas as pd
import subprocess
import random

import time
import sys

# make it so there's no warning about chained assignment
pd.set_option('mode.chained_assignment', None)

# Function to read pedigree and cohort files and return the pedigree and cohort data frames
def load_files(ped_filename, cohort_filename):
    # read the pedigree file into a data frame
    pedigree = pd.read_csv(ped_filename, sep=' ', names=None, header=None)
    # read the cohort file into a data frame
    cohort = pd.read_csv(cohort_filename, sep=' ', names=None, header=None)
    # return the two data frames
    return pedigree, cohort

# Function to perform biased coinflips to assign the sex of unsexed individuals
def assign_sex(prop_males):
	# get a random number between 0 and 1
	random_number = random.random()
	# assign male if random number is less than the proportion of males, otherwise assign female
	if random_number < prop_males:
		sex = 1 # 1 = male
	else:
		sex = 2 # 2 = female
	return sex

# Function to make a dictionary of alleles
def check_alleles_v2(p_all, founders, a1, a2):
    # get a numpy ndarray of the unique values of the a1 and a2 alleles 
    # .values tells it to operate on just the values, not axes lables
    # .ravel('K') tells it to read them in the order they're stored in memory
    alleles = pd.unique(p_all[[a1, a2]].values.ravel('K'))
    # convert each item in the array to string
    alleles_str = map(str, alleles)
    # make a dictionary of alleles
    # keys: allele in string format
    # values: allele in string format...
    allele_mapping = dict(zip(alleles_str,alleles))
    return allele_mapping
        

# get command line arguments
ped_fname = sys.argv[1] ## ped file
coh_fname = sys.argv[2] ## cohort file
out_prefix = sys.argv[3] ## prefix of file
nDrops = sys.argv[4] ## number of drops for each SNP
useUngenotypedIndFilter = sys.argv[5] ## whether to count ungenotyped individuals when doing drops (default no)
outputType = sys.argv[6] ## l for long, s for short
chromType = sys.argv[7] ## X for sex chromosomes from an XY system, Z for sex chromosomes from a ZW system, A for autosomes
unsexedPropMale = float(sys.argv[8]) ## the proportion of unsexed individuals that should be assigned as males

# put a dot after the out_prefix so the filenames will work
out_prefix = out_prefix + '.'

# print names of the pedigree and cohort files to the terminal
print ped_fname, coh_fname
# load pedigree and cohort data usign load_files function
p, c = load_files(ped_fname, coh_fname)
# convert the ID numbers of nestlings to string format
c[0] = c[[0]].applymap(str)

# assign unsexed individuals a sex based on user-provided proportion male
for i in range(len(p[4])): # go through all individuals
	if p[4][i] == 0: # if an individual has a sex of 0 (i.e., unsexed)
		p[4][i] = assign_sex(unsexedPropMale) # assign a sex based on a biased coinflip that produces males in the proportion specified by the user

## make file containing allele mapping ##

allele_mapping_str = ''
## columns representing first of two loci in the ped file (columns 6, 8, 10, ..)
loci = list(p)[6::2]

# need to convert from int to str because allele can be either letter or number
p[[1,2,3,4]] = p[[1,2,3,4]].applymap(str)

for locus in loci:
    a1 = locus
    a2 = locus + 1
    # need to convert from int to str because allele can be either letter or number
    p[[a1,a2]] = p[[a1,a2]].applymap(str)
    # use loc indexer to pull out ungenotyped founders
    founders_ng = p.loc[( (p[2] == '0') | (p[3] == '0') ) & (p[a1] == '0') & (p[a2] == '0') ]
    # select just the IDs (column 1) of these non-genotyped founders
    ids = founders_ng[1]
    # also put the IDs in a variable called `all_ng_founders`
    all_ng_founders = ids
    # make a copy of the pedigree (p) so that we can modify the copy without changing the original
    pt = p.copy()    
    # Now trim the pedigree
    # while there are IDs in the `ids` variable (.size gets the total number of items in the array)
    while(ids.size > 0):
    # make all founder_ng kids founders
        # if an individual's mom or dad is in the list of ungenotyped founders, 
        # mark this individual as a founder by changing the IDs of its parents to 0
        pt[2].loc[(pt[2].isin(ids)) | (pt[3].isin(ids))] = '0'
        pt[3].loc[(pt[2].isin(ids)) | (pt[3].isin(ids)) | (pt[2] == '0')] = '0'
        # remove all founder_ng from previous generation
        pt = pt.loc[ ~pt[1].isin(ids)]
        # pull out all *new* non-genotyped founders
        ids = pt[1].loc[( (pt[2] == '0') | (pt[3] == '0') ) & (pt[a1] == '0') & (pt[a2] == '0') ]
        # append this latest group of IDs to the list of all founder_ngs
        all_ng_founders = all_ng_founders.append(ids)    
        # this while loop ends when it can't find any more ungenotyped founders
    # founders of trimmed pedigree (both parents unknown in trimmed pedigree, genotype known)
    founders = pt.loc[( (pt[2] == '0') | (pt[3] == '0') ) ]
    # all (non-founder) individuals in trimmed pedigree (descendents of trimmed pedigree founders)
    in_pedigree = pt.loc[ ~pt[1].isin(founders[1])]
    # combine founders and in_pedigree to get the full list of who is in the pedigree
    p_all = pd.concat([founders, in_pedigree])

## Now we have a trimmed pedigree! 
## Next, make the input files for geneDrop_final_NC.c

    # this sums up the number of individuals that have each allele (0, 1, or 2) but it doesn't save anything so I don't know what it's doing
    p_all[a1].value_counts().add(p_all[a2].value_counts(),fill_value=0)
    # this does the same for just the founders
    founders[a1].value_counts().add(founders[a2].value_counts(),fill_value=0)
    # now run check_alleles_v2 to get a dictionary of the unique alleles
    allele_mapping = check_alleles_v2(p_all, founders, a1, a2)
    # this is a duplicate of line 125 above and I still don't think it's doing anything!
    p_all[a1].value_counts().add(p_all[a2].value_counts(),fill_value=0)
    # for the C output, we want to reorder the columns to have ID, a1, a2, dad, mom, sex
    # and renumber the rows in order
    c_output = p_all[[1,a1,a2,2,3,4]].reset_index(drop=True)
    # if an individual is a founder (parent IDs 0) then change the parent IDs to -1
    c_output.replace({2:'0',3:'0'},-1,inplace=True)
    # add the column name "id" to the first column
    c_output.rename(index=int, columns={1:'id'},inplace=True)
    # in the cohort data, add the column name "id" to the first column
    c.rename(index=int, columns={0:'id'},inplace=True)
    # reshape cohort data to have individuals along the left side and cohort years across the top
    cp = c.pivot(index='id', columns=1, values=1)
    # get rid of the index row that just says "id"
    cp.index.name = None
    # add row numbers to the left of the IDs and rename the ID column "index"
    cp.reset_index(inplace=True)
    # add a 'c' to the start of each cohort year 
    # this way we will never have an ID and a year that match (that could mess up the merge later)
    cp.columns = ['c'+str(i) for i in list(cp)]
    cp.rename(columns={'cindex':'id'},inplace=True)
    # create a column at the far left of the data frame that lists the cohort year(s) for each individual
    # I am not entirely clear on how this works
    cp['cohort'] = cp[cp.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(int).astype(str)),axis=1)
    # also create a column that says how many cohorts each individual is in
    cp['cohortN'] = cp.shape[1] - cp.isnull().sum(axis=1) - 2
    # now add in the allele and parent ID data from c_output, to the left of all of the cohort year columns but to the right of the id column
    result = pd.merge(c_output, cp, how='left', on='id')
    # find the total number of cohorts based on the difference between the most recent year and the earliest year in the data
    nCohorts = int(result[list(result)[6:-2]].max().max() - result[list(result)[6:-2]].min().min() + 1)
    # if an individual doesn't have a cohort year listed (i.e., if cohort is NaN) then make both cohort and cohortN 0 
    result['cohortN'].fillna('0', inplace=True)
    result['cohort'].fillna('0', inplace=True)
    # pull out just ID, allele 1, allele 2, dad ID, mom ID, sex, cohortN, and cohort
    result = result[['id', a1, a2, 2, 3, 4, 'cohortN', 'cohort']]
    # redo indexes so that they go from 0 to length of individuals
    # make a dictionary of every individual's ID and index
    d = dict(zip(result['id'], result.index))
    # replace the dad IDs with indices
    result.replace({2: d}, inplace=True)
    # replace the mom IDs with indices
    result.replace({3: d}, inplace=True)
    # replace the individual IDs with indices
    result['id'] = result.index

## create c file
    # get list of all alleles
    xAlleles = pd.unique(result[[a1, a2]].values.ravel('K'))
    ## special case: if everyone is genotyped, you still need to include the '0' allele
    # check if '0' is already in the list of alleles:
    if '0' in xAlleles:
        # If it is, you have all the alleles
        # calculate how many alleles this is.
        nAlleles = pd.unique(result[[a1, a2]].values.ravel('K')).size
    else:
        # If '0' isn't in your list of alleles, add 1 to the number of alleles you have to get total number of alleles
        nAlleles = pd.unique(result[[a1, a2]].values.ravel('K')).size + 1
    # count how many individuals are in the tree
    nIndividuals = result.shape[0]
    # get a non-redundant list of cohort years
    existingCohorts = result['cohort'].unique()[1:].astype(int)
    # sort the list of years
    existingCohorts.sort()
    # count up how many years/cohorts that is
    existingCohortsNumber = str(existingCohorts.size)
    # make C header with num cohorts in first column, num alleles in second, num individuals in third, and num drops in fourth
    c_header = str(nCohorts) + ',' + str(nAlleles) + ',' + str(nIndividuals) + ',' + str(nDrops) + '\n'
    # second line is the number of cohorts followed by the cohort years
    c_line2 = existingCohortsNumber + ',' + ','.join(existingCohorts.astype(str)) + '\n'
    # make filename: out_prefix followed by allele number plus .drop suffix
    fname = out_prefix + str( (int(a1) - 4)//2 ) +'.drop'
    # write c_header and c_line2 to file, followed by the `result` data frame in csv format
    with open(fname,'w') as f:
        f.write(c_header)
        f.write(c_line2)
        result.to_csv(f,columns=None,index=False,header=None,sep=',')
    f.close()

## Finally, actually call the C script

    # if user has requested that output be in long format:
    if outputType == 'l':
        # call a.out (the compiled C code) on this file with long format option
        # tell the user you're getting the simulation results
        print 'creating .drop.sim.txt file for .pedped file columns ', a1,a2, ' corresponding to locus ', str((int(a1) - 4)//2) 
        # assemble the command to run the gene dropping simulation and calculate allele frequencies
        callstr = './a.out ' + fname + ' ' + useUngenotypedIndFilter + ' R L ' + chromType + ' > ' + fname + '.sim.txt\n'
        # print to stdout
        print(callstr)
        # actually call the command
        subprocess.call(callstr, shell=True)
        # tell the user you're getting the data results
        print 'creating .drop.data.txt file for .pedped file columns ', a1,a2, ' corresponding to locus ', str((int(a1) - 4)//2) 
        # assemble the command to just calculate allele frequencies based on data, not run gene dropping simulation
        callstr = './a.out ' + fname + ' ' + useUngenotypedIndFilter + ' D L ' + chromType + ' > ' + fname + '.data.txt\n'
        # print to stdout
        print(callstr)
        # actually call the command
        subprocess.call(callstr, shell=True)
    # if user has requested output in short format or not requested anything:
    else:
        # call a.out (the compiled C code) on this file without long format option
        # tell the user you're getting the simulation results
        print 'creating .drop.sim.txt file for .pedped file columns ', a1,a2, ' corresponding to locus ', str((int(a1) - 4)//2) 
        # assemble the command to run the gene dropping simulation and calculate allele frequencies
        callstr = './a.out ' + fname + ' ' + useUngenotypedIndFilter + ' R S ' + chromType + ' > ' + fname + '.sim.txt\n'
        # print to stdout
        print callstr
        # actually call the command
        subprocess.call(callstr, shell=True)
        # tell the user you're getting the data results
        print 'creating .drop.data.txt file for .pedped file columns ', a1,a2, ' corresponding to locus ', str((int(a1) - 4)//2) 
        # assemble the command to just calculate allele frequencies based on data, not run gene dropping simulation
        callstr = './a.out ' + fname + ' ' + useUngenotypedIndFilter + ' D S ' + chromType + ' > ' + fname + '.data.txt\n'
        # print to stdout
        print callstr
        # actually call the command
        subprocess.call(callstr, shell=True)

    # add the names of the two files and list of the alleles and the allele mapping to the allele_mapping_str variable
    allele_mapping_str = allele_mapping_str + fname + '.sim.txt ' + fname + '.data.txt ' +\
    (','.join("{},{}".format(k, v) for k, v in allele_mapping.items())) + '\n'

# save everything from allele_mapping_str in a .info.txt file
allele_mapping_fname = out_prefix + 'info.txt'
with open(allele_mapping_fname, 'w') as f:
    f.write(allele_mapping_str)
f.close()
