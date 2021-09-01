## Script to randomly assign sexes to unsexed individuals for the allele frequency variance models
## Rose Driscoll
## 20201015


# get existing sex data from pedigree
ped<-read.table('FSJpedgeno_Zsexlinked.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
sex_data <- ped[,c(2,5)]
colnames(sex_data) <- c("USFWS", "Sex")

# check for unsexed indivs & assign them a sex
unsexed_indivs <- sex_data[sex_data$Sex==0,]
simulated_sexes <- sample(x = c(1,2), size = nrow(unsexed_indivs), prob = c(0.5,0.5), replace = TRUE)
sex_data[sex_data$Sex==0,"Sex"] <- simulated_sexes

# write to file for later use
write.table(sex_data, "FSJ_sex_data_real_and_simulated_20201015.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
