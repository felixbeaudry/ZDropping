load(ldhatFile) #percontig_LDhat_rho.Rdata

#subsample
LD.list.ordered_sub <- LD.list.ordered[[1]][seq(1, nrow(LD.list.ordered[[1]]), 100), ]

#convert into df
for(i in c(2:length(LD.list.ordered))){
  LD.list.ordered_sub <- rbind.data.frame(LD.list.ordered_sub,
                                          LD.list.ordered[[i]][seq(1, nrow(LD.list.ordered[[i]]), 100), ]
  )
}

save(LD.list.ordered_sub,file="FSJ2.rho.rdata")
