rm(list = ls())
setwd("/scratch/user/cesar.medinaculma/20211108_131650/nph/")
library(parallel)
library(doParallel)
library(iterators)
library(foreach)
library(sommer)
library(AGHmatrix)
library(GROAN)

load("WGBLUP_nph.RData")
list2 <- c("nph_38849", "temp_38849", "nph_80177", "temp_80177", "nph_140646", "temp_140646")
# L2 data_8
# L3 Geno

#############

L6_nph <- list()
for (i in 1:length(L3)) {
  for (j in 1:ncol(L2[[i]])) {   
    m <- as.numeric(t(L2[[i]][,j]))
    l <- Gmatrix(L3[[i]], method="VanRaden", ploidy=4, ploidy.correction = T, weights = m)
    L6_nph[[length(L6_nph)+1]] <- l
  }
}

names(L6_nph) <- list2
save(L6_nph, file = "L6_nph.RData")

##############
D1[1:5,1:5]
