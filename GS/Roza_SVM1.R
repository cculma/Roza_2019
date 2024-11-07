
rm(list = ls())
setwd("/home/samac/medin297/msi/Roza2019/WGBLUP")
# setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/GP1/GROAN/")
library(GROAN)
library(AGHmatrix)
library(data.table)
library(tidyverse)

load("MT_GBLUP.RData")

y <- Y2 %>% dplyr::select(may_20, INDIV)

#--------- TUNE ---------
#tuning the SVR on a grid of hyperparameters
results.tune = phenoRegressor.SVR(
  phenotypes = y[,1],
  genotypes = G2,
  covariances = NULL,
  extraCovariates = NULL,
  mode = 'tune',
  kernel = 'linear', cost = 10^(-3:+3) #SVR-specific parameters
)

save(results.tune, file = "results.tune.RData")
