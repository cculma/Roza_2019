# ld Roza2019

library(ldsep)
library(updog)

setwd("~/Documents/git/big_files/")
load("mout_Roza2019.RData")

msub <- filter_snp(mout, prop_mis < 0.5 & od < 0.010726 & bias < 1.11057) 

ploidy <- 4
gp <- format_multidog(x = msub, varname = paste0("Pr_", 0:ploidy))
class(gp)
dim(gp) # 39611   499     5

ldout <- ldfast(gp = gp, type = "r2")

save(ldout, file = "ldout_Roza2019.RData")
# load("ldout_Roza2019.RData")