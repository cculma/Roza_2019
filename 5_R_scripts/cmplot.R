rm(list = ls())
library(tidyverse)
library(rMVP)
library(CMplot)
library(goeveg)
setwd("~/Documents/git/big_files/")
a1 <- read.csv("Roza2019_06_GWASPoly.txt")

a1[1:5,1:5]
a2 <- a1[,c(1:3)]

a3 <- subset(a2, 
             Chrom == "Chr1" |
               Chrom == "Chr2" |
               Chrom == "Chr3" |
               Chrom == "Chr4" |
               Chrom == "Chr5" |
               Chrom == "Chr6" |
               Chrom == "Chr7" |
               Chrom == "Chr8")

a4 <- subset(a2, 
             Chrom == "Chr6")

CMplot(a3,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       main="illumilla_60K",file.output=TRUE,verbose=TRUE,width=9,height=6)

CMplot(a4,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       main="illumilla_60K",file.output=F,verbose=TRUE,width=9,height=6)

nrow(a3)-nrow(a2)


CMplot(a2,type="p",plot.type="d",bin.size=1e6,
       file.output=F,verbose=T)

head(a2)
cc <- dplyr::count(a2, Chrom)

# difference between markers

a3 <- a2 %>%
  group_by(Chrom) %>%
  mutate(pos1 = Position - lag(Position, default = Position[1]))

a4 <- a3 %>% dplyr::filter(pos1 != 0)

a4 <- a4 %>% group_by(Chrom) %>% summarise_at(vars(pos1), list(min = min, max = max, mean = mean))

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/VCF_Monoploid/")

write.csv(a4, "gap.csv", quote = F, row.names = F)                                              
