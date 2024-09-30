library(AGHmatrix)

# setwd("~/Documents/git/big_files/")
setwd("/home/samac/medin297/msi/Roza2019/WGBLUP")

pheno <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", row.names = 1)
head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1

G1 <- read.table("Roza2019_06_GS.txt", header = TRUE, row.names = 1, check.names = F, sep = "\t")
# G1[1:5,1:5]
dim(G1)
have.both = intersect(rownames(pheno), rownames(G1))
G2 <- G1[have.both,]
# G2[1:5,1:5]
dim(G2)
pheno1 <- pheno[have.both,]
colnames(pheno1)
pheno1 <- pheno1[1:(length(pheno)-1)]

order1 <- match(rownames(pheno), rownames(G2))
G2  <- G2[order1,]
# rownames(G2)[1:10]
# rownames(pheno1)[1:10]

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/GP1/")
Fst <- read.csv("Fst.csv", row.names = 1)
cor(Fst)
order2 <- match(rownames(Fst), colnames(G2))
Fst  <- Fst[order2,]
# rownames(Fst)[1:10]
# colnames(G2)[1:10]
G2 <- as.matrix(G2)
# G2[1:5,1:5]

m1 <- as.numeric(t(Fst$F_it))

D1 <- Gmatrix(G2, method="VanRaden", ploidy=4, ploidy.correction = T, weights = m1)

m2 <- as.numeric(t(Fst$F_st))
D2 <- Gmatrix(G2, method="VanRaden", ploidy=4, ploidy.correction = T, weights = m2)

m3 <- as.numeric(t(Fst$F_st1))
D3 <- Gmatrix(G2, method="VanRaden", ploidy=4, ploidy.correction = T, weights = m3)

save(D1,D2,D3, file = "WGBLUP_SFT.RData")

