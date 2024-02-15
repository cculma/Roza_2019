
setwd("~/Documents/git/big_files/")
rm(list = ls())

library(VariantAnnotation)
library(vcfR)
library(updog)
library(future)


vcf1 <- read.vcfR("Roza2019_04.vcf")
vcf2 <- readVcf("Roza2019_04.vcf")

any(elementNROWS(rowRanges(vcf2)$ALT) > 1)

BSDP <- geno(vcf2)$BSDP
dim(BSDP) # 62839   499     4
dimnames(BSDP)[[3]] <- c("A","C","G","T")

altvec <- unlist(CharacterList(rowRanges(vcf2)$ALT))
refvec <- as.character(rowRanges(vcf2)$REF)
stopifnot(length(altvec) == length(refvec))

refmat <- matrix(NA_real_, nrow = dim(BSDP)[[1]], ncol = dim(BSDP)[[2]])
altmat <- matrix(NA_real_, nrow = dim(BSDP)[[1]], ncol = dim(BSDP)[[2]])
dimnames(refmat) <- dimnames(BSDP)[1:2]
dimnames(altmat) <- dimnames(BSDP)[1:2]

for (i in seq_len(nrow(refmat))) {
  refmat[i, ] <- BSDP[i, , refvec[[i]]]
  altmat[i, ] <- BSDP[i, , altvec[[i]]]
}

sizemat <- refmat + altmat

mout <- multidog(refmat = refmat,
                 sizemat = sizemat,
                 ploidy = 4,
                 model = "norm",
                 nc = 50)

setwd("~/Documents/git/big_files/")
save(mout, file = "mout_Roza2019.RData")

vcf <- vcf1
ADmat <- extract.gt(vcf, "AD")
refmat <- strsplit(ADmat, ",")
refmat <- sapply(refmat, "[[", 1)
refmat <- matrix(refmat, nrow=nrow(ADmat))
refmat <- apply(refmat, 2, as.numeric)
altmat <- strsplit(ADmat, ",")
altmat <- sapply(altmat, "[[", 2)
altmat <- matrix(altmat, nrow=nrow(ADmat))
altmat <- apply(altmat, 2, as.numeric)
colnames(refmat) <- colnames(ADmat)
colnames(altmat) <- colnames(ADmat)
rownames(refmat) <- rownames(ADmat)
rownames(altmat) <- rownames(ADmat)
sizemat <- refmat + altmat
print(refmat[1:5,1:5])
print(sizemat[1:5,1:5])
mout_UMN3097 <- multidog(refmat = refmat,
                         sizemat = sizemat,
                         ploidy = 4,
                         model = "norm",
                         nc = 60)
save(mout_UMN3097, file = "mout_UMN3097.Rdata")
