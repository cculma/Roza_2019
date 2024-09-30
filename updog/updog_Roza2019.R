
setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/BI_Project_11/")
setwd("~/Documents/git/big_files/")
rm(list = ls())

library(VariantAnnotation)
library(vcfR)
library(updog)
library(future)
library(ldsep)

vcf2 <- readVcf("Roza_2019.vcf")

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
# save(mout, file = "mout_Roza2019.RData")
load("mout_Roza2019.RData")

# ldsep -------------------------------------------------------------------

nrow(mout$snpdf) # 62839
mout$snpdf[1:5,1:5]
class(mout)
hist(mout$snpdf$prop_mis)
summary(mout[["snpdf"]][["prop_mis"]]) # 0.16575
summary(mout[["snpdf"]][["od"]]) # 0.008906
summary(mout[["snpdf"]][["bias"]]) # 1.11057
msub <- filter_snp(mout, prop_mis < 0.16575 & od < 0.008906 & bias < 1.11057) 

msub$snpdf[1:5,1:5]
dim(msub$snpdf) # 39611    20

ploidy <- 4
# gp <- format_multidog(x = msub, varname = paste0("Pr_", 0:ploidy))
# class(gp)
# dim(gp) # 39611   499     5
# ldout <- ldfast(gp = gp, type = "r2")

setwd("~/Documents/git/big_files/")
# save(ldout, file = "ldout_Roza2019.RData")
load("ldout_Roza2019.RData")

dim(ldout$ldmat)
ldout$ldmat[1:5,1:5]

lev1 <- msub$snpdf$snp
length(lev1)
lev1[1:5]

lev2 <- as.data.frame(lev1)
lev2 <- lev2 %>% separate(1, c("Chr", "position"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev2$Chr <- as.factor(lev2$Chr)
levels(lev2$Chr)

# lev2$marker <- gsub(":", "_", lev2$marker)
# lev2 <- lev2 %>% separate(2, c("Chr", "position"), sep = "_", remove = F, convert = FALSE, extra = "merge")

C1 <- lev2 %>% dplyr::filter(Chr == "chr1.1")
C2 <- lev2 %>% dplyr::filter(Chr == "chr2.1")
C3 <- lev2 %>% dplyr::filter(Chr == "chr3.1")
C4 <- lev2 %>% dplyr::filter(Chr == "chr4.1")
C5 <- lev2 %>% dplyr::filter(Chr == "chr5.1")
C6 <- lev2 %>% dplyr::filter(Chr == "chr6.1")
C7 <- lev2 %>% dplyr::filter(Chr == "chr7.1")
C8 <- lev2 %>% dplyr::filter(Chr == "chr8.1")


R2 <- ldout$ldmat
dim(R2)
length(lev2$lev1)
colnames(R2) <- rownames(R2) <- lev2$lev1

R2[1:5,1:5]
R2.1 <- R2[C1$lev1, C1$lev1]
R2.2 <- R2[C2$lev1, C2$lev1]
R2.3 <- R2[C3$lev1, C3$lev1]
R2.4 <- R2[C4$lev1, C4$lev1]
R2.5 <- R2[C5$lev1, C5$lev1]
R2.6 <- R2[C6$lev1, C6$lev1]
R2.7 <- R2[C7$lev1, C7$lev1]
R2.8 <- R2[C8$lev1, C8$lev1]

R2.1[upper.tri(R2.1)] <- NA
diag(R2.1) <- NA
R2.1 <- reshape2::melt(R2.1, na.rm = T)

R2.2[upper.tri(R2.2)] <- NA
diag(R2.2) <- NA
R2.2 <- reshape2::melt(R2.2, na.rm = T)

R2.3[upper.tri(R2.3)] <- NA
diag(R2.3) <- NA
R2.3 <- reshape2::melt(R2.3, na.rm = T)

R2.4[upper.tri(R2.4)] <- NA
diag(R2.4) <- NA
R2.4 <- reshape2::melt(R2.4, na.rm = T)

R2.5[upper.tri(R2.5)] <- NA
diag(R2.5) <- NA
R2.5 <- reshape2::melt(R2.5, na.rm = T)

R2.6[upper.tri(R2.6)] <- NA
diag(R2.6) <- NA
R2.6 <- reshape2::melt(R2.6, na.rm = T)

R2.7[upper.tri(R2.7)] <- NA
diag(R2.7) <- NA
R2.7 <- reshape2::melt(R2.7, na.rm = T)

R2.8[upper.tri(R2.8)] <- NA
diag(R2.8) <- NA
R2.8 <- reshape2::melt(R2.8, na.rm = T)

R3 <- rbind(R2.1,R2.2,R2.3,R2.4,R2.5,R2.6,R2.7,R2.8)
head(R3)

R3 <- R3 %>% separate(2, c("Locus2", "Position2"), sep = "_", remove = T, convert = FALSE, extra = "merge") %>% separate(1, c("Locus1", "Position1"), sep = "_", remove = T, convert = FALSE, extra = "merge")

R3$Position1 <- as.numeric(R3$Position1)
R3$Position2 <- as.numeric(R3$Position2)
R3$dist <- as.numeric(R3$Position1 - R3$Position2)
str(R3)
dim(R3)
colnames(R3)[5] <- "rsq"
R3 <- R3[R3$dist != "NaN",]


N = 1502 * 2
Cstart <- c(C=0.1)
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), data=R3, start=Cstart, control=nls.control(maxiter=100))

rho <- summary(modelC)$parameters[1]

newrsq <- ( (10+rho*R3$dist) / ( (2+rho*R3$dist) * (11+rho*R3$dist) ) ) *
  ( 1 + ( (3+rho * R3$dist) * (12+12*rho*R3$dist + (rho*R3$dist)^2) ) / 
      (2*N*(2+rho*R3$dist) * (11+rho*R3$dist) ) )

newfile <- data.frame(R3$dist, newrsq)
head(newfile)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$R3.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$R3.dist),]
newfile1 <- newfile


R3$dist
R3$rsq

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")

pdf("LD_DAI1.pdf", height=5, width = 5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(R3$dist, R3$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$R3.dist, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizontal line
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
dev.off()

head(newfile)
newfile$Dist_kb <- newfile$R3.dist/1000
head(R3)

ggplot(newfile, aes(x = R3.dist, y = newrsq)) + geom_line() + xlim(0, 1000000) + theme_bw()

# + geom_hline(yintercept=0.1, linetype="dashed", color = "gray") + geom_vline(xintercept = halfdecaydist)

plot1 <- ggplot(newfile, aes(x = Dist_kb, y = newrsq)) + geom_line() + ylim(0, 0.2) + scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, by = 200)) + theme_classic(base_family = "Arial", base_size = 12) + geom_hline(yintercept=0.1, linetype="dashed", color = "gray") + geom_vline(xintercept = (halfdecaydist/1000), linetype="dashed", color = "blue") + labs(y = expression(LD ~ (r^2)), x = "Distance (kb)")

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/7_LD/")
ggsave(filename = "LD_DAI1.jpg", plot = plot1, width = 4, height = 4)

ggsave(filename = "LD_DAI1.pdf", plot = plot1, width = 4, height = 4, device = cairo_pdf)



lev3 <- mout_2$inddf$ind
