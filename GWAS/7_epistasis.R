load("~/Documents/git/big_files/Yi_st2_51081_sqrt.RData")
data_4 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)
cc01 <- dplyr::count(QTL_04, Trait)

QTL_04.1 <- QTL_04 %>% dplyr::filter(Trait == "yi_hrv")
QTL_06.1 <- QTL_04.1 %>% distinct(Marker, .keep_all = T)

QTL_06.1 <- QTL_06.1 %>% unite(col = "marker1", c("Chrom", "Position"), sep = "_", remove = F)

setwd("~/Documents/git/big_files/")
Y <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", header = T)
colnames(Y)
head(Y)
Y <- Y %>% dplyr::select(Roza2019_VCF_ID, yi_hrv)
Y1 <- na.omit(Y) 
Y1$yi_hrv <- (Y1$yi_hrv * 453.592)
head(Y1)

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks) # 424 51081
class(marks)

marks.1 <- marks %>% dplyr::select(QTL_06.1$marker1)
dim(marks.1) # 424  11

Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2) # 424  2
Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2

data <- merge(as.data.frame(Y3), as.data.frame(marks.1), by = 'row.names', all = TRUE) %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
colnames(data)[1] <- "PHENOTYPE"
data <- na.omit(data)
dim(data) # 424  138 : 424 genotypes and 137 markers

setwd("~/Documents/git/big_files/")
write.csv(data, "epistasis_yi_hrv.csv", quote = F, row.names = T)
