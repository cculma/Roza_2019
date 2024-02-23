



# Load the required packages
library(Rtsne)
library(ggplot2)


library(heplots)
library(candisc)
library(car)
library(PopGenReport)
library(dartR)
library(adegenet)
library(hierfstat)
library(tidyverse)
library(RColorBrewer)
library(janitor)
library(caret)
library(ASRgenomics)





dim(G4)
class(G4)

?Rtsne
# perform dimensionality redcution from 64D to 2D
tsne <- Rtsne(as.matrix(trn[,1:64]), check_duplicates = FALSE, pca = FALSE, perplexity=30, theta=0.5, dims=2)

tsne <- Rtsne(G4, check_duplicates = T, pca = T, perplexity=30, theta=0.5, dims=2)
# display the results of t-SNE
cols <- rainbow(30)

# Conversion of matrix to dataframe
tsne_plot <- data.frame(x = tsne$Y[,1], 
                        y = tsne$Y[,2])

Y3 <- inner_join(Y1, Y2, by = "Sample_ID")
Y3 <- cbind(Y3, tsne_plot)
head(Y3)

# Plotting the plot using ggplot() function
ggplot(tsne_plot) + geom_point(aes(x=x,y=y))

ggplot(Y3, aes(x = x, y = y, color = Susceptible_parent)) + geom_point() + theme_minimal() + ggrepel::geom_text_repel(aes(label = Susceptible_parent, size = 4), nudge_x = .75)


# DAPC

b2 <- G4
rownames(b2)
b2[1:5,1:5]
dim(b2) # 424 51081

a2.5 <- as.list(as.data.frame(t(b2)))

Y3 <- inner_join(Y1, Y2, by = "Sample_ID")
head(Y3)
Y3 <- Y3 %>% column_to_rownames("Roza2019_VCF_ID")

order1 <- match(rownames(G4), rownames(Y3))
Y3  <- Y3[order1,]

head(rownames(G4))
head(rownames(Y3))

ind <- as.character(rownames(Y3))
population <- as.character(Y3$Susceptible_parent)

a1.7 <- as.data.frame(colnames(b2))
colnames(a1.7) <- "Marker"
a1.7 <- a1.7 %>% separate(col = Marker, into = c("Chrom","ChromPos"), sep = "_", remove = F)

x <- new("genlight", a2.5, ploidy = 4, ind.names = ind, pop = population, chromosome = a1.7$Chrom1, position = a1.7$ChromPos) 
x

class(x)
pop(x)

grp1 <- find.clusters(x, choose.n.clust = T, method = "ward", criterion = "diffNgroup", scale = F, pca.select = "percVar", perc.pca = 80, max.n.clust = 10, stat = "BIC") 

dapc1 <- dapc(x = x, grp1$grp, var.loadings = T, var.contrib = T, pca.select = "percVar", perc.pca = 80, n.da=5)

# scatter
my.cols <- brewer.pal(5, "Set1")
my.cols <- c("1"= "#FF7F00", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#E41A1C", "5" = "#984EA3")
col <- rainbow(length(levels(pop(x))))

sc1 <- scatter(dapc1,
               bg = "white", solid = 1, cex = 1, # cex circle size
               col = my.cols,
               pch = 20, # shapes
               cstar = 1, # 0 or 1, arrows from center of cluster
               cell = 2, # size of elipse
               scree.da = T, # plot da
               scree.pca = T, # plot pca
               posi.da = "topright", 
               posi.pca="bottomright",
               mstree = F, # lines connecting clusters
               lwd = 1, lty = 2, 
               leg = F, clab = 1) # legend and label of legend clusters. clab 0 or 1


?find.clusters
grp2 <- find.clusters(x, choose.n.clust = T, method = "ward", criterion = "goodfit", scale = F, pca.select = "percVar", perc.pca = 80, max.n.clust = 12, stat = "BIC") 

dapc2 <- dapc(x = x, grp2$grp, var.loadings = T, var.contrib = T, pca.select = "percVar", perc.pca = 80, n.da=3)

grp3 <- find.clusters(x, choose.n.clust = T, max.n.clust = 12, stat = "BIC")
grp3$grp
# Choose the number PCs to retain (>=1): 
#   300
# Choose the number of clusters (>=2): 
#   3
dapc3 <- dapc(x, grp3$grp, n.pca= NULL, n.da= NULL, var.contrib = T, scale = F)

my.cols <- brewer.pal(3, "Set1")
my.cols


sc2 <- scatter(dapc2,
               bg = "white", solid = 1, cex = 1, # cex circle size
               col = my.cols,
               pch = 20, # shapes
               cstar = 1, # 0 or 1, arrows from center of cluster
               cell = 2, # size of elipse
               scree.da = T, # plot da
               scree.pca = T, # plot pca
               posi.da = "topright", 
               posi.pca="bottomright",
               mstree = F, # lines connecting clusters
               lwd = 1, lty = 2, 
               leg = F, clab = 1) # legend and label of legend clusters. clab 0 or 1

Y8 <- as.data.frame(grp2$grp)
Y8 <- as.data.frame(grp3$grp)

colnames(Y8) <- "DAPC"
Y9  <- merge(Y3, Y8, by = "row.names")
colnames(Y9)[1] <- "Roza2019_VCF_ID"

setwd("~/Documents/git/Roza_2019/pheno_data/")
write.csv(Y9, "PCA_Roza2019.csv", row.names = F, quote = F)


cc <- dplyr::count(Y9, DAPC, Susceptible_parent1)
cc <- cc %>% spread(key = DAPC, value = n)
colnames(cc)[2:4] <- c("DAPC1", "DAPC2","DAPC3")
cc1 <- dplyr::count(Y9, kmeans, Susceptible_parent1)
cc1 <- cc1 %>% spread(key = kmeans, value = n)
colnames(cc1)[2:4] <- c("kmeans1", "kmeans2","kmeans3")
cc2 <- inner_join(cc, cc1, by = "Susceptible_parent1")

setwd("~/Documents/git/Roza_2019/pheno_data/")
write.csv(cc2, "DAPC_groups.csv", row.names = F, quote = F)


# set.seed(999)
# pramx <- xvalDapc(tab(x, NA.method = "mean"), pop(x),
#                               n.pca = 200:300, n.rep = 1000,
#                               parallel = "multicore", ncpus = 40)

