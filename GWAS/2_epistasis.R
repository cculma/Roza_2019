# MAEIF1

rm(list = ls())

library(GWASpoly)
library(tidyverse)
library(vcfR)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggthemes)
library(hrbrthemes)
library(plotly)
library(GenomicRanges)
library(genomation)
library(plyranges)
library(Repitools)
library(devtools)
library(sommer)
library(GWASpoly)
# run epiMEIF in git 


source("~/Documents/git/epiMEIF/codes/Interactions_Score_Age1.R")
# source("~/Documents/git/epiMEIF/codes/Interactions_Score_Ageing.R")
# source("~/Documents/git/epiMEIF/codes/MEIF.R")
# source("~/Documents/git/epiMEIF/codes/Interaction_Score_RandomForest.R")

setwd("~/Documents/git/big_files/")
load("Yi_st1_51081_sqrt.RData")
data_3 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_03 <- get.QTL(data_3)

cc <- dplyr::count(QTL_03, Trait)
sep21 <- QTL_03 %>% dplyr::filter(Trait == "sep_21")
sep21 <- sep21 %>% distinct(Marker, .keep_all = T)
sep21 <- sep21 %>% unite(col = "marker1", c("Chrom", "Position"), sep = "_", remove = F)

setwd("~/Documents/git/big_files/")
Y <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", header = T)
colnames(Y)
head(Y)
Y <- Y %>% dplyr::select(Roza2019_VCF_ID, sep_21)
Y1 <- na.omit(Y) 

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks) # 424 51081
class(marks)

marks.1 <- marks %>% dplyr::select(sep21$marker1)
dim(marks.1) # 424  48
marks.1[1:5,1:5]

Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2) # 424  19
Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2
# Y3 <- Y2 %>% dplyr::select(sep_21)
data <- merge(as.data.frame(Y3), as.data.frame(marks.1), by = 'row.names', all = TRUE) %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
colnames(data)[1] <- "PHENOTYPE"
data <- na.omit(data)
dim(data) # 424  49 markers
data[1:5,1:5]

plot(density(data$PHENOTYPE))

Single_Marker_Significance <- sapply( setdiff( colnames(data), c("PHENOTYPE")), 
                                      function(snp){
                                        fit <- lm( as.formula( paste("PHENOTYPE~", snp, sep="")),data=data);
                                        t <- anova( fit);                       
                                        return( t$`Pr(>F)`[1])})

Single_Marker_Significance <- unlist(Single_Marker_Significance)

ntree = 48
n = length(Single_Marker_Significance)
q = round(ntree/3,0)

Single_Marker_Significance <- data.frame( pvalue=Single_Marker_Significance, Rank_Markers=rank(Single_Marker_Significance))

hist(-log10(Single_Marker_Significance$pvalue))
hist(Single_Marker_Significance$pvalue)

Single_Marker_Significance <- Single_Marker_Significance[ order( Single_Marker_Significance$Rank_Markers),]

max_ind <- which(Single_Marker_Significance$pvalue > 0.1)[1]

Single_Marker_Significance$weights <- Single_Marker_Significance$weights1 <-  rep(1, nrow(Single_Marker_Significance))  

Single_Marker_Significance$weights[1:(n-q+1)]<- sapply(1:(n-q+1),
                                                       function(i) as.numeric(comboCount(n-i, q-1)/comboCount(n,q)))

Single_Marker_Significance$weights[(max_ind+1):n]<- Single_Marker_Significance$weights[(max_ind)]

Single_Marker_Significance$weights1[1:max_ind] <- (log(1/Single_Marker_Significance$weights[1:(max_ind)])/max(log(1/Single_Marker_Significance$weights[1:(max_ind)])))

Single_Marker_Significance$weights1[(max_ind+1):n] <- Single_Marker_Significance$weights1[(max_ind+1)]

weight_variable <- Single_Marker_Significance$weights1
Single_Marker_Significance$weights <- weight_variable


# plot 1 ------------------------------------------------------------------
library(hrbrthemes)
#Plotting the weights and the p-values
temperatureColor <- "#69b3a2"
priceColor <- rgb(0.2, 0.6, 0.9, 1)

ggplot(Single_Marker_Significance, aes(x=Rank_Markers)) +
  geom_line( aes(y=-log10(pvalue)),linewidth=2, color=temperatureColor) + 
  geom_line( aes(y=weights1*6.57), linewidth=2, color=priceColor) +
  scale_y_continuous(
    # Features of the first axis
    name = "-log10(p-value)", limit=c(0,7),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*(1/6.57), name="Weights")) + 
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = temperatureColor),
    axis.title.y.right = element_text(color = priceColor)) +
  ggtitle("P-values down, weight of sampling up")


# cforests ----------------------------------------------------------------
# restart session

names(weight_variable)<- rownames(Single_Marker_Significance)
weight_variable <- weight_variable[match(colnames(data)[-1], names(weight_variable))]

set.seed(100)
Seeds=sample(1:1000000, 10)
Importance_Score <- list()

## Parameters we use for the random forest. This can be trained prior to the forest application.
#ntree=50
#q=round(ntree/3,0)
#mincriterion=0.95, minsplit = 30, minbucket = 30

addq <- function(x) paste0("`", x, "`")
fit.rf.formula <-  as.formula( paste(" PHENOTYPE ~(", paste( addq(setdiff(colnames(data), c("PHENOTYPE"))), collapse = "+"), ")"))

for(seed in Seeds)
{
  Importance_Score[[match(seed,Seeds)]] <- list()
  set.seed(seed)
  fit.rf <-cforest_gen( fit.rf.formula, data = data, weight_variable=c(0,weight_variable), mtry = round(ntree/3,0), ntree = ntree, control = ctree_control(teststat="quad", testtype="Univ", mincriterion=0.95, minsplit = 30, minbucket = 30)) 
  
  
  ##computing the SNP interaction metric for each forest run using the getInteractionMatrix function.
  Importance_Score[[match(seed,Seeds)]] <-getInteractionMatrix(fit.rf)
}

NiterChild=10
Importance_Score_SRF <- Interaction_List_SRF <- list()
Interaction_List_RF <- list()

##Recording the binary interaction score 
for(i in 1:NiterChild)
  Importance_Score_SRF[[i]] <- Importance_Score[[i]][[1]]

##Recording the interaction sets and their corresponding score
for(i in 1:NiterChild)
{if(nrow(Importance_Score[[i]][[2]])!=0)
  Interaction_List_SRF[[i]] <- Importance_Score[[i]][[2]][,-match("Forest_Score", colnames(Importance_Score[[i]][[2]]))]
}

##Summarizing over the SNP Interaction Metrix across the 10 forest using the GenerateInteractionList.
All_Interactions_Stats_SRF <- GenerateInteractionList(Interaction_List_SRF, Importance_Score)


library(viridis)
library(gridExtra)

plotSNPInteraction(data, c("chr5.1_42131042", "chr6.1_49779667"))

data4 <- data %>% mutate_if(is.integer, as.factor)
plotSNPInteraction(data4, c("chr5.1_42131042", "chr6.1_49779667"))
plotSNPInteraction(data2, c("chr5.1_42131042", "chr6.1_67638734"))

