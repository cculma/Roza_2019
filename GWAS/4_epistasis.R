rm(list = ls())
# MAEIF1 epistasis 137 markers

setwd("~/Documents/git/big_files/")
# data <- read.csv("epistasis01.csv", row.names = 1)
data <- read.csv("epistasis_model2f.csv", row.names = 1)
# epistasis_model2f.csv is a file generated using a single step MLM using ASReml


plot(density(data$PHENOTYPE))

Single_Marker_Significance <- sapply( setdiff( colnames(data), c("PHENOTYPE")), 
                                      function(snp){
                                        fit <- lm( as.formula( paste("PHENOTYPE~", snp, sep="")),data=data);
                                        t <- anova( fit);                       
                                        return( t$`Pr(>F)`[1])})

Single_Marker_Significance <- unlist(Single_Marker_Significance)

ntree = 20
n = length(Single_Marker_Significance)
q = round(ntree/3,0)

Single_Marker_Significance <- data.frame( pvalue=Single_Marker_Significance, Rank_Markers=rank(Single_Marker_Significance))

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

source("~/Documents/git/epiMEIF/codes/Interactions_Score_Age1.R")

names(weight_variable)<- rownames(Single_Marker_Significance)
weight_variable <- weight_variable[match(colnames(data)[-1], names(weight_variable))]

set.seed(100)
Seeds=sample(1:1000000, 10)
Importance_Score <- list()

dim(data) # 424  21

## Parameters we use for the random forest. This can be trained prior to the forest application.
ntree=400
q=round(ntree/3,0)
#mincriterion=0.95, minsplit = 30, minbucket = 30

addq <- function(x) paste0("`", x, "`")
fit.rf.formula <-  as.formula( paste(" PHENOTYPE ~(", paste( addq(setdiff(colnames(data), c("PHENOTYPE"))), collapse = "+"), ")"))
fit.rf.formula

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
dim(All_Interactions_Stats_SRF) # 185  16


# All1 <- All_Interactions_Stats_SRF[order(All_Interactions_Stats_SRF$Sum_Forest_Score, decreasing = T), ]
# All2 <- All_Interactions_Stats_SRF %>% dplyr::filter(Sum_Forest_Score >= 17)
All2 <- All_Interactions_Stats_SRF %>% dplyr::filter(Median_Forest_Score > 1)
All2 <- All2[order(All2$Sum_Forest_Score, decreasing = T), ]

setwd("~/Documents/git/Roza_2019/epiMEIF/")
write.csv(All_Interactions_Stats_SRF, "logitudinal_model2f.csv", quote = F, row.names = T)
write.csv(All2, "logitudinal_model2f.1.csv", quote = F, row.names = T)


library(viridis)
library(gridExtra)
library(ggpubr)
library(svglite)

data4 <- data %>% mutate_if(is.integer, as.factor)
plot1 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240")) # 6
plot2 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434")) # 3
plot3 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr6.1_27950240")) # 5
plot4 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434","chr6.1_27950240")) # 7
plot5 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_11234756")) # 9
plot6 <- plotSNPInteraction(data4, c("chr6.1_27950240", "chr7.1_11234756")) # 27
plot7 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240","chr7.1_11234756")) # 60
plot8 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr6.1_27950240")) # 18
plot9 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr6.1_27950240")) # 17
plot10 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr5.1_1429434","chr6.1_27950240")) # 19
plot12 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr7.1_11234756")) # 10
plot11 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434","chr7.1_11234756")) # 11
plot13 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr5.1_54813818")) # 22
plot14 <- plotSNPInteraction(data4, c("chr4.1_20948923", "chr6.1_27950240")) # 24


ggboxplot(data4, x = "chr2.1_41048095", y = "PHENOTYPE", width = 0.8, add = "jitter", color = "chr2.1_41048095")
ggboxplot(data4, x = "chr5.1_1429434", y = "PHENOTYPE", width = 0.8, add = "jitter", color = "chr5.1_1429434")
ggboxplot(data4, x = "chr6.1_27950240", y = "PHENOTYPE", width = 0.8, add = "jitter", color = "chr6.1_27950240")


lm1 <- lm(PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434 + chr6.1_27950240, data = data4)
anova(lm1)

lm2 <- lm(PHENOTYPE ~ chr2.1_41048095:chr5.1_1429434:chr6.1_27950240, data = data4)
anova(lm2)

anova(lm1, lm2)


ggline(data4, x = "chr2.1_41048095", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr5.1_1429434", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr6.1_27950240", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr7.1_11234756", y = "PHENOTYPE", add = "mean", linetype = "chr6.1_27950240", shape = "chr6.1_27950240", color = "chr6.1_27950240")

data[1:5,1:5]
data5 <- data %>% dplyr::select(PHENOTYPE, cc05$node)
data5 <- data5 %>% gather(key = "marker", value = "AlleleDos", 2:8)
head(data5)

ggline(data4, x = "AlleleDos", y = "PHENOTYPE", add = "mean")
ggline(data5, x = "AlleleDos", y = "PHENOTYPE",  add = c("mean_se", "jitter"), linetype = "marker", shape = "marker", color = "marker", alpha = 0.5) + facet_wrap(vars(marker))

ggline(data5, x = "AlleleDos", y = "PHENOTYPE",  add = c("mean_se"), linetype = "marker", color = "marker", alpha = 0.5) + facet_grid(, vars(marker))

ggline(data5, x = "AlleleDos", y = "PHENOTYPE",  add = c("mean_se"), linetype = "marker") + facet_grid(, vars(marker)) # this

ggboxplot(data5, x = "AlleleDos", y = "PHENOTYPE", width = 0.8, color = "marker")

plot06 <- ggline(data4, x = "chr5.1_1429434", y = "PHENOTYPE", add = "mean_se", color = "chr6.1_27950240", palette = cols1) + ylab("g x plant-1")

plot03 <- ggline(data4, x = "chr2.1_41048095", y = "PHENOTYPE", add = "mean_se", color = "chr5.1_1429434", palette = cols1) + ylab("g x plant-1")

plot03 <- ggline(data4, x = "chr5.1_1429434", y = "PHENOTYPE", add = "mean_se", color = "chr2.1_41048095", palette = cols1) + ylab("g x plant-1")

levels(data4$chr5.1_1429434)

library(RColorBrewer)
myColors <- brewer.pal(5, "Set1")
cols1 <- c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3", "4" = "#FF7F00" )
cols1 <- c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3", "4" = "#FF7F00" )



str_stack <- function(x) {
  x %>% str_split("") %>% map(~ .x %>% paste(collapse = "\n"))
}


# Interaction_6 <- data %>% dplyr::select(PHENOTYPE, chr5.1_1429434, chr6.1_27950240)
# Interaction_6$all <- paste0(Interaction_6$chr5.1_1429434, Interaction_6$chr6.1_27950240)
# Interaction_6$all <- as.factor(Interaction_6$all)
# ggboxplot(Interaction_6, x = "all", y = "PHENOTYPE", width = 0.8, add = "jitter", color = "all") + scale_x_discrete(label = str_stack) + theme(legend.position="none")

Interaction_7 <- data %>% dplyr::select(PHENOTYPE, chr2.1_41048095, chr5.1_1429434, chr6.1_27950240)
Interaction_7$all <- paste0(Interaction_7$chr2.1_41048095, Interaction_7$chr5.1_1429434, Interaction_7$chr6.1_27950240)
Interaction_7$all <- as.factor(Interaction_7$all)
summary(Interaction_7$all)
plot07 <- ggboxplot(Interaction_7, x = "all", y = "PHENOTYPE", width = 0.8, color = "all") + scale_x_discrete(label = str_stack) + theme(legend.position="none") + ylab("g x plant-1")

# Interaction_60 <- data %>% dplyr::select(PHENOTYPE, chr5.1_1429434, chr6.1_27950240, chr7.1_11234756)
# Interaction_60$all <- paste0(Interaction_60$chr5.1_1429434, Interaction_60$chr6.1_27950240, Interaction_60$chr7.1_11234756)
# Interaction_60$all <- as.factor(Interaction_60$all)
# ggboxplot(Interaction_60, x = "all", y = "PHENOTYPE", width = 0.8, add = "jitter", color = "all") + scale_x_discrete(label = str_stack) + theme(legend.position="none")

library(svglite)

g1 <- ggplotGrob(lev5[[1]])
g2 <- ggplotGrob(lev5[[2]])

setwd("~/Documents/git/Roza_2019/Figures/")
svglite(filename = "Interaction_6.svg", width = 3, height = 2.5, fix_text_size = F)
plot(plot06)
invisible(dev.off())

svglite(filename = "Interaction_3.svg", width = 3, height = 2.5, fix_text_size = F)
plot(plot03)
invisible(dev.off())

svglite(filename = "Interaction_7.svg", width = 5, height = 2.5, fix_text_size = F)
plot(plot07)
invisible(dev.off())

