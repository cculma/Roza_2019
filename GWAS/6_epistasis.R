rm(list = ls())
# MAEIF1 epistasis sep_21 markers

setwd("~/Documents/git/big_files/")
data <- read.csv("epistasis_sep_21.csv", row.names = 1)


plot(density(data$PHENOTYPE))

Single_Marker_Significance <- sapply( setdiff( colnames(data), c("PHENOTYPE")), 
                                      function(snp){
                                        fit <- lm( as.formula( paste("PHENOTYPE~", snp, sep="")),data=data);
                                        t <- anova( fit);                       
                                        return( t$`Pr(>F)`[1])})

Single_Marker_Significance <- unlist(Single_Marker_Significance)

# ntree = 300
# n = length(Single_Marker_Significance)
# q = round(ntree/3,0)

ntree = 48
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

## Parameters we use for the random forest. This can be trained prior to the forest application.
ntree=300
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
dim(All_Interactions_Stats_SRF) # 2873   18

# setwd("~/Documents/git/Roza_2019/epiMEIF/")
# write.csv(All_Interactions_Stats_SRF, "logitudinal_sep_21.csv", quote = F, row.names = T)


library(viridis)
library(gridExtra)
library(ggpubr)
library(svglite)

data4 <- data %>% mutate_if(is.integer, as.factor)
plot01_sep21 <- plotSNPInteraction(data4, c("chr3.1_85101832", "chr6.1_39491281"))
# PDF 7 X 4

plot1_sep21.1 <- ggline(data4, x = "chr6.1_49779667", y = "PHENOTYPE", add = "mean", linetype = "chr5.1_42131042", shape = "chr5.1_42131042", color = "chr5.1_42131042")
plot1_sep21.1


plot02_sep21 <- plotSNPInteraction(data4, c("chr6.1_27950240", "chr7.1_11234756"))
plot02_sep21.1 <- ggline(data4, x = "chr5.1_42131042", y = "PHENOTYPE", add = "mean", linetype = "chr6.1_67638734", shape = "chr6.1_67638734", color = "chr6.1_67638734")

library(svglite)
setwd("~/Documents/git/Roza_2019/Figures/")
svglite(filename = "plot02_sep21.1.svg", width = 4, height = 3, fix_text_size = F)
plot(plot02_sep21.1)
invisible(dev.off())

ggline(data4, x = "chr7.1_11234756", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr6.1_27950240", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr7.1_11234756", y = "PHENOTYPE", add = "mean", linetype = "chr6.1_27950240", shape = "chr6.1_27950240", color = "chr6.1_27950240")


"chr3.1_85101832"
"chr6.1_39491281"


ggline(data4, x = "chr3.1_85101832", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr6.1_39491281", y = "PHENOTYPE", add = "mean")
ggline(data4, x = "chr3.1_85101832", y = "PHENOTYPE", add = "mean", linetype = "chr6.1_39491281", shape = "chr6.1_39491281", color = "chr6.1_39491281")

ggscatter(data, x = "chr7.1_11234756", y = "PHENOTYPE",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "loess",
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))


plot3 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_24014763"))