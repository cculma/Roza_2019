
rm(list = ls())
# MAEIF1 in logitudinal dataset X23

setwd("~/Documents/git/big_files/")
data <- read.csv("epistasis02.csv", row.names = 1)
data[1:5,1:5]
data$Stress <- as.factor(data$Stress)
ggplot(data, aes(x = PHENOTYPE, fill = Stress)) + geom_density(alpha = 0.9) + theme_classic()

plot(density(data$PHENOTYPE))

source("~/Documents/git/epiMEIF/codes/Interactions_Score_Age1.R")
##Parameters we use for the random forest
ntree=300
q=round(ntree/3,0)
mincriterion=0.95
minsplit = 30
minbucket = 30

addq <- function(x) paste0("`", x, "`")
fit.rf.formula <-  as.formula( paste(" PHENOTYPE ~(", paste( addq(setdiff(colnames(data), c("PHENOTYPE"))), collapse = "+"), ")"))
fit.rf.formula

set.seed(100)
Seeds=sample(1:1000000, 10)
Importance_Score <- list()

for(seed in Seeds)
{
  Importance_Score[[match(seed,Seeds)]] <- list()
  set.seed(seed)
  
  ##Running the epiforest
  fit.rf <- partykit::cforest( fit.rf.formula, data = data, mtry = round(ntree/3,0), ntree = ntree, control = ctree_control(teststat="quad", testtype="Univ", mincriterion=0.95, minsplit = 30, minbucket = 30)) 
  
  
  ##computing the SNP interaction metric for each forest run using the getInteractionMatrix function.
  Importance_Score[[match(seed,Seeds)]] <- getInteractionMatrix(fit.rf)
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
dim(All_Interactions_Stats_SRF) # 7904   19

# setwd("~/Documents/git/Roza_2019/epiMEIF/")
# write.csv(All_Interactions_Stats_SRF, "logitudinal_23.csv", quote = F, row.names = T)


library(viridis)
library(gridExtra)
str(data)
data4 <- data %>% mutate_if(is.integer, as.factor)
plot1 <- plotSNPInteraction(data4, c("Stress", "chr6.1_71935085"))
plot2 <- plotSNPInteraction(data4, c("Stress", "chr4.1_79960120", "chr6.1_71935085"))
plot3 <- plotSNPInteraction(data4, c("chr4.1_60953306", "chr8.1_84752525"))
plot4 <- plotSNPInteraction(data4, c("Stress", "chr5.1_1429434"))


# PDF 600 X 400
"chr6.1_71935085"
"chr4.1_79960120"
class(plot1)

colnames(data4)

ggline(data4, x = "chr4.1_79960120", y = "PHENOTYPE", add = "mean", linetype = "chr6.1_71935085", shape = "chr6.1_71935085")

ggline(data4, x = "chr6.1_71935085", y = "PHENOTYPE", add = "mean_se", linetype = "Stress", shape = "Stress")

plot5 <- ggline(data4, x = "chr4.1_79960120", y = "PHENOTYPE", add = "mean_se", linetype = "chr6.1_71935085", shape = "chr6.1_71935085") +
  stat_compare_means() +
  stat_compare_means(aes(group = chr6.1_71935085), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5))  + ylim(1, 5)

data4[1:5,1:5]

ggline(data4, x = "chr4.1_79960120", y = "PHENOTYPE", add = "mean_se", linetype = "chr6.1_71935085", shape = "chr6.1_71935085") + facet_wrap(~ Stress, scales = "free")


ggline(data4, x = "chr6.1_71935085", y = "PHENOTYPE", add = "mean_se", linetype = "chr4.1_79960120", shape = "chr4.1_79960120") + facet_wrap(~ Stress, scales = "free")


svglite(filename = "Epistasis2.svg", width = 6, height = 2.5, fix_text_size = F)
plot(plot5)
invisible(dev.off())