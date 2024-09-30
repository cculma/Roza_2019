# count GWAS markers

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
library(ggvenn)
library(svglite)

# load("~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")
data_3 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_03 <- get.QTL(data_3)

# load("~/Documents/git/big_files/Yi_st2_51081_sqrt.RData")
data_4 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)

# load("~/Documents/git/big_files/data_Yi_DS_20.RData")
data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_05 <- get.QTL(data_5)

QTL_03.1 <- QTL_03 %>% distinct(Marker, .keep_all = T)
QTL_04.1 <- QTL_04 %>% distinct(Marker, .keep_all = T)
QTL_05.1 <- QTL_05 %>% distinct(Marker, .keep_all = T)

x <- list(ST1 = QTL_03.1$Marker,
     ST2 = QTL_04.1$Marker,
     DS = QTL_05.1$Marker)


plot01 <- ggvenn(x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  fill_alpha = 0.3, text_size = 6,
  show_percentage = F,
  stroke_size = 0.5, set_name_size = 6)

plot01
class(plot01)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Venn1.svg", width = 3, height = 4, fix_text_size = F)
plot(plot01)
invisible(dev.off())


QTL_06 <- rbind(QTL_03, QTL_04, QTL_05)
QTL_06.1 <- QTL_06 %>% distinct(Marker, .keep_all = T)


QTL_06.1$Chrom <- as.factor(QTL_06.1$Chrom)
levels(QTL_06.1$Chrom)
QTL_06.1$Chrom <- recode_factor(QTL_06.1$Chrom,
                            chr1.1 = "Chr1", chr2.1 = "Chr2", chr3.1 = "Chr3",
                            chr4.1 = "Chr4", chr5.1 = "Chr5", chr6.1 = "Chr6",
                            chr7.1 = "Chr7", chr8.1 = "Chr8")


df2 <- count(QTL_06.1, Chrom)
colnames(df2) <- c("Model", "n")



# Plotting the piechart using plot_ly() function
i_6 <- plotly::plot_ly(data=df2, values=~n, labels=~ Model, sort = F,
                type="pie", textposition = 'outside',textinfo = 'label+value') %>%
  layout(legend = list(font = list(size = 20)),
         hoverlabel = list(font = list(size = 20)))

class(i_6)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
orca(i_6, "pie_marker_chrom.svg", width = 1, height = 1)



# Create a dataframe for Plot data
data <- data.frame(year = c(2011, 2012, 2013, 2014, 2015),
                   point = c(10, 20, 30, 40, 50))

# Plot the scatter plot with horizontal 
# line at Y=20
ggplot(data, aes(year, point)) +    
  ggtitle("Point Chart")+
  geom_point(aes(size = 1.0), col = "green")+
  geom_hline(yintercept = 20)

setwd("~/Documents/git/big_files/")
GS1 <- read.table("Roza2019_06_GS.txt")
GS2 <- read.csv("Roza2019_06_GWASPoly.txt")

GS1[1:5,1:5]
GS2[1:5,1:5]
GS2 <- GS2[,c(1:3)]
str(GS2)

library(CMplot)
CMplot(GS2,type="h",plot.type="d", bin.size=1e6,file.output=F,verbose=TRUE, chr.border = F)

cc2 <- dplyr::count(GS2, Chrom)
sum(cc2$n)
GS2 %>% group_by(Chrom) %>% summarise_at(Position, c(min = min))


GS2.1 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr1")
GS2.2 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr2")
GS2.3 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr3")
GS2.4 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr4")
GS2.5 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr5")
GS2.6 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr6")
GS2.7 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr7")
GS2.8 <- QTL_06.1 %>% dplyr::filter(Chrom == "Chr8")


#1
chr1.1 = 82459472
chr2.1 = 76462061
chr3.1 = 93164418
chr4.1 = 90245664
chr5.1 = 81211777
chr6.1 = 80303593
chr7.1 = 88407277
chr8.1 = 87242343


markers1 <- ggplot()+
  geom_segment(aes(x=1,xend=chr1.1,y=8,yend=8)) + 
  geom_segment(aes(x=1,xend=chr2.1,y=7,yend=7)) +
  geom_segment(aes(x=1,xend=chr3.1,y=6,yend=6)) +
  geom_segment(aes(x=1,xend=chr4.1,y=5,yend=5)) +
  geom_segment(aes(x=1,xend=chr5.1,y=4,yend=4)) +
  geom_segment(aes(x=1,xend=chr6.1,y=3,yend=3)) +
  geom_segment(aes(x=1,xend=chr7.1,y=2,yend=2)) +
  geom_segment(aes(x=1,xend=chr8.1,y=1,yend=1)) +
  annotate("pointrange", x = GS2.1$Position, y = 8, ymin = 7.8, ymax = 8.2,
             colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.2$Position, y = 7, ymin = 6.8, ymax = 7.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.3$Position, y = 6, ymin = 5.8, ymax = 6.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.4$Position, y = 5, ymin = 4.8, ymax = 5.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.5$Position, y = 4, ymin = 3.8, ymax = 4.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.6$Position, y = 3, ymin = 2.8, ymax = 3.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.7$Position, y = 2, ymin = 1.8, ymax = 2.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + 
  annotate("pointrange", x = GS2.8$Position, y = 1, ymin = 0.8, ymax = 1.2,
           colour = "blue", size = 0.1, linewidth = 0.5, alpha = 0.5) + theme_classic()

markers1

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Markers_pos.svg", width = 5, height = 3, fix_text_size = F)
plot(markers1)
invisible(dev.off())

df2 <- count(QTL_06.1, Chrom)


ggplot(data=df2, aes(x=Chrom, y=n)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=n), vjust=-0.3, size=3.5) +
  theme_classic()


df2 <- as.data.frame(t(data.frame(coding = 85, no_coding = 46, promoter = 6))) %>% rownames_to_column("annot")


# Plotting the piechart using plot_ly() function
i_6 <- plotly::plot_ly(data=df2, values=~V1, labels=~ annot, sort = F,
                       type="pie", textposition = 'outside',textinfo = 'label+value') %>%
  layout(legend = list(font = list(size = 20)),
         hoverlabel = list(font = list(size = 20)))
i_6
class(i_6)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
orca(i_6, "pie_marker_annot.svg", width = 1, height = 1)
