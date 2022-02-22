library(plotly)
stocks <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/stockdata2.csv", stringsAsFactors = FALSE)

p <- ggplot(stocks, aes(sample=change)) + geom_qq() + geom_qq_line()

p <- ggplot(stocks, aes(sample=change))+ geom_qq_line() + geom_qq(aes(colour=stock), alpha=0.3, size=0.1) + ylim(-10,10)

p <- ggplot(stocks, aes(sample=change)) + geom_qq_line() + geom_qq(aes(colour=stock), alpha=0.3, size=0.1) + facet_wrap(~stock) + ylim(-10,10)

ggplotly(p)


data_7 <- na.omit(data_6)

p <- ggplot(data_6, aes(sample=may_20_1stage.general)) + geom_qq() + geom_qq_line()

# p <- ggplot(stocks, aes(x=change)) + geom_density(aes(color=stock))

ggplotly(p)


############
# qqplot 
library(car)


data_5 <- set.threshold(data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)

data_6 <- as.data.frame(data_5@scores)
# data_6[is.na(data_6)] <- 0
rownames(data_6) <- paste(data_5@map$Chrom, data_5@map$Position, sep = "_")

# qq.plot(data = data_5, trait = "MET_jul", chrom = "Chr1")
p <- qqPlot(data_6$may_20_1stage.general)
qqPlot(data_6$FA1_all.general)
qqPlot(data_6$FA1_all.additive)

