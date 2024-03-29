# generate manhattan plot by trait

rm(list = ls())
# library(plotly)
# library(dplyr)
# library(ggpubr)
library(GWASpoly)
library(ggplot2)

# load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
load("data_3_80177_year.RData")

load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/myplot1.RData")
load("~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/myplot1.RData")
myplot1 <- LD.plot(data_5) + theme_light(base_family = "Arial", base_size = 12) + ggtitle("LD plot")
myplot2 <- myplot1 + theme_light(base_family = "Arial", base_size = 12) + ggtitle("LD plot")

ggsave(filename = "~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/myplot2.pdf", plot = myplot2, width = 3, height = 3)

ggsave(filename = "myplot1.jpg", plot = myplot1, width = 6, height = 6)
save(myplot1, file = "myplot1.RData")
save(myplot1, file = "myplot1.RData")

load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/myplot1.RData")
myplot2 <- myplot1 + theme_light()

ggsave(filename = "myplot2.pdf", plot = myplot2, width = 3, height = 3)
