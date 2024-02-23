rm(list = ls())
library(plotly)
library(dplyr)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")

QTL_02 <- QTL_01 %>% filter(Trait != "MET_all")
QTL_02$Chrom <- as.character(QTL_02$Chrom)
QTL_02$Chrom[grepl("^contig*", QTL_02$Chrom)] <- "Contig"
df1 <- count(QTL_02, Model)
df2 <- count(QTL_02, Chrom)
colnames(df2) <- c("Model", "n")

i_6 <- plot_ly(labels = ~Model, values = ~n, legendgroup = ~Model,
               textposition = 'outside',textinfo = 'label+value') %>%
  add_pie(data = df1, name = "DF1", domain = list(x = c(0, 0.3), y = c(0.1, 0.6)))%>%
  add_pie(data = df2, name = "DF2", domain = list(x = c(0.7, 1), y = c(0.1, 0.6))) %>%
  layout(grid=list(rows=1, columns=2), showlegend = FALSE, autosize = F,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


orca(i_6, "pie1.svg")
