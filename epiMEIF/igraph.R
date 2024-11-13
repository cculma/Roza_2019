# install.packages("igraph")
# install.packages("igraphdata")
# install.packages('qgraph')

library(igraph)
library(igraphdata)
require(qgraph)
qgraph(g)
qgraph(m,edge.labels=TRUE)  #

g1 <- make_graph( edges=c(1,2, 2,3, 3,1), n=3, directed=F ) 
plot(g1) 

All1 <- All_Interactions_Stats_SRF[order(All_Interactions_Stats_SRF$Sum_Forest_Score, decreasing = T), ]
 
cc01 <- dplyr::count(All_Interactions_Stats_SRF, Node1)
cc02 <- dplyr::count(All_Interactions_Stats_SRF, Node2)
cc03 <- dplyr::count(All_Interactions_Stats_SRF, Node3)
cc04 <- dplyr::count(All_Interactions_Stats_SRF, Node4)
colnames(cc01)[1] <- "node"
colnames(cc02)[1] <- "node"
colnames(cc03)[1] <- "node"
colnames(cc04)[1] <- "node"

cc05 <- rbind(cc01, cc02, cc03,cc04)


g1 <- make_graph( edges=c(All1[1,1], All1[1,2], All1[2,1], All1[2,2],All1[2,1], All1[2,2]), n=6, directed=T ) 
plot(g1)



## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(g)
## The opposite operation
as_data_frame(g, what="vertices")
as_data_frame(g, what="edges")

hist(All1$Sum_Forest_Score)
hist(All1$Median_Forest_Score)

All2 <- All1 %>% dplyr::filter(Median_Forest_Score > 1)
g <- graph_from_data_frame(All2, directed = T)
E(g)$weight <- All2$Median_Forest_Score


net2 <- simplify(g, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum"))
net2 <- simplify(g, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum", "ignore"))


plot( net2, layout = layout.reingold.tilford,
      edge.width=(E(g)$weight) * 0.4,
      edge.arrow.width = 0.08,
      vertex.size = 5,
      edge.arrow.size = 1,
      vertex.size2 = 1,
      vertex.label.cex = 0.5,
      asp = 0.35,
      margin = 0.2)

plot01 <- plot(g, layout = layout.reingold.tilford,
      edge.width= log(E(g)$weight),
      edge.arrow.width = 0.5,
      vertex.size = 5,
      edge.arrow.size = 0.5,
      vertex.size2 = 3,
      vertex.label.cex = 1,
      asp = 0.35,
      margin = 0.2)


library(ggplot2)
setwd("~/Documents/git/Roza_2019/epiMEIF/")
ggsave("fig01.png",
        plot(g, layout = layout.reingold.tilford,
      edge.width= log(E(g)$weight),
      edge.arrow.width = 0.5,
      vertex.size = 5,
      edge.arrow.size = 0.5,
      vertex.size2 = 3,
      vertex.label.cex = 1,
      asp = 0.35,
      margin = 0.2),
        height = 5, width = 5, dpi = 1000)


class(plot01)

plot( net2, layout = layout_as_star,
      edge.width=(E(g)$weight) * 0.4,
      edge.arrow.width = 0.5,
      vertex.size = 10,
      edge.arrow.size = 0.3,
      vertex.size2 = 0.5,
      vertex.label.cex = 1,
      asp = 0.8,
      margin = 0.2)

plot( net2, layout = layout_with_sugiyama,
      edge.width=(E(g)$weight) * 0.025,
      edge.arrow.width = 0.5,
      vertex.size = 5,
      edge.arrow.size = 0.3,
      vertex.size2 = 0.5,
      vertex.label.cex = 1,
      asp = 1,
      margin = 0.2)



lay2 <- layout_with_sugiyama(g, attributes="all",  
                              hgap=1, vgap=1) 

log(E(g)$weight)

plot(lay2$extd_graph, edge.width= log(E(g)$weight),
     edge.arrow.width = 0.5,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     vertex.size2 = 0.5,
     vertex.label.cex = 1,
     asp = 0.5,
     margin = 0.2)



lay2 <- layout_with_sugiyama(net2, attributes="all",  
                             hgap=1, vgap=1) 

plot(lay2$extd_graph, edge.width= log(E(g)$weight),
     edge.arrow.width = 0.5,
     vertex.size = 5,
     edge.arrow.size = 0.3,
     vertex.size2 = 0.5,
     vertex.label.cex = 1,
     asp = 0.5,
     margin = 0.2)
