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

All2
 

cc01 <- dplyr::count(All2, Node1)
cc02 <- dplyr::count(All2, Node2)
cc03 <- dplyr::count(All2, Node3)
cc04 <- dplyr::count(All2, Node4)
colnames(cc01)[1] <- "node"
colnames(cc02)[1] <- "node"
colnames(cc03)[1] <- "node"
colnames(cc04)[1] <- "node"

cc05 <- rbind(cc01, cc02, cc03,cc04)
cc05 <- cc05 %>% distinct(node, .keep_all = T)
cc05 <- cc05[-8,]
cc05$node

All2 <- All1 %>% dplyr::filter(Median_Forest_Score > 1)
g <- graph_from_data_frame(All2, directed = T)
E(g)$weight <- All2$Median_Forest_Score


net2 <- simplify(g, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum"))
net2 <- simplify(g, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum", "ignore"))



plot01 <- plot(g, layout = layout.reingold.tilford,
      edge.width= log(E(g)$weight),
      edge.arrow.width = 0.5,
      vertex.size = 5,
      edge.arrow.size = 0.3,
      vertex.size2 = 3,
      vertex.label.cex = 1,
      asp = 0.35,
      margin = 0.2)



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
