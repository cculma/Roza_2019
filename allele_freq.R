library(sommer)


# MPP_Ms2_GWASPoly.txt
G <- read.csv('~/Documents/Cesar/blup_data/Roza2019/GWASPOLY/MPP_Ms2_GWASPoly.txt', header = TRUE, row.names = 1, check.names = F)
G[1:5,1:5]
dim(G)
G1 <- G[1:1000,]
G1[1:5,1:5]
dim(G1)

G1 <- G %>% unite(Chrom1, 1:2, remove = T)
G3 <- G[,c(1,2)]
G1[1:5,1:5]
G1 <- as.matrix(G1 %>% remove_rownames() %>% column_to_rownames(var = "Chrom1"))

G2 <- t(G1)
G2[1:5,1:5]
G2 <- as.data.frame(G2)
str(G2)
numo <- atcg1234(data=G2, ploidy=4); 
numo$M[1:5,1:5] 
numo$ref.allele[,1:5]

pheno <- read.csv("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/pheno_nph.csv", row.names = 1)
head(pheno)
common <- intersect(rownames(pheno),rownames(numo$M))

marks <- numo$M[common,]
marks[1:5,1:5]
dim(marks)
str(marks)
marks.1 <- as.data.frame(t(marks))
marks.1[1:5,1:5]
dim(marks.1)
marks.2 <- tibble::rownames_to_column(marks.1, "Marker1")
marks.2[1:5,1:5]
marks.3 <- inner_join(QTL_06, marks.2, by = "Marker1") %>% dplyr::select(-c(1,3)) %>% column_to_rownames(var = "Marker1")
marks.3[1:5,1:5]
marks.3 <- as.data.frame(t(marks.3)) %>% mutate_if(is.integer,as.factor)
str(marks.3)

marks.3$Chr1_10270033
a0 <- marks.3[marks.3$Chr1_10270033 %in% c(0), ]
a1 <- marks.3[marks.3$Chr1_10270033 %in% c(1), ]
a2 <- marks.3[marks.3$Chr1_10270033 %in% c(2), ]
a3 <- marks.3[marks.3$Chr1_10270033 %in% c(3), ]
a4 <- marks.3[marks.3$Chr1_10270033 %in% c(4), ]
rownames(a0)


c0 <- intersect(rownames(pheno), rownames(a0))
p0 <- pheno[c0, ]
p0 <- pheno[,1, drop = F]
p0 <- merge(p0, marks.3, by = "row.names")
p1 <- merge(as.data.frame(p0), as.data.frame(marks.3), by = 'row.names', all = TRUE)
rownames(p0)


pheno$may_20_1stage
rownames(pheno)
boxplot()
G6 <- as.data.frame(summary(marks.3))
G6 <- G6[,-1]
G6 <- separate(data = G6, col = 2, into = c("MARKER", "count"), remove = T, sep = ":")
colnames(G6)
G6$count <- as.numeric(G6$count)
G6 <- na.omit(G6)
G7 <- spread(data = G6, key = MARKER, value = count, )
colnames(G7) <- c("Marker1", "0", "1", "2", "3", "4")

G8 <- QTL_01 %>% dplyr::select(Marker, Chrom, Position, Ref, Alt, Effect) %>% distinct(Marker, .keep_all = TRUE) %>% unite(col = "SNP", 5:4, sep = "/", remove = T) %>% unite(col = "Marker1", 2:3, sep = "_", remove = T)

G9 <- inner_join(G8, G7, by = "Marker1")
