#top_GO
rm(list = ls())
# BiocManager::install("topGO")
library(topGO)

#######################
# This part of code allows to use GO_id from allele_aware to ZhongmuNo1

setwd("~/Documents/Cesar/RNA/globus/lordec_reports/tama_merge_2/fusion_all2")
a1 <- read.table('~/Documents/Cesar/RNA/globus/lordec_reports/tama_merge_2/fusion_all2/background.csv', header = F, sep = '\t')
colnames(a1) <- c('uniprot', 'GO_id')

load("/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/i_5.2.3.RData")
i_5.2.4 <- i_5.2.3 %>% separate(4, col_headings_1, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% separate(4, col_headings_2, sep = "\\.", remove = TRUE, convert = FALSE, extra = "warn") %>% dplyr::select(4,6) %>% distinct(gene_id, .keep_all = TRUE)
head(i_5.2.4)
a2 <- inner_join(i_5.2.4, a1, by = "uniprot")  %>% dplyr::select(-1)
head(a2)
dim(a2)
write.table(a2, '/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/background_3.tsv', sep = '\t', col.names = F, row.names = F, quote = F)
a2 <- read.table('/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/background_3.tsv', header = F, sep = '\t')

##########################
# this part of the code allows to generate QTL_10 table

# load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")

QTL_08 <- inner_join(QTL_06, QTL_03, by = "Marker") %>% inner_join(., QTL_04, by = "Marker") %>% inner_join(., QTL_05, by = "Marker") %>% left_join(., df3, by = "Marker1") %>% dplyr::select(-1)
QTL_09 <- QTL_08 %>% group_by(gene_id) %>% summarise(Marker1 = paste(Marker1, collapse = ";")) 
QTL_10 <- inner_join(QTL_09, QTL_10, by = "gene_id")
QTL_10 <- QTL_08 %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(6:8)
QTL_10 <- na.omit(QTL_10)
head(QTL_10)

##########################
# this part of code is the GO enrichment

GOesByID.1 <- readMappings(file = '/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/background_4.tab')
GOesByID.2 <- readMappings(file = '/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/background_3.tsv')

head(GOesByID.1)
head(GOesByID.2)
bg_genes <- names(GOesByID.1)
bg_genes <- names(GOesByID.2)
bg_genes[1:5]

compared_genes <- factor(as.integer(bg_genes %in% QTL_10$uniprot))
names(compared_genes) <- bg_genes
head(compared_genes)

myGOdata_BP <- new("topGOdata", description="My project", ontology="BP", allGenes=compared_genes, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID.2)

myGOdata_CC <- new("topGOdata", description="My project", ontology="CC", allGenes=compared_genes, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID.2)

myGOdata_MF <- new("topGOdata", description="My project", ontology="MF", allGenes=compared_genes, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = GOesByID.2)

resultFisher_BP <- runTest(myGOdata_BP, algorithm = "weight01", statistic = "fisher")
resultFisher_CC <- runTest(myGOdata_CC, algorithm = "weight01", statistic = "fisher")
resultFisher_MF <- runTest(myGOdata_MF, algorithm = "weight01", statistic = "fisher")

tab_BP <- GenTable(myGOdata_BP, raw.p.value = resultFisher_BP, numChar = 120, topNodes = 20)
tab_CC <- GenTable(myGOdata_CC, raw.p.value = resultFisher_CC, numChar = 120, topNodes = 20)
tab_MF <- GenTable(myGOdata_MF, raw.p.value = resultFisher_MF, numChar = 120, topNodes = 40)

tab_BP <- subset(tab_BP, Significant > 0)
tab_CC <- subset(tab_CC, Significant > 0)
tab_MF <- subset(tab_MF, Significant > 0)

tab_BP$GO_process <- "BP"
tab_CC$GO_process <- "CC"
tab_MF$GO_process <- "MF"

# this part requires to change c(tab_BP, tab_CC, tab_MF) and c(myGOdata_BP, myGOdata_CC, myGOdata_MF)
myterms =tab_MF$GO.ID
mygenes = genesInTerm(myGOdata_MF, myterms)
tab_MF$genes <- sapply(tab_MF$GO.ID, function(x)
{
  genes <- genesInTerm(myGOdata_MF, x) 
  genes[[1]][genes[[1]] %in% sigGenes(myGOdata_MF)] # significant genes
})
# ~~~~~~~~
GO_associated <- rbind(tab_BP, tab_CC, tab_MF)

a1 <- unnest(GO_associated, genes) %>% dplyr::select(1,2,7,8)
head(a1)
colnames(a1)[4] <- "uniprot"
QTL_11 <- a1 %>% group_by(uniprot) %>% summarise(GO.ID = paste(GO.ID, collapse = ";"))
head(QTL_10)
QTL_07 <- left_join(QTL_10, QTL_11, by = "uniprot")

write.table(QTL_07, "~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/QTL_07.tsv", row.names = F, quote = F, sep = "\t")
