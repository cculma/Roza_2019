# annotate associated markers

library(tidyverse)
load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")

# # load gene annotation Medicago sativa Zhongmu No1
# load("~/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/i_5.2.8.RData")
# i_5.2.8 <- i_5.2.8 %>% dplyr::select(1,3)
head(i_5.2.8)

# # GRanges
# file <- ("~/Documents/Cesar/RNA/globus/lordec_reports/lordec_trim/bed_Shen/ORF_NMD/blast_corrected_shen.bed")
# txdb <- readBed(file, track.line = FALSE, remove.unusual = FALSE,
#                 zero.based = TRUE)

col_headings_1 <- c('gene_id',	'uniprot', 'gene_name',	'trans_length_flag',	'blastp_match_flag',	'nmd_flag',	'frame')
col_headings_2 <- c('gene_id',	'isoform')

QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T) 
data_4 <- QTL_02
gr5 <- GRanges(seqnames = data_4$Chrom,
               ranges = IRanges(data_4$Position, width = 1))

overlaps <- join_overlap_left(gr5, txdb)

df2 <- annoGR2DF(overlaps)
df2 <- unite(data = df2, col = "Marker1", 1:2, sep = "_", remove = F) %>% distinct(Marker1, .keep_all = TRUE) %>% dplyr::select(1:4,7) 
head(df2)
df3 <- df2 %>% separate(5, col_headings_1, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% separate(5, col_headings_2, sep = "\\.", remove = TRUE, convert = FALSE, extra = "warn") %>% dplyr::select(1,5,7) %>% inner_join(., i_5.2.8, by = "gene_id")

QTL_03 <- QTL_01 %>% group_by(Marker) %>% top_n(1, abs(Score)) %>% dplyr::select(Marker, Score)%>% distinct(Marker, .keep_all = TRUE)
QTL_04 <- QTL_01 %>% group_by(Marker) %>% summarise(Trait = paste(Trait, collapse = ";")) 
QTL_05 <- QTL_01 %>% group_by(Marker) %>% summarise(Model = paste(Model, collapse = ";")) 
QTL_06 <- QTL_01 %>% dplyr::select(Marker, Chrom, Position, Ref, Alt) %>% distinct(Marker, .keep_all = TRUE) %>% unite(col = "SNP", 5:4, sep = "/", remove = T) %>% unite(col = "Marker1", 2:3, sep = "_", remove = T) 

QTL_08 <- inner_join(QTL_06, QTL_03, by = "Marker") %>% inner_join(., QTL_04, by = "Marker") %>% inner_join(., QTL_05, by = "Marker") %>% left_join(., df3, by = "Marker1") %>% dplyr::select(-1)

nrow(QTL_08 %>% distinct(gene_id, .keep_all = TRUE))
sum(!is.na(QTL_08$gene_id))

QTL_09 <- QTL_08 %>% group_by(gene_id) %>% summarise(Marker1 = paste(Marker1, collapse = ";")) 
QTL_10 <- inner_join(QTL_09, QTL_10, by = "gene_id")
QTL_10 <- QTL_08 %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(6:8)
QTL_10 <- na.omit(QTL_10)
head(QTL_10)
write.table(QTL_10, "~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/QTL_10.tsv", row.names = F, quote = F, sep = "\t")
