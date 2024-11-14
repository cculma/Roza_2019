

setwd("~/Documents/git/big_files/")

PCA <- read.csv("PCA_Roza2019.csv")
head(PCA)

PCA1 <- PCA %>% dplyr::select( Susceptible_parent1, Roza2019_VCF_ID)
head(PCA1)
# 
# PCA$pop2 <- recode_factor(PCA$Susceptible_parent, 
#                           "1613-1"  = "1", "1613-10" = "2", "1613-11" = "3", 
#                           "1613-12" = "4", "1613-13" = "5", "1613-14" = "6",
#                           "1613-15" = "7", "1613-16" = "8", "1613-17" = "9",
#                           "1613-18" = "10", "1613-19" = "11", "1613-2" = "12",
#                           "1613-20" = "13", "1613-21" = "14", "1613-22" = "15",
#                           "1613-23" = "16", "1613-25" = "17", "1613-26" = "18", 
#                           "1613-28" = "19", "1613-29" = "20", "1613-3" = "21", 
#                           "1613-30" "1613-31" "1613-32" , 
#                           "1613-4"  "1613-5" "1613-6" ,
#                           "1613-7"  "1613-8"  "1613-9" )
                          
cc06 <- dplyr::count(PCA1, Susceptible_parent1)
cc06$Susceptible_parent1

G1 <- read.table("Roza2019_06_GS.txt", sep = "\t", row.names = 1)
G1[1:5,1:5]
str(G1)
dim(G1)
G1 <- G1 %>% mutate_all(as.character)
G1 <- as.matrix(G1)

G1[G1 == "0"] <- "001001001001"
G1[G1 == "1"] <- "001001001002"
G1[G1 == "2"] <- "001001002002"
G1[G1 == "3"] <- "001002002002"
G1[G1 == "4"] <- "002002002002"

df4 <- as.data.frame(G1)
df4 <- df4 %>% rownames_to_column("Roza2019_VCF_ID")
df4[1:5,1:5]

col9.2 <- inner_join(PCA1, df4, by = "Roza2019_VCF_ID")
colnames(col9.2)[1] <- "pop"
colnames(col9.2)[2] <- "ind"
col9.2[1:10, 1:5]

setwd("~/Documents/git/big_files/")
write.table(col9.2, "Roza2019_06_genodive.txt", sep = "\t", col.names = T, row.names = F, quote = F) # 2247 2999
