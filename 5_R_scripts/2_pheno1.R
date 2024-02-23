head(ST1)
head(ST2)
head(ST3)
head(ST4)
hist(ST4$predicted.value)

ST2$month <- paste0("ST2_", ST2$month)
ST3$year <- paste0("ST3_", ST3$year)


colnames(ST4)[2] <- "ST4_Yi"
ST4 <- ST4 %>% select(1:2)

ST2 <- ST2 %>% select(1:3) %>% spread(key = month, value = predicted.value) 
ST3 <- ST3 %>% select(1:3) %>% spread(key = year, value = predicted.value)
ST2 <- ST2[,c(1,5,4,3,2,6)]

# BLUP5 <- inner_join(ST2, ST3, by = "gen") %>% inner_join(., ST4, by = "gen") %>% column_to_rownames("gen")
BLUP4 <- inner_join(ST2, ST3, by = "gen") %>% inner_join(., ST4, by = "gen")

BLUP5 <- BLUP4 %>% gather(key = "env", value = "BLUE", 2:10)
BLUP5$env <- factor(BLUP5$env, levels = c("ST2_may","ST2_jun","ST2_jul","ST2_aug","ST2_sep","ST3_20","ST3_21","ST3_22","ST4_Yi"))
levels(BLUP5$env)

ggplot(BLUP5, aes(x = env, y = BLUE)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial", base_size = 12) + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST0", y = "BLUE Yield (lb * plot * 100)", x = "") + ylim(0, 100)

BLUP5 <- BLUP4 %>% column_to_rownames("gen")
BLUP5 <- cor(BLUP5, use = "complete")

ggcorrplot(BLUP5[,ncol(BLUP5):1], hc.order = F, type = "full", lab = T, lab_col = "grey3", lab_size = 3, show.diag = T) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title.x=element_blank(), axis.title.y = element_blank()) + labs(title = "Single-Stage vs Stage-Wise")

Y2.1 <- Y2
Y2.1$env <- paste0("ST0_", Y2.1$env)
ST0 <- Y2.1 %>% select(1:3) %>% spread(key = env, value = predicted.values)

head(ST1)
ST1$env <- paste0("ST1_", ST1$env)
ST1 <- ST1 %>% select(1:3) %>% spread(key = env, value = predicted.value)

BLUP5 <- inner_join(ST0, ST1, by = "gen") %>% inner_join(., ST2, by = "gen") %>% inner_join(., ST3, by = "gen")  %>% inner_join(., ST4, by = "gen") %>% inner_join(., PCA, by = "gen") 

a3 <- read.csv("~/Documents/git/Roza_2019/pheno_data/pheno_nph.csv")
a4 <- read.csv("~/Documents/git/Roza_2019/pheno_data/VCF_IDs.csv")

a3 <- a3[,c(1,25:27)]
a4 <- a4[,c(3,8)]
colnames(a4) <- c("gen","all")
a2 <- inner_join(a3,a4, by = "all")

BLUP3 <- inner_join(BLUP5, a2, by = "gen")
colnames(BLUP3)
BLUP3 <- BLUP3[,c(45,2:44,46:48)]
# write.csv(BLUP3, "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/pheno_yi.csv", quote = F, row.names = F)
# setwd("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/")