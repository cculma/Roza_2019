> plot1 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_24014763"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr7.1_24014763
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    417 326594                      
2    413 324269  4    2324.4   0.5645

> plot2 <- plotSNPInteraction(data4, c("chr4.1_36675681", "chr4.1_60953306"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr4.1_36675681 + chr4.1_60953306
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    418 329187                      
2    414 327542  4    1645.1   0.7212

> plot3 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_24014763"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr7.1_24014763
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    417 326594                      
2    413 324269  4    2324.4   0.5645