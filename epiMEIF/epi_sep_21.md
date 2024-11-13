> plot1_sep21 <- plotSNPInteraction(data4, c("chr5.1_42131042", "chr6.1_49779667"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_42131042 + chr6.1_49779667
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    415 319758                           
2    407 299955  8     19803 0.0007445 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> plot02_sep21 <- plotSNPInteraction(data4, c("chr5.1_42131042", "chr6.1_67638734"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_42131042 + chr6.1_67638734
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    415 328892                           
2    402 300843 13     28050 0.0003481 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 