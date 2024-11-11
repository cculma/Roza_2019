Longitudinal X23

> plot1 <- plotSNPInteraction(data4, c("Stress", "chr6.1_71935085"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ Stress + chr6.1_71935085
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df     RSS Df Sum of Sq  Pr(>Chi)    
1    842 1913699                           
2    838 1849438  4     64261 7.401e-06 ***

plot2 <- plotSNPInteraction(data4, c("Stress", "chr4.1_79960120", "chr6.1_71935085"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ Stress + chr4.1_79960120 + chr6.1_71935085
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df     RSS Df Sum of Sq  Pr(>Chi)    
1    839 1877804                           
2    818 1775065 21    102739 0.0008437 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1