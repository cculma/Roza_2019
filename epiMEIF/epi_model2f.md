> plot1 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)   
1    416 574972                         
2    410 548209  6     26763 0.002752 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot2 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    417 603479                           
2    412 563655  5     39824 2.207e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot3 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr6.1_27950240"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    416 598191                      
2    410 593839  6    4351.9   0.8083
> plot4 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434","chr6.1_27950240"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    413 572661                           
2    392 485607 21     87053 3.177e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot5 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_11234756"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    416 580048                        
2    411 559875  5     20172  0.01121 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot6 <- plotSNPInteraction(data4, c("chr6.1_27950240", "chr7.1_11234756"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr6.1_27950240 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    415 577333                        
2    407 556622  8     20711  0.05641 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot7 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240","chr7.1_11234756"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr6.1_27950240 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)   
1    412 562568                         
2    388 497424 24     65143 0.001115 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot8 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr5.1_1429434"))
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr1.1_18947323 + chr5.1_1429434
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    417 583226                        
2    413 567487  4     15739  0.02191 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
