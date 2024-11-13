> plot1 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240")) # 6
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)   
1    416 574972                         
2    410 548209  6     26763 0.002752 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot2 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434")) # 3
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    417 603479                           
2    412 563655  5     39824 2.207e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot3 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr6.1_27950240")) # 5
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    416 598191                      
2    410 593839  6    4351.9   0.8083
> plot4 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434","chr6.1_27950240")) # 7
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    413 572661                           
2    392 485607 21     87053 3.177e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot5 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr7.1_11234756")) # 9
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    416 580048                        
2    411 559875  5     20172  0.01121 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot6 <- plotSNPInteraction(data4, c("chr6.1_27950240", "chr7.1_11234756")) # 27
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr6.1_27950240 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    415 577333                        
2    407 556622  8     20711  0.05641 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot7 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr6.1_27950240","chr7.1_11234756")) # 60
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr6.1_27950240 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)   
1    412 562568                         
2    388 497424 24     65143 0.001115 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot8 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr6.1_27950240")) # 18
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr1.1_18947323 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    416 581400                      
2    411 576398  5    5001.7   0.6134
> plot9 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr6.1_27950240")) # 17
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr1.1_18947323 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    416 581400                      
2    411 576398  5    5001.7   0.6134
> plot10 <- plotSNPInteraction(data4, c("chr1.1_18947323", "chr5.1_1429434","chr6.1_27950240")) # 19
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr1.1_18947323 + chr5.1_1429434 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
1    413 565906                           
2    394 501644 19     64261 0.0001117 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot12 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr7.1_11234756")) # 10
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    416 588161                      
2    409 582708  7      5453   0.7994
> plot11 <- plotSNPInteraction(data4, c("chr2.1_41048095", "chr5.1_1429434","chr7.1_11234756")) # 11
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr2.1_41048095 + chr5.1_1429434 + chr7.1_11234756
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)  
1    413 571602                        
2    392 528160 21     43442  0.05532 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot13 <- plotSNPInteraction(data4, c("chr5.1_1429434", "chr5.1_54813818")) # 22
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr5.1_1429434 + chr5.1_54813818
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)   
1    417 601373                         
2    411 570406  6     30967 0.001062 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> plot14 <- plotSNPInteraction(data4, c("chr4.1_20948923", "chr6.1_27950240")) # 24
Analysis of Variance Table

Model 1: PHENOTYPE ~ chr4.1_20948923 + chr6.1_27950240
Model 2: PHENOTYPE ~ Genotype_Combination
  Res.Df    RSS Df Sum of Sq Pr(>Chi)
1    415 581261                      
2    407 572033  8    9227.9   0.5841