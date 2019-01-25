---
title: "Practicals GWAS. Association analysis and GWAS"
author: "Sodbo Sharapov"
date: "January 22, 2019"
output:
  html_document: 
    keep_md: true
---



# Outline of the practicals

The practical session consists of three steps:

 1. Exploring phenotype data
 2. Exploring genotype data and quality control
 4. Association analysis and GWAS

## Loading QCed ge03d2clean data

Load GenABEL package
 

```r
library(GenABEL)
```

```
## Loading required package: MASS
```

```
## Loading required package: GenABEL.data
```

Load phenotype and genotype data


```r
load('data/ge02d2clean.RData')
```

Let's check number of IDs and SNPs


```r
nids(ge03d2clean)
```

```
## [1] 893
```

```r
nsnps(ge03d2clean)
```

```
## [1] 7090
```

## Single SNP association analysis

Let's run association analysis for a single SNP.
For this we will store genotypes for this SNP in separate object.


```r
snpnames(ge03d2clean)[1:10]
```

```
##  [1] "rs1646456" "rs7950586" "rs4785242" "rs4435802" "rs2847446"
##  [6] "rs946364"  "rs299251"  "rs2456488" "rs1292700" "rs3712159"
```

```r
my_snp_num <- as.numeric(gtdata(ge03d2clean[,'rs1646456']))

my_snp_cha <- as.character(gtdata(ge03d2clean[,'rs1646456']))
```

Let's check the distribution of genotypes across samples


```r
table(my_snp_cha)
```

```
## my_snp_cha
## C/C C/G G/G 
## 409 381  91
```

```r
table(my_snp_num)
```

```
## my_snp_num
##   0   1   2 
## 409 381  91
```

Let's run logistic regresion with SNP as predictor and dm2 as outcome.


```r
log_reg_snp <- glm(phdata(ge03d2clean)$dm2~my_snp_num)

summary(log_reg_snp)
```

```
## 
## Call:
## glm(formula = phdata(ge03d2clean)$dm2 ~ my_snp_num)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -0.476  -0.472  -0.468   0.528   0.532  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  0.47578    0.02342   20.32   <2e-16
## my_snp_num  -0.00385    0.02547   -0.15     0.88
## 
## (Dispersion parameter for gaussian family taken to be 0.2498)
## 
##     Null deviance: 219.62  on 880  degrees of freedom
## Residual deviance: 219.62  on 879  degrees of freedom
##   (12 observations deleted due to missingness)
## AIC: 1282
## 
## Number of Fisher Scoring iterations: 2
```

As we can see, there is no association at all.
Let's check the distribution of cases and controls across genotypes.


```r
table(phdata(ge03d2clean)$dm2,my_snp_cha)
```

```
##    my_snp_cha
##     C/C C/G G/G
##   0 214 202  48
##   1 195 179  43
```

Indeed, the distribution of cases and controls across genotypes is random.

## Genome-wide association scan

Now let's run tiny GWAS :). We will use qtscore function of GenABEL package.
Basically, we will run 7000 regressions. We will run GWAS for dm2 with age and sex as covariates


```r
qt_simple <- mlreg(dm2~sex+age, data=ge03d2clean, trait="binomial")
```

Let's have a look at the top 10 signals


```r
descriptives.scan(qt_simple , sort="Pc1df")
```

```
## Summary for top 10 results, sorted by Pc1df
```

```
##           Chromosome Position Strand A1 A2   N    effB se_effB chi2.1df
## rs7903146          1  1047389      +  A  T 887 -0.7342  0.1507    23.72
## rs289981           1  1043860      -  T  G 882  0.4750  0.1103    18.56
## rs70099            2  8857747      +  C  A 883  0.8369  0.1984    17.80
## rs3436694          2  8921418      -  C  G 885  0.7212  0.1714    17.70
## rs7064741          1  1044233      -  C  G 885 -0.4617  0.1186    15.16
## rs2975760          3 10518480      +  A  T 888  0.3981  0.1059    14.14
## rs3074653          2  8915495      -  G  C 885  0.5156  0.1405    13.48
## rs5743183          1   648911      +  C  T 882 -0.6869  0.1976    12.09
## rs8541156          3 10397400      -  T  C 885 -0.4621  0.1332    12.04
## rs9386314          1   756075      -  C  A 889 -0.3941  0.1143    11.88
##                P1df     Pc1df effAB effBB chi2.2df P2df
## rs7903146 1.111e-06 4.389e-06    NA    NA       NA   NA
## rs289981  1.648e-05 4.879e-05    NA    NA       NA   NA
## rs70099   2.460e-05 6.978e-05    NA    NA       NA   NA
## rs3436694 2.593e-05 7.314e-05    NA    NA       NA   NA
## rs7064741 9.864e-05 2.416e-04    NA    NA       NA   NA
## rs2975760 1.696e-04 3.922e-04    NA    NA       NA   NA
## rs3074653 2.417e-04 5.386e-04    NA    NA       NA   NA
## rs5743183 5.071e-04 1.046e-03    NA    NA       NA   NA
## rs8541156 5.208e-04 1.071e-03    NA    NA       NA   NA
## rs9386314 5.664e-04 1.154e-03    NA    NA       NA   NA
```

As you can see, the top associated SNP has P-value of association of 3e-7,
which is not significant at the genome-wide significance level of 5e-8.

But what would happend if we will bmi as covariate into the model?


```r
qt_sex_age_bmi <- mlreg(dm2 ~ sex + age + bmi, ge03d2clean, trait="binomial")

descriptives.scan(qt_sex_age_bmi , sort="Pc1df")
```

```
## Summary for top 10 results, sorted by Pc1df
```

```
##           Chromosome Position Strand A1 A2   N    effB se_effB chi2.1df
## rs7903146          1  1047389      +  A  T 885 -0.9617  0.1689    32.41
## rs289981           1  1043860      -  T  G 880  0.5529  0.1191    21.55
## rs7064741          1  1044233      -  C  G 883 -0.5648  0.1281    19.45
## rs2975760          3 10518480      +  A  T 886  0.4824  0.1142    17.85
## rs6083205          3 10500286      +  T  G 885 -0.4570  0.1133    16.27
## rs3075920          1  1870260      -  G  T 883  0.5869  0.1552    14.29
## rs9468596          3 10506420      +  G  A 878  0.4338  0.1223    12.59
## rs6079246          2  7048058      +  A  T 882 -1.0377  0.2945    12.41
## rs8398809          1  1039030      +  A  C 880  0.3701  0.1099    11.33
## rs6227551          1  1044797      -  C  A 884 -0.9186  0.2735    11.28
##                P1df     Pc1df effAB effBB chi2.2df P2df
## rs7903146 1.250e-08 2.916e-08    NA    NA       NA   NA
## rs289981  3.454e-06 6.112e-06    NA    NA       NA   NA
## rs7064741 1.035e-05 1.737e-05    NA    NA       NA   NA
## rs2975760 2.392e-05 3.853e-05    NA    NA       NA   NA
## rs6083205 5.486e-05 8.488e-05    NA    NA       NA   NA
## rs3075920 1.563e-04 2.299e-04    NA    NA       NA   NA
## rs9468596 3.874e-04 5.456e-04    NA    NA       NA   NA
## rs6079246 4.260e-04 5.973e-04    NA    NA       NA   NA
## rs8398809 7.629e-04 1.040e-03    NA    NA       NA   NA
## rs6227551 7.838e-04 1.067e-03    NA    NA       NA   NA
```

Now we see one SNP with P-value of 1.056e-08, that is genome-wide significantly associated with T2D.

## Manhattan plot

Let's plot the Manhattan plot

```r
plot(qt_sex_age_bmi,df = 1)
```

![](GWAS_practicals.Part_3_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

## Regional assoc plot

Now let's zoom Manhattan plot around found signal. This plot is called regional association plot

```r
gwas_sum <- summary(qt_sex_age_bmi)
```

```
## Summary for top 10 results, sorted by P1df
```

```r
gwas_sum_sm <- gwas_sum[gwas_sum$Chromosome==1 & gwas_sum$Position>1047389-250000 & gwas_sum$Position<1047389+250000,]

plot(-log10(gwas_sum_sm$P1df)~gwas_sum_sm$Position)
```

![](GWAS_practicals.Part_3_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

## Table of genotypes and c/c status for the best SNP

Now let's have a look at the distribution of cases and controls across genotypes of the
top associated SNP.


```r
table(phdata(ge03d2clean)$dm2,as.character(gtdata(ge03d2clean[,'rs7903146'])))
```

```
##    
##     A/A A/T T/T
##   0 323 127  18
##   1 346  71   2
```

As you can see, the effect size is not big. Let's estimate OR for allele T.


```r
exp(0.7898)
```

```
## [1] 2.203
```

Let's compare results of mlreg with glm


```r
summary(glm(dm2~sex+age+bmi+as.numeric(gtdata(ge03d2clean[,'rs7903146'])),family='binomial', data=phdata(ge03d2clean)))
```

```
## 
## Call:
## glm(formula = dm2 ~ sex + age + bmi + as.numeric(gtdata(ge03d2clean[, 
##     "rs7903146"])), family = "binomial", data = phdata(ge03d2clean))
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -2.120  -0.996  -0.499   1.006   2.233  
## 
## Coefficients:
##                                                Estimate Std. Error z value
## (Intercept)                                    -3.96879    0.43055   -9.22
## sex                                             0.36877    0.14969    2.46
## age                                             0.01245    0.00582    2.14
## bmi                                             0.11187    0.01138    9.83
## as.numeric(gtdata(ge03d2clean[, "rs7903146"])) -0.96167    0.16893   -5.69
##                                                Pr(>|z|)
## (Intercept)                                     < 2e-16
## sex                                               0.014
## age                                               0.032
## bmi                                             < 2e-16
## as.numeric(gtdata(ge03d2clean[, "rs7903146"]))  1.3e-08
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1224.4  on 884  degrees of freedom
## Residual deviance: 1047.6  on 880  degrees of freedom
##   (8 observations deleted due to missingness)
## AIC: 1058
## 
## Number of Fisher Scoring iterations: 4
```

As you can see, the results between mlreg and glm are concordant.

## Saving summary stats


```r
sum_stats <- results(qt_sex_age_bmi)

snp_info <- summary(gtdata(ge03d2clean))

sum_stats2write <- sum_stats[,c('Chromosome','Position','A2','A1','N','effB','se_effB','P1df','Strand')]

sum_stats2write <- cbind(rownames(snp_info),rownames(snp_info),sum_stats2write)

sum_stats2write <- cbind(sum_stats2write[,c(1:7)],snp_info$Q.2,sum_stats2write[,c(8:11)])

colnames(sum_stats2write) <- c('rsid','snpid','chr','pos','a1','a0','n','freq1','beta1','se','p','strand')

write.table(sum_stats2write,file='results/sum_stats.tsv', sep='\t', quote = FALSE,row.names = FALSE)
```
## Additional literature: 

Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738

Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564–73. http://doi.org/10.1038/nprot.2010.116
## Saving data
