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
Basically, we will run 7000 regressions.


```r
qt_simple <- qtscore(dm2, ge03d2clean, trait="binomial")

descriptives.scan(qt_simple , sort="Pc1df")
```

```
## Summary for top 10 results, sorted by Pc1df
```

```
##           Chromosome Position Strand A1 A2   N   effB se_effB chi2.1df
## rs7903146          1  1047389      +  A  T 887 0.4732 0.09238    26.24
## rs289981           1  1043860      -  T  G 882 1.6084 0.35859    20.12
## rs3436694          2  8921418      -  C  G 885 2.0983 0.47455    19.55
## rs70099            2  8857747      +  C  A 883 2.3164 0.52491    19.47
## rs7064741          1  1044233      -  C  G 885 0.6176 0.14989    16.98
## rs2975760          3 10518480      +  A  T 888 1.4799 0.38810    14.54
## rs3074653          2  8915495      -  G  C 885 1.6536 0.44798    13.63
## rs5743183          1   648911      +  C  T 882 0.5026 0.14067    12.76
## rs9386314          1   756075      -  C  A 889 0.6716 0.19058    12.42
## rs1801282          2  8931192      +  C  T 882 1.6314 0.46993    12.05
##                P1df  effAB  effBB chi2.2df      P2df     Pc1df
## rs7903146 3.008e-07 0.5219 0.1037    26.80 1.512e-06 1.550e-06
## rs289981  7.280e-06 1.3125 2.3067    21.00 2.751e-05 2.591e-05
## rs3436694 9.791e-06 1.9610    Inf    20.49 3.552e-05 3.368e-05
## rs70099   1.020e-05 2.1524    Inf    20.04 4.449e-05 3.491e-05
## rs7064741 3.784e-05 0.5769 0.4686    17.83 1.343e-04 1.114e-04
## rs2975760 1.372e-04 1.1014 3.3707    22.84 1.097e-05 3.484e-04
## rs3074653 2.231e-04 1.5341 4.0687    14.31 7.798e-04 5.361e-04
## rs5743183 3.533e-04 0.4646 0.6684    13.97 9.240e-04 8.059e-04
## rs9386314 4.250e-04 0.5848 0.6154    15.00 5.542e-04 9.493e-04
## rs1801282 5.174e-04 1.5987 3.0346    12.10 2.359e-03 1.130e-03
```


```r
qt_sex_age_bmi <- qtscore(dm2 ~ sex + age + bmi, ge03d2clean, trait="binomial")
```

```
## Warning in qtscore(dm2 ~ sex + age + bmi, ge03d2clean, trait = "binomial"):
## 2 observations deleted due to missingness
```

```r
descriptives.scan(qt_sex_age_bmi , sort="Pc1df")
```

```
## Summary for top 10 results, sorted by Pc1df
```

```
##           Chromosome Position Strand A1 A2   N   effB se_effB chi2.1df
## rs7903146          1  1047389      +  A  T 885 0.7898  0.1380    32.74
## rs289981           1  1043860      -  T  G 880 1.1503  0.2515    20.91
## rs7064741          1  1044233      -  C  G 883 0.8648  0.2004    18.63
## rs2975760          3 10518480      +  A  T 886 1.1306  0.2749    16.92
## rs6083205          3 10500286      +  T  G 885 0.8900  0.2325    14.65
## rs3075920          1  1870260      -  G  T 883 1.1634  0.3165    13.51
## rs9195215          1  1464966      -  G  C 880 0.8854  0.2529    12.26
## rs9468596          3 10506420      +  G  A 878 1.1185  0.3208    12.16
## rs8889070          3 10220702      -  G  A 879 0.8559  0.2463    12.08
## rs6308770          1  1464229      -  T  A 883 0.8897  0.2643    11.33
##                P1df  effAB  effBB chi2.2df      P2df     Pc1df
## rs7903146 1.056e-08 0.8020 0.5825    32.98 6.879e-08 3.297e-08
## rs289981  4.804e-06 1.0807 1.2777    21.84 1.811e-05 1.005e-05
## rs7064741 1.585e-05 0.8527 0.7776    19.01 7.430e-05 3.070e-05
## rs2975760 3.906e-05 1.0489 1.4017    23.32 8.649e-06 7.139e-05
## rs6083205 1.293e-04 0.8516 0.7745    15.16 5.100e-04 2.189e-04
## rs3075920 2.374e-04 1.1526 1.4399    13.73 1.045e-03 3.865e-04
## rs9195215 4.630e-04 0.8912 0.7712    12.31 2.126e-03 7.223e-04
## rs9468596 4.881e-04 1.0424 1.4192    18.48 9.730e-05 7.589e-04
## rs8889070 5.100e-04 0.8553 0.7375    12.08 2.378e-03 7.908e-04
## rs6308770 7.634e-04 0.8963 0.7765    11.40 3.354e-03 1.154e-03
```

## Manhattan plot

## Regional assoc plot

## Table of genotypes and c/c status for the best SNP

## Saving summary stats

## Additional literature: 

Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738

Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564–73. http://doi.org/10.1038/nprot.2010.116
## Saving data
