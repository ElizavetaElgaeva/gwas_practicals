Outline of the practicals
=========================

The practical session consists of XXX steps:

1.  Exploring phenotype data
2.  Exploring genotype data
3.  Quality control of genotype data
4.  Association analysis and GWAS
5.  Vizualization of the results

Exploring phenotype data
========================

Load GenABEL package

``` r
library(GenABEL)
```

    ## Loading required package: MASS

    ## Loading required package: GenABEL.data

Load phenotype and genotype data

``` r
data("ge03d2")
```

Let's get the Structure of ge03d2 object

``` r
str(ge03d2)
```

    ## Formal class 'gwaa.data' [package "GenABEL"] with 2 slots
    ##   ..@ phdata:'data.frame':   950 obs. of  8 variables:
    ##   .. ..$ id    : chr [1:950] "id4" "id10" "id25" "id33" ...
    ##   .. ..$ sex   : int [1:950] 0 1 0 0 0 0 1 1 0 1 ...
    ##   .. ..$ age   : num [1:950] 51.6 53.7 66 44.7 49.9 ...
    ##   .. ..$ dm2   : int [1:950] 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..$ height: num [1:950] 152 170 165 174 157 ...
    ##   .. ..$ weight: num [1:950] 83 65.6 69.6 69.7 84.3 ...
    ##   .. ..$ diet  : int [1:950] 0 0 0 0 0 0 0 0 0 0 ...
    ##   .. ..$ bmi   : num [1:950] 36 22.6 25.6 22.9 34.1 ...
    ##   ..@ gtdata:Formal class 'snp.data' [package "GenABEL"] with 11 slots
    ##   .. .. ..@ nbytes    : num 238
    ##   .. .. ..@ nids      : int 950
    ##   .. .. ..@ nsnps     : int 7589
    ##   .. .. ..@ idnames   : chr [1:950] "id4" "id10" "id25" "id33" ...
    ##   .. .. ..@ snpnames  : chr [1:7589] "rs1646456" "rs7950586" "rs4785242" "rs4435802" ...
    ##   .. .. ..@ chromosome: Factor w/ 4 levels "1","2","3","X": 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. .. .. ..- attr(*, "names")= chr [1:7589] "rs1646456" "rs7950586" "rs4785242" "rs4435802" ...
    ##   .. .. ..@ map       : Named num [1:7589] 653 849 1766 5291 5555 ...
    ##   .. .. .. ..- attr(*, "names")= chr [1:7589] "rs1646456" "rs7950586" "rs4785242" "rs4435802" ...
    ##   .. .. ..@ coding    :Formal class 'snp.coding' [package "GenABEL"] with 1 slot
    ##   .. .. .. .. ..@ .Data: raw [1:7589] 11 07 09 0f ...
    ##   .. .. ..@ strand    :Formal class 'snp.strand' [package "GenABEL"] with 1 slot
    ##   .. .. .. .. ..@ .Data: raw [1:7589] 01 02 02 01 ...
    ##   .. .. ..@ male      : Named int [1:950] 0 1 0 0 0 0 1 1 0 1 ...
    ##   .. .. .. ..- attr(*, "names")= chr [1:950] "id4" "id10" "id25" "id33" ...
    ##   .. .. ..@ gtps      :Formal class 'snp.mx' [package "GenABEL"] with 1 slot
    ##   .. .. .. .. ..@ .Data: raw [1:238, 1:7589] 99 57 9a 6a ...

Let's save phenotype data into separate object named *pheno*. And look at the first rows of the *pheno* object

``` r
pheno <- phdata(ge03d2)

head(pheno)
```

    ##        id sex   age dm2 height weight diet   bmi
    ## id4   id4   0 51.64   1  151.8  83.03    0 36.04
    ## id10 id10   1 53.74   1  170.3  65.60    0 22.61
    ## id25 id25   0 66.01   1  164.8  69.55    0 25.61
    ## id33 id33   0 44.67   1  174.5  69.73    0 22.91
    ## id35 id35   0 49.88   1  157.2  84.28    0 34.11
    ## id58 id58   0 37.00   1  165.1  89.54    0 32.85

Description of the columns:

1.  ID: id code of sample
2.  sex: sex, '0' is female, '1' is male
3.  age: age, in years
4.  dm2: diabetes mellitus type 2, '0' is unaffected, '1' is affected
5.  height: height, in cm
6.  wieght: wieght, in kg
7.  diet: following the diet, '0' is false, '1' is true
8.  bmi: BMI, body mass index, in *c**m*/*k**g*<sup>2</sup>

Now let's have a look at the summary of the phenotype data

``` r
summary(pheno)
```

    ##       id                 sex             age            dm2       
    ##  Length:950         Min.   :0.000   Min.   :18.6   Min.   :0.000  
    ##  Class :character   1st Qu.:0.000   1st Qu.:40.0   1st Qu.:0.000  
    ##  Mode  :character   Median :0.000   Median :49.7   Median :0.000  
    ##                     Mean   :0.489   Mean   :50.0   Mean   :0.487  
    ##                     3rd Qu.:1.000   3rd Qu.:59.3   3rd Qu.:1.000  
    ##                     Max.   :1.000   Max.   :87.6   Max.   :1.000  
    ##                                                                   
    ##      height        weight           diet             bmi      
    ##  Min.   :140   Min.   : 34.1   Min.   :0.0000   Min.   :13.5  
    ##  1st Qu.:161   1st Qu.: 66.7   1st Qu.:0.0000   1st Qu.:24.3  
    ##  Median :168   Median : 78.8   Median :0.0000   Median :27.9  
    ##  Mean   :168   Mean   : 84.0   Mean   :0.0442   Mean   :29.7  
    ##  3rd Qu.:175   3rd Qu.: 96.2   3rd Qu.:0.0000   3rd Qu.:34.1  
    ##  Max.   :200   Max.   :161.2   Max.   :1.0000   Max.   :64.2  
    ##  NA's   :2     NA's   :2                        NA's   :2

For the categorial traits **table** function can be applied to get summary of the trait.

``` r
table(pheno$sex)
```

    ## 
    ##   0   1 
    ## 485 465

``` r
table(pheno$dm2)
```

    ## 
    ##   0   1 
    ## 487 463

``` r
table(pheno$diet)
```

    ## 
    ##   0   1 
    ## 908  42

You can plot a historgram for quantitaive traits

``` r
hist(pheno$age, main='Histogram of age', xlab = 'Age')
```

![](GWAS_practicals_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
hist(pheno$height, main='Histogram of height', xlab = 'Height')
```

![](GWAS_practicals_files/figure-markdown_github/unnamed-chunk-7-2.png)

``` r
hist(pheno$weight, main='Histogram of weight', xlab = 'Weight')
```

![](GWAS_practicals_files/figure-markdown_github/unnamed-chunk-7-3.png)

``` r
hist(pheno$bmi, main='Histogram of bmi', xlab = 'BMI')
```

![](GWAS_practicals_files/figure-markdown_github/unnamed-chunk-7-4.png)

It is important to check whether your trait of interest is correlated with covariates. To eheck this you can estimate the Peasron correlation coefficient and/or perform linear regression analysis. First, let's examine the correlation structure for phenotype data

``` r
cor_matrix <- cor(pheno[,c('dm2','sex','age','height','weight','bmi','diet')], use='complete.obs')

cor_matrix
```

    ##             dm2      sex       age   height   weight      bmi      diet
    ## dm2     1.00000  0.11345  0.100623  0.05504 0.372321  0.36080 -0.046280
    ## sex     0.11345  1.00000  0.029226  0.61193 0.326259  0.08408 -0.025769
    ## age     0.10062  0.02923  1.000000 -0.24960 0.007096  0.10903 -0.150986
    ## height  0.05504  0.61193 -0.249596  1.00000 0.367637 -0.04506  0.038821
    ## weight  0.37232  0.32626  0.007096  0.36764 1.000000  0.90628  0.005094
    ## bmi     0.36080  0.08408  0.109029 -0.04506 0.906279  1.00000 -0.004970
    ## diet   -0.04628 -0.02577 -0.150986  0.03882 0.005094 -0.00497  1.000000

Now let's vizualize correlation matrix using *corplot* package

``` r
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
corrplot(cor_matrix)
```

![](GWAS_practicals_files/figure-markdown_github/unnamed-chunk-9-1.png)

Heatmap clearly shows that dm2 is correlated with weight and bmi. Let's run logistic regression for dm2 vs weight

``` r
dm2_weight <- glm(dm2~weight,data = pheno,family = 'binomial')
summary(dm2_weight)
```

    ## 
    ## Call:
    ## glm(formula = dm2 ~ weight, family = "binomial", data = pheno)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.896  -1.033  -0.698   1.085   1.988  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  -3.0583     0.2878   -10.6   <2e-16
    ## weight        0.0361     0.0034    10.6   <2e-16
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1313.7  on 947  degrees of freedom
    ## Residual deviance: 1172.4  on 946  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 1176
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
dm2_bmi <- glm(dm2~bmi,data = pheno,family = 'binomial')
summary(dm2_bmi)
```

    ## 
    ## Call:
    ## glm(formula = dm2 ~ bmi, family = "binomial", data = pheno)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -2.215  -1.024  -0.692   1.086   1.971  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  -3.2657     0.3171   -10.3   <2e-16
    ## bmi           0.1096     0.0107    10.2   <2e-16
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1313.7  on 947  degrees of freedom
    ## Residual deviance: 1177.2  on 946  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: 1181
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
dm2_sex <- glm(dm2~sex,data = pheno,family = 'binomial')
summary(dm2_sex)
```

    ## 
    ## Call:
    ## glm(formula = dm2 ~ sex, family = "binomial", data = pheno)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ##  -1.25   -1.07   -1.07    1.10    1.29  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  -0.2697     0.0916   -2.94  0.00325
    ## sex           0.4465     0.1306    3.42  0.00063
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1316.4  on 949  degrees of freedom
    ## Residual deviance: 1304.6  on 948  degrees of freedom
    ## AIC: 1309
    ## 
    ## Number of Fisher Scoring iterations: 3

betaOR for sex is 0.4465, that corresponds to `OR=exp(0.4465)=1.563`, which concordant with that, calculated from 2x2 table:

``` r
table(pheno$dm2,pheno$sex)
```

    ##    
    ##       0   1
    ##   0 275 212
    ##   1 210 253

``` r
275*253/210/212
```

    ## [1] 1.563
