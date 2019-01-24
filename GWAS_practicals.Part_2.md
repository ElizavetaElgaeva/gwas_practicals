---
title: "Practicals GWAS. Exploring genotype data and QC"
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

## Loading genotype data

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
data("ge03d2")
```

Let's get the Structure of ge03d2 object


```r
str(ge03d2)
```

```
## Formal class 'gwaa.data' [package "GenABEL"] with 2 slots
##   ..@ phdata:'data.frame':	950 obs. of  8 variables:
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
```

Genotypes are stored as a huge N x M table, where N is number of IDs and M is a number of SNPs.


```r
as.numeric(gtdata(ge03d2[1:5,1:5]))
```

```
##      rs1646456 rs7950586 rs4785242 rs4435802 rs2847446
## id4          1         1         1         1         0
## id10         0         0         1         0         1
## id25         1         0         2         1         1
## id33         0         0         1         0         1
## id35         0         0         2         0         2
```

```r
as.character(gtdata(ge03d2[1:5,1:5]))
```

```
##      rs1646456 rs7950586 rs4785242 rs4435802 rs2847446
## id4  "C/G"     "T/A"     "T/C"     "C/A"     "T/T"    
## id10 "C/C"     "T/T"     "T/C"     "C/C"     "T/A"    
## id25 "C/G"     "T/T"     "C/C"     "C/A"     "T/A"    
## id33 "C/C"     "T/T"     "T/C"     "C/C"     "T/A"    
## id35 "C/C"     "T/T"     "C/C"     "C/C"     "A/A"
```

How many IDs and SNPs do we have in our data?


```r
nids(ge03d2)
```

```
## [1] 950
```

```r
nsnps(ge03d2)
```

```
## [1] 7589
```

What information about SNPs do we have?


```r
snp_info <- summary(gtdata(ge03d2))

head(snp_info)
```

```
##           Chromosome Position Strand A1 A2 NoMeasured CallRate      Q.2
## rs1646456          1      653      +  C  G        938   0.9874 0.318230
## rs7950586          1      849      -  T  A        938   0.9874 0.036780
## rs4785242          1     1766      -  T  C        936   0.9853 0.627137
## rs4435802          1     5291      +  C  A        943   0.9926 0.080594
## rs2847446          1     5555      +  T  A        937   0.9863 0.317503
## rs9308393          1     6739      +  T  C        937   0.9863 0.002134
##           P.11 P.12 P.22  Pexact       Fmax    Plrt
## rs1646456  436  407   95 1.00000  0.0000406 0.99901
## rs7950586  873   61    4 0.03242  0.0821844 0.04059
## rs4785242  140  418  378 0.18400  0.0450984 0.16844
## rs4435802  799  136    8 0.37870  0.0268318 0.42778
## rs2847446  445  389  103 0.20028  0.0420746 0.19955
## rs9308393  933    4    0 1.00000 -0.0021390 0.92630
```

Description of columns:
  
  1. ...
  2. ...
  3. ...

What is the distribution of number of SNPs per chromosome?

```r
table(chromosome(ge03d2))
```

```
## 
##    1    2    3    X 
## 3555 1999 1466  569
```

## Summary for a single SNP

Let's extract information about particular SNP


```r
snp_info['rs6212914',]
```

```
##           Chromosome Position Strand A1 A2 NoMeasured CallRate     Q.2
## rs6212914          1  1099988      +  A  T        934   0.9832 0.02248
##           P.11 P.12 P.22 Pexact   Fmax   Plrt
## rs6212914  892   42    0      1 -0.023 0.3256
```

What is the distribution of genotypes for this SNP across samples


```r
table(as.character(gtdata(ge03d2[,'rs6212914'])))
```

```
## 
## A/A A/T 
## 892  42
```

Let's estiamte effective allele frequency for this SNP using different ways


```r
(snp_info['rs6212914','P.12'] + snp_info['rs6212914','P.22'] * 2) / snp_info['rs6212914','NoMeasured'] / 2
```

```
## [1] 0.02248
```

```r
mean(as.numeric(gtdata(ge03d2[,'rs6212914'])), na.rm=TRUE) / 2
```

```
## [1] 0.02248
```
Let's perform HWE test for this SNP


```r
HWE.show(ge03d2[,'rs6212914'])
```

```
## HWE summary for 1 :
##                A/A      A/B    B/B
## observed 8.920e+02 42.00000 0.0000
## expected 8.925e+02 41.05567 0.4722
## chi2+    2.498e-04  0.02172 0.4722
## Chi2 = 0.4941 ; P = 0.4821; exact P = 1
```

## Per-SNP Summary of genotype data

Let's have a look at the summary of the genotype data.

First, let's add EAF (effective allele frequency) column to the *snp_info* object


```r
snp_info$EAF <- (snp_info[,'P.12'] + snp_info[,'P.22'] * 2) / snp_info[,'NoMeasured'] / 2

head(snp_info)
```

```
##           Chromosome Position Strand A1 A2 NoMeasured CallRate      Q.2
## rs1646456          1      653      +  C  G        938   0.9874 0.318230
## rs7950586          1      849      -  T  A        938   0.9874 0.036780
## rs4785242          1     1766      -  T  C        936   0.9853 0.627137
## rs4435802          1     5291      +  C  A        943   0.9926 0.080594
## rs2847446          1     5555      +  T  A        937   0.9863 0.317503
## rs9308393          1     6739      +  T  C        937   0.9863 0.002134
##           P.11 P.12 P.22  Pexact       Fmax    Plrt      EAF
## rs1646456  436  407   95 1.00000  0.0000406 0.99901 0.318230
## rs7950586  873   61    4 0.03242  0.0821844 0.04059 0.036780
## rs4785242  140  418  378 0.18400  0.0450984 0.16844 0.627137
## rs4435802  799  136    8 0.37870  0.0268318 0.42778 0.080594
## rs2847446  445  389  103 0.20028  0.0420746 0.19955 0.317503
## rs9308393  933    4    0 1.00000 -0.0021390 0.92630 0.002134
```

Let's look at the distribution of EAF - effective allele frequency.


```r
hist(snp_info$EAF, breaks=50)
```

![](GWAS_practicals.Part_2_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

There are a lot of SNPs with EAF below 0.05 and above 0.95.

What about call rate?


```r
hist(snp_info$CallRate, breaks=50)
```

![](GWAS_practicals.Part_2_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

What about HWE test?


```r
catable(snp_info$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)
```

```
##      X<=6.58848333113717e-06 X<=0.01 X<=0.05   X<=0.1 all X
## No                   290.000 479.000 807.000 1113.000  7589
## Prop                   0.038   0.063   0.106    0.147     1
```

What if we will test HWE using controls?


```r
control_ids <- idnames(ge03d2)[phdata(ge03d2)$dm2==0]

catable(summary(gtdata(ge03d2[control_ids,]))$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)
```

```
##      X<=6.58848333113717e-06 X<=0.01 X<=0.05  X<=0.1 all X
## No                    41.000 113.000 321.000 543.000  7589
## Prop                   0.005   0.015   0.042   0.072     1
```
As you can see, now much less SNPs show deviation from HWE.

## Per-ID Summary of genotype data

Now let's check the per-ID SNP call rate and heterozygosity rate.

First let's get per-ID summary of the data.

```r
idsummary <- perid.summary(ge03d2)

head(idsummary)
```

```
##      NoMeasured NoPoly    Hom E(Hom)    Var         F CallPP    Het
## id4        7470   7469 0.7193 0.7419 0.4771 -0.087557 0.9843 0.2807
## id10       7454   7453 0.7509 0.7419 0.4723  0.034716 0.9822 0.2491
## id25       7467   7466 0.7348 0.7421 0.4875 -0.028142 0.9839 0.2652
## id33       7470   7469 0.7427 0.7420 0.4602  0.002546 0.9843 0.2573
## id35       7457   7456 0.7452 0.7422 0.4817  0.011761 0.9826 0.2548
## id58       7468   7467 0.7600 0.7419 0.4617  0.070203 0.9841 0.2400
```

Description of columns:

  1. ...
  2. ...
  3. ...
  
Let's check the distribution of per-ID call rate. It can give information about samples with low callrate, that may indicate problems with these samples.


```r
hist(idsummary$CallPP, breaks=100)
```

![](GWAS_practicals.Part_2_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

Few samples have low per-ID call rate.

What about heterozigosity rate?


```r
hist(idsummary$Het, breaks=100)
```

![](GWAS_practicals.Part_2_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

Few samples have outlying heterozigosity rate, that may indicate DNA sample contamination.

## Summary of exploring genotype data

Let's see the descriptive statistics for all SNP markes in our data.


```r
descriptives.marker(ge03d2)
```

```
## $`Minor allele frequency distribution`
##      X<=0.01 0.01<X<=0.05 0.05<X<=0.1 0.1<X<=0.2    X>0.2
## No   334.000     1335.000    1257.000   1699.000 2964.000
## Prop   0.044        0.176       0.166      0.224    0.391
## 
## $`Cumulative distr. of number of SNPs out of HWE, at different alpha`
##      X<=1e-04 X<=0.001 X<=0.01 X<=0.05 all X
## No    331.000  367.000 479.000 807.000  7589
## Prop    0.044    0.048   0.063   0.106     1
## 
## $`Distribution of proportion of successful genotypes (per person)`
##      X<=0.9 0.9<X<=0.95 0.95<X<=0.98 0.98<X<=0.99 X>0.99
## No    4.000           0        1.000      945.000      0
## Prop  0.004           0        0.001        0.995      0
## 
## $`Distribution of proportion of successful genotypes (per SNP)`
##      X<=0.9 0.9<X<=0.95 0.95<X<=0.98 0.98<X<=0.99   X>0.99
## No   82.000           0       95.000     4923.000 2489.000
## Prop  0.011           0        0.013        0.649    0.328
## 
## $`Mean heterozygosity for a SNP`
## [1] 0.2579
## 
## $`Standard deviation of the mean heterozygosity for a SNP`
## [1] 0.1607
## 
## $`Mean heterozygosity for a person`
## [1] 0.2458
## 
## $`Standard deviation of mean heterozygosity for a person`
## [1] 0.02703
```

```r
descriptives.marker(ge03d2, ids=control_ids)
```

```
## $`Minor allele frequency distribution`
##      X<=0.01 0.01<X<=0.05 0.05<X<=0.1 0.1<X<=0.2    X>0.2
## No   362.000     1354.000    1260.000    1672.00 2941.000
## Prop   0.048        0.178       0.166       0.22    0.388
## 
## $`Cumulative distr. of number of SNPs out of HWE, at different alpha`
##      X<=1e-04 X<=0.001 X<=0.01 X<=0.05 all X
## No     46.000   52.000 113.000 321.000  7589
## Prop    0.006    0.007   0.015   0.042     1
## 
## $`Distribution of proportion of successful genotypes (per person)`
##      X<=0.9 0.9<X<=0.95 0.95<X<=0.98 0.98<X<=0.99 X>0.99
## No    1.000           0        1.000      485.000      0
## Prop  0.002           0        0.002        0.996      0
## 
## $`Distribution of proportion of successful genotypes (per SNP)`
##      X<=0.9 0.9<X<=0.95 0.95<X<=0.98 0.98<X<=0.99   X>0.99
## No   82.000           0      251.000     4064.000 3192.000
## Prop  0.011           0        0.033        0.536    0.421
## 
## $`Mean heterozygosity for a SNP`
## [1] 0.2554
## 
## $`Standard deviation of the mean heterozygosity for a SNP`
## [1] 0.1619
## 
## $`Mean heterozygosity for a person`
## [1] 0.247
## 
## $`Standard deviation of mean heterozygosity for a person`
## [1] 0.02645
```


## Quality contron of genotype data

Now let's summarize our observations and perform quality control for genotype data.

Basic QC consists of following steps:

  1. Identification of individuals with discordant sex information
    1. X-Heterozygose males
    2. Y-Heterozygose females
  
  2. Identification of individuals with elevated missing data rates or outlying heterozygosity rate
    1. IDs with low call rate (<97%, 98% or 99%)
    2. IDs with ourlying heterozygosity rate
      1. Too high - DNA sample contamination 
      2. Too low - inbreeding
  
  3. Identification of duplicated or related individuals
    1. Estimation of kinship matrix
    2. IBS=0.5 -> monozygos twins or duplicated samples
    3. IBS > 0.125 -> close relatives
  
  4. Identification of individuals of divergent ancestry
    1. PCA analysis using kinship matrix
    2. 1000 Genomes project data can be used as reference sample
  
  5. Removal of all individuals failing QC
  
  6. Identification of all markers with an excessive missing data rate
    1. SNP call rate below 97%
  
  7. Removal of all markers failing QC
  

Let's perform QC using checkmarker function.


```r
qc1 <- check.marker(ge03d2, p.level = 0, maf = 0.01, perid.call = 0.97, callrate = 0.97)
```

```
## Excluding people/markers with extremely low call rate...
## 7589 markers and 950 people in total
## 0 people excluded because of call rate < 0.1 
## 7 markers excluded because of call rate < 0.1 
## Passed: 7582 markers and 950 people
## 
## Running sex chromosome checks...
## 1934 heterozygous X-linked male genotypes found
## 2 X-linked markers are likely to be autosomal (odds > 1000 )
## 10 male are likely to be female (odds > 1000 )
## 6 female are likely to be male (odds > 1000 )
## 0 people have intermediate X-chromosome inbreeding (0.5 > F > 0.5)
## If these people/markers are removed, 8 heterozygous male genotypes are left
## ... these will be considered missing in analysis.
## ... Use Xfix() to fix these problems.
## Passed: 7580 markers and 934 people
## 
## ... 8 X/Y/mtDNA ( 8 0 0 ) impossible heterozygotes and female Ys set as missing
## 
## 
## RUN 1 
## 7580 markers and 934 people in total
## 360 (4.749%) markers excluded as having low (<1%) minor allele frequency
## 75 (0.9894%) markers excluded because of low (<97%) call rate
## 0 (0%) markers excluded because they are out of HWE (P <0)
## 4 (0.4283%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.2657 (s.e. 0.02132)
## 4 (0.4283%) people excluded because too high autosomal heterozygosity (FDR <1%)
## Excluded people had HET >= 0.4789
## Mean IBS is 0.7761 (s.e. 0.01639), as based on 2000 autosomal markers
## 8 (0.8565%) people excluded because of too high IBS (>=0.95)
## In total, 7145 (94.26%) markers passed all criteria
## In total, 918 (98.29%) people passed all criteria
## 
## RUN 2 
## 7145 markers and 918 people in total
## 34 (0.4759%) markers excluded as having low (<1%) minor allele frequency
## 0 (0%) markers excluded because of low (<97%) call rate
## 0 (0%) markers excluded because they are out of HWE (P <0)
## 0 (0%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.266 (s.e. 0.01573)
## 0 people excluded because too high autosomal heterozygosity (FDR <1%)
## Mean IBS is 0.7772 (s.e. 0.01765), as based on 2000 autosomal markers
## 0 (0%) people excluded because of too high IBS (>=0.95)
## In total, 7111 (99.52%) markers passed all criteria
## In total, 918 (100%) people passed all criteria
## 
## RUN 3 
## 7111 markers and 918 people in total
## 0 (0%) markers excluded as having low (<1%) minor allele frequency
## 0 (0%) markers excluded because of low (<97%) call rate
## 0 (0%) markers excluded because they are out of HWE (P <0)
## 0 (0%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.266 (s.e. 0.01573)
## 0 people excluded because too high autosomal heterozygosity (FDR <1%)
## Mean IBS is 0.7754 (s.e. 0.0166), as based on 2000 autosomal markers
## 0 (0%) people excluded because of too high IBS (>=0.95)
## In total, 7111 (100%) markers passed all criteria
## In total, 918 (100%) people passed all criteria
```

Now let's define IDs which will be used for HWE test.


```r
hweids <- intersect(control_ids, qc1$idok)

qc2 <- check.marker(ge03d2, 
                    p.level = 1e-6, 
                    maf = 0.01, 
                    perid.call = 0.97, 
                    callrate = 0.97,
                    hweidsubset = hweids)
```

```
## Excluding people/markers with extremely low call rate...
## 7589 markers and 950 people in total
## 0 people excluded because of call rate < 0.1 
## 7 markers excluded because of call rate < 0.1 
## Passed: 7582 markers and 950 people
## 
## Running sex chromosome checks...
## 1934 heterozygous X-linked male genotypes found
## 2 X-linked markers are likely to be autosomal (odds > 1000 )
## 10 male are likely to be female (odds > 1000 )
## 6 female are likely to be male (odds > 1000 )
## 0 people have intermediate X-chromosome inbreeding (0.5 > F > 0.5)
## If these people/markers are removed, 8 heterozygous male genotypes are left
## ... these will be considered missing in analysis.
## ... Use Xfix() to fix these problems.
## Passed: 7580 markers and 934 people
## 
## ... 8 X/Y/mtDNA ( 8 0 0 ) impossible heterozygotes and female Ys set as missing
## 
## 
## RUN 1 
## 7580 markers and 934 people in total
## 360 (4.749%) markers excluded as having low (<1%) minor allele frequency
## 75 (0.9894%) markers excluded because of low (<97%) call rate
## 37 (0.4881%) markers excluded because they are out of HWE (P <1e-06)
## 4 (0.4283%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.2657 (s.e. 0.02132)
## 4 (0.4283%) people excluded because too high autosomal heterozygosity (FDR <1%)
## Excluded people had HET >= 0.4789
## Mean IBS is 0.7785 (s.e. 0.01697), as based on 2000 autosomal markers
## 8 (0.8565%) people excluded because of too high IBS (>=0.95)
## In total, 7145 (94.26%) markers passed all criteria
## In total, 918 (98.29%) people passed all criteria
## 
## RUN 2 
## 7145 markers and 918 people in total
## 34 (0.4759%) markers excluded as having low (<1%) minor allele frequency
## 0 (0%) markers excluded because of low (<97%) call rate
## 0 (0%) markers excluded because they are out of HWE (P <1e-06)
## 0 (0%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.266 (s.e. 0.01573)
## 0 people excluded because too high autosomal heterozygosity (FDR <1%)
## Mean IBS is 0.7779 (s.e. 0.01602), as based on 2000 autosomal markers
## 0 (0%) people excluded because of too high IBS (>=0.95)
## In total, 7111 (99.52%) markers passed all criteria
## In total, 918 (100%) people passed all criteria
## 
## RUN 3 
## 7111 markers and 918 people in total
## 0 (0%) markers excluded as having low (<1%) minor allele frequency
## 0 (0%) markers excluded because of low (<97%) call rate
## 0 (0%) markers excluded because they are out of HWE (P <1e-06)
## 0 (0%) people excluded because of low (<97%) call rate
## Mean autosomal HET is 0.266 (s.e. 0.01573)
## 0 people excluded because too high autosomal heterozygosity (FDR <1%)
## Mean IBS is 0.7794 (s.e. 0.0162), as based on 2000 autosomal markers
## 0 (0%) people excluded because of too high IBS (>=0.95)
## In total, 7111 (100%) markers passed all criteria
## In total, 918 (100%) people passed all criteria
```

QC is done. Now let's save QCed data


```r
ge03d2clean <- ge03d2[qc2$idok,qc2$snpok]
```


## Additional literature: 

Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738

Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564–73. http://doi.org/10.1038/nprot.2010.116
## Saving data
