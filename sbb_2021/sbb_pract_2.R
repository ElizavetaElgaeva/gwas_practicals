#title: "��������� �� ������������ �������������� ��������. ����� 2: ������ � ������� ���������, �������� ��������."
#author: "����� �������, ������� ���������"
#date: "6 ������� 2021 �."

# ���� ����������

#������������ ������� ������� �� ���� ������:

#1. ������ � ��������������� �������
#2. ������ � ������� ���������, �������� ��������
#3. ���������� GWAS

## �������� ����������� �������

library(GenABEL)
library(GenABEL.data)

## �������� ��������� ������ ������

data("ge03d2")

# ����� ��������� ������

str(ge03d2)
gtdata(ge03d2) #  �� �� �����, ��� � summary(ge03d2@gtdata)
ge03d2@gtdata@chromosome

# ������ � ����������

nids(ge03d2)
nsnps(ge03d2)
chromosome(ge03d2)
coding(ge03d2)
as.character(gtdata(ge03d2[1:5, 1:5]))
as.numeric(gtdata(ge03d2[1:5, 1:5]))
effallele(ge03d2@gtdata)
refallele(ge03d2@gtdata)

## ���������� �� SNP
  
snp_info <- summary(gtdata(ge03d2))
head(snp_info)

## �������� �������:
  
# 1. Q.2 - frequency of allele 2 (A2)
# 2. Pexact - p-value of HWE test
# 3. P.11, P.12, P.22 - number of genotypes

## ������������� ���������� SNP �� ����������
  
table(chromosome(ge03d2))

## ���������� �� ��������� SNP

snp_info['rs6212914',]

## ������������� ��������� ��� ���������� SNP

table(as.character(gtdata(ge03d2[,'rs6212914'])))


## ������ ������� ������������ ������

(snp_info['rs6212914','P.12'] + snp_info['rs6212914','P.22'] * 2) / snp_info['rs6212914','NoMeasured'] / 2
mean(as.numeric(gtdata(ge03d2[,'rs6212914'])), na.rm=TRUE) / 2

## �������� �� ���������� �����-��������� (HWE test)

HWE.show(ge03d2[,'rs6212914'])

## ���������� �� ���������� SNP

# ������ ������� ������� EAF (effective allele frequency). �� �� ���������� ���� � 'Q.2',
# �� �� ��������� ��������� �� ��������������

snp_info$EAF <- (snp_info[,'P.12'] + snp_info[,'P.22'] * 2) / snp_info[,'NoMeasured'] / 2
head(snp_info)

# ������������� EAF

hist(snp_info$EAF, breaks=50)

# ��������� SNP c EAF ������ 0.05

# �������� ������� ��������������� �� �������� (call rate)
  
hist(snp_info$CallRate, breaks=50)

# �������� HWE test
  
catable(snp_info$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)

# �������� ��� ����������� �������
  
control_ids <- idnames(ge03d2)[phdata(ge03d2)$dm2==0]
catable(summary(gtdata(ge03d2[control_ids,]))$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)

# � ����������� ������� ������ SNP ����������� �� ����������


## ������ ���������, ��� ��������� ���������� �� ������

# ���������� �� ID
idsummary <- perid.summary(ge03d2)
head(idsummary)

# �������� �������:
  
# 1. NoMeasured - number of measured SNPs
# 2. NoPoly - number of polymorfic SNPs (SNPs with minor allele frequency > 0)
# 3. Het - heterozigosity rate

# �������� ������� ��������������� �� SNP

hist(idsummary$CallPP, breaks=100)

# �������� ���������������� (heterozigosity rate)

hist(idsummary$Het, breaks=100)

# ��������� �������� � ������� �����������������, ��� ������� � ��������� ������������ ��������

descriptives.marker(ge03d2)
descriptives.marker(ge03d2, ids=control_ids)

## �������� �������� ������ ���������������

# ������� ����� QC:

# 1. �������� �������������� �� ����
# a. X-�������������� �������
# b. Y-�������������� �������
 
# 2. �������� �� �������� ��������������� �� SNP � ����� ��������� �� ����������������
# a. ID � ������ call rate (<95%, 97% or 99%)
# b. ID � �������������� ���������� ����������������
#  ������� - ������������ ��� 
#  ������ - ���������

# 3. ����� ���������� � �������������
# a. ������ ������� ������� (kinship matrix)
#  IBS = 1 -> ������������ �������� ��� ���������
#  IBS > 0.125 -> ������� ������������

# 4. ����������� �������� ������ ���������� ��������������
# a. PCA ������ � �������������� ������� �������
# b. 1000 Genomes project � �������� ���������

# 5. �������� ���� ��������, �� ��������� QC
 
# 6. �������� SNP �� ������� ���������������� �� ��������
# a. SNP � call rate ������ 97% (95%, 99%)

# 7. �������� ������ �������
# MAF > 0.01

# 8. �������� �� ���������� �����-���������
 
# 9. �������� ���� SNP, �� ��������� QC


# ������ QC c �������� **check.marker**

qc1 <- check.marker(ge03d2, p.level = 0, maf = 0.01, perid.call = 0.97, callrate = 0.97)

# ��������� ID, ������� ����� ��������� �� HWE

hweids <- intersect(control_ids, qc1$idok)
qc2 <- check.marker(ge03d2, 
p.level = 1e-6, 
maf = 0.01, 
perid.call = 0.97, 
callrate = 0.97,
hweidsubset = hweids)

# QC �������, �������� ���������

ge03d2clean <- ge03d2[qc2$idok,qc2$snpok]

## �������� ������ �� �������������e

# ���������� ��������, ������������ Price at al.
# ������� �������� ������� ������� �� ���������� ��������:
  
gkin <- ibs(ge03d2clean[, autosomal(ge03d2clean)], weight="freq")
gkin[1:5, 1:5]

# ����� ��� ���������� �������� ������ ��������� �������
# (�genomic kinship� ��� �genome-wide IBS�), ����� �� ��������� �������������
# 0.5 + �������������� ������, ��������������� �������� ������� � ���,
# ������� SNP ���� �������������� ��� ����� �������� (IBD ����������� �� ���� SNP)

# ��������� ������� � ������� ����������

data.dist <- as.dist(0.5 - gkin)

# ����������� ������������ ������������� (Classical Multidimensional Scaling)

data.mds <- cmdscale(data.dist)

# ������������ ����������

plot(data.mds)

# ������ ����� �� ������� - ��� ��������� �������, 2D ���������� ����� ���� 
# ���������� ���, ����� ���� ����������� ������� � ���������� �� �������� ������� ����������
# ������ ���� �������������� � ��� ������ 
# ������� id �������� �� ������� �������� � �������� ��, ����� �������� QC

clean_ids <- rownames(data.mds)[data.mds[,1]<0.1]
controls_clean <- intersect(clean_ids,idnames(ge03d2clean)[phdata(ge03d2clean)$dm2==0])
qc3 <- check.marker(ge03d2[clean_ids,], 
p.level = 1e-6, 
maf = 0.01, 
perid.call = 0.97, 
callrate = 0.97,
hweidsubset = controls_clean)
ge03d2clean <- ge03d2[qc3$idok,qc3$snpok]


## �������� ������ ��� GWAS

save(ge03d2clean, file='.../ge03d2clean.RData')

## ����������: 

# Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738
 
# Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564�73. http://doi.org/10.1038/nprot.2010.116

