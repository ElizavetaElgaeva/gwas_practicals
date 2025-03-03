#title: "��������� �� ������������ �������������� ��������. ����� 3: ���������� GWAS."
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

load('.../ge03d2clean.RData')

# ��������

nids(ge03d2clean)
nsnps(ge03d2clean)

## Single SNP association analysis

# �������� ���������� ��� ������ SNP
# �������� ������ ��������� � ���������

snpnames(ge03d2clean)[1:10]
my_snp_num <- as.numeric(gtdata(ge03d2clean[,'rs1646456']))
my_snp_cha <- as.character(gtdata(ge03d2clean[,'rs1646456']))

# �������� ������������� ��������� �� ��������

table(my_snp_cha)
table(my_snp_num)

# �������� ������������ ��������� � SNP � �������� ����������� ��������� � dm2 � �������� ���������

log_reg_snp <- glm(phdata(ge03d2clean)$dm2~my_snp_num)
summary(log_reg_snp)

# ���������� �� ���������� (p-value>0.05).
# �������� ������������� ������� � ��������� �� ���������

table(phdata(ge03d2clean)$dm2,my_snp_cha)

# ������������� � ���������������� ��������

## Genome-wide association scan

# �������� ��������� GWAS, ��������� ������� *mlreg* �� ������ GenABEL
# ����������, �� �������� 7000 ���������. �� �������� GWAS ��� dm2 � ����� � ��������� � �������� ��������

qt_simple <- mlreg(dm2~sex+age, data=ge03d2clean, trait="binomial")

# �������� �� ���-10 ��������

descriptives.scan(qt_simple , sort="Pc1df")

# ��� ����������, ���� �� ������� BMI ��� ���������?

qt_sex_age_bmi <- mlreg(dm2 ~ sex + age + bmi, ge03d2clean, trait="binomial")
descriptives.scan(qt_sex_age_bmi , sort="Pc1df")

## Manhattan plot

plot(qt_sex_age_bmi,df = 1)

## Regional assoc plot

gwas_sum <- summary(qt_sex_age_bmi)
gwas_sum_sm <- gwas_sum[gwas_sum$Chromosome==1 & gwas_sum$Position>1047389-250000 & gwas_sum$Position<1047389+250000,]
plot(-log10(gwas_sum_sm$P1df)~gwas_sum_sm$Position)

## �������� ������������� ��������� � �������/��������� ��� ���-SNP

table(phdata(ge03d2clean)$dm2,as.character(gtdata(ge03d2clean[,'rs7903146'])))

# ������ ������� �������. �������� OR ��� ������ T � ������� ������� exp()


# ������ ������� ���������� mlreg � glm

summary(glm(dm2~sex+age+bmi+as.numeric(gtdata(ge03d2clean[,'rs7903146'])),family='binomial', data=phdata(ge03d2clean)))


# ��� �����, ���������� ����������� ���� � ������

## �������� ��������� ����������

sum_stats <- results(qt_sex_age_bmi)
snp_info <- summary(gtdata(ge03d2clean))
sum_stats2write <- sum_stats[,c('Chromosome','Position','A2','A1','N','effB','se_effB','P1df','Strand')]
sum_stats2write <- cbind(rownames(snp_info),rownames(snp_info),sum_stats2write)
sum_stats2write <- cbind(sum_stats2write[,c(1:7)],snp_info$Q.2,sum_stats2write[,c(8:11)])
colnames(sum_stats2write) <- c('rsid','snpid','chr','pos','a1','a0','n','freq1','beta1','se','p','strand')
write.table(sum_stats2write,file='.../sum_stats.tsv', sep='\t', quote = FALSE, row.names = FALSE)

## ����������: 

# Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738

# Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564�73. http://doi.org/10.1038/nprot.2010.116
