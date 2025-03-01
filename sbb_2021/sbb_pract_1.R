#title: "��������� �� ������������ �������������� ��������. ����� 1: ������ � ��������������� �������."
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
phdata(ge03d2) # �� �� �����, ��� � ge03d2@phdata
ge03d2@phdata$dm2

gtdata(ge03d2) #  �� �� �����, ��� � summary(ge03d2@gtdata)
ge03d2@gtdata@chromosome

# ������ � ������� ���������

pheno <- phdata(ge03d2)
head(pheno)


## �������� ������� � pheno:

# 1. ID: id code of sample 
# 2. sex: sex, '0' is female, '1' is male
# 3. age: age, in years
# 4. dm2: diabetes mellitus type 2, '0' is unaffected, '1' is affected
# 5. height: height, in cm
# 6. wieght: wieght, in kg
# 7. diet: following the diet, '0' is false, '1' is true
# 8. bmi: BMI, body mass index, in $cm/kg^2$

## C�������� �������������� ������

summary(pheno)

## ��� �������������� ��������� �������� ������� **table**

table(pheno$sex)
table(pheno$dm2)
table(pheno$diet)

## �������� ������������� ��������� ��� �������������� ���������

hist(pheno$age, main='Histogram of age', xlab = 'Age')
hist(pheno$height, main='Histogram of height', xlab = 'Height')
hist(pheno$weight, main='Histogram of weight', xlab = 'Weight')
hist(pheno$bmi, main='Histogram of bmi', xlab = 'BMI')

## �������� �� ������������� ������������� � ������� ������ �����
## ������� �������� - ������������� �������� ����������.
## ���� p-value <= 0.05/n, ��� n - ���������� ������, �� ������������� ������������
 
shapiro.test(pheno$height)
shapiro.test(pheno$weight)
shapiro.test(pheno$bmi)

## ��� ���������� � ����������� ������������� ����� ������������ ����������������

hist(log(pheno$bmi), main='Histogram of logBMI', xlab = 'logBMI')
shapiro.test(log(pheno$bmi))

# ��� �������� �������������

hist(rntransform(pheno$bmi), main='Histogram of BMI after rank transformation', xlab = 'BMI, rank transformed')
shapiro.test(rntransform(pheno$bmi))

## ������ �������������� ��������� ��� ����������� �������� 
## ����� ������������ ����������� ���������� ������� �/��� �������� ������ �������� ���������

cor_matrix <- cor(pheno[,c('dm2','sex','age','height','weight','bmi','diet')], use='complete.obs')
cor_matrix

## ������������ � ������� ������ *corplot*

library(corrplot)
corrplot(cor_matrix)

## �������� �������� ��������� �� ���������

dm2_weight <- glm(dm2~weight,data = pheno,family = 'binomial')
summary(dm2_weight)

dm2_bmi <- glm(dm2~bmi,data = pheno,family = 'binomial')
summary(dm2_bmi)

dm2_sex <- glm(dm2~sex,data = pheno,family = 'binomial')
summary(dm2_sex)

# beta ��� ���� 0.4465, ��� ������������� OR=exp(0.4465)=1.563. �������� �� ������� 2x2:
table(pheno$dm2,pheno$sex)
(275*253)/(210*212)

dm2_cov <- glm(dm2~sex + age,data = pheno,family = 'binomial')
summary(dm2_cov)

## ���������� �������������� ������ ����� �� ��������� ���:

table(is.na(pheno))
table(is.na(pheno$id))
table(is.na(pheno$sex))
table(is.na(pheno$age))
table(is.na(pheno$dm2))
table(is.na(pheno$height))
table(is.na(pheno$weight))
table(is.na(pheno$diet))
table(is.na(pheno$bmi))

dm2_adj <- residuals(glm(dm2~sex + age,data = pheno,family = 'binomial'))
shapiro.test(dm2_adj)
shapiro.test(rntransform(dm2_adj))
dm2_adj_norm <- rntransform(dm2_adj)

