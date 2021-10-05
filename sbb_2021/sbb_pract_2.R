#title: "Практикум по Компьютерной статистической геномике. Часть 2: Работа с данными генотипов, контроль качества."
#author: "Содбо Шарапов, Елгаева Елизавета"
#date: "6 октября 2021 г."

# План практикума

#Практическое занятие состоит из трех частей:

#1. Работа с фенотипическими данными
#2. Работа с данными генотипов, контроль качества
#3. Проведение GWAS

## Загрузка программных пакетов

library(GenABEL)
library(GenABEL.data)

## Загрузка тестового набора данных

data("ge03d2")

# Общая структура данных

str(ge03d2)
gtdata(ge03d2) #  то же самое, что и summary(ge03d2@gtdata)
ge03d2@gtdata@chromosome

# Работа с генотипами

nids(ge03d2)
nsnps(ge03d2)
chromosome(ge03d2)
coding(ge03d2)
as.character(gtdata(ge03d2[1:5, 1:5]))
as.numeric(gtdata(ge03d2[1:5, 1:5]))
effallele(ge03d2@gtdata)
refallele(ge03d2@gtdata)

## Информация по SNP
  
snp_info <- summary(gtdata(ge03d2))
head(snp_info)

## Описание колонок:
  
# 1. Q.2 - frequency of allele 2 (A2)
# 2. Pexact - p-value of HWE test
# 3. P.11, P.12, P.22 - number of genotypes

## Распределение количества SNP по хромосомам
  
table(chromosome(ge03d2))

## Информация по отдельным SNP

snp_info['rs6212914',]

## Распределение генотипов для отдельного SNP

table(as.character(gtdata(ge03d2[,'rs6212914'])))


## Оценка частоты эффекторного аллеля

(snp_info['rs6212914','P.12'] + snp_info['rs6212914','P.22'] * 2) / snp_info['rs6212914','NoMeasured'] / 2
mean(as.numeric(gtdata(ge03d2[,'rs6212914'])), na.rm=TRUE) / 2

## Проверка на равновесие Харди-Вайнберга (HWE test)

HWE.show(ge03d2[,'rs6212914'])

## Статистика по отдельному SNP

# Сперва добавим колонку EAF (effective allele frequency). Та же информация есть в 'Q.2',
# но мы попробуем посчитать ее самостоятельно

snp_info$EAF <- (snp_info[,'P.12'] + snp_info[,'P.22'] * 2) / snp_info[,'NoMeasured'] / 2
head(snp_info)

# Распределение EAF

hist(snp_info$EAF, breaks=50)

# Множество SNP c EAF меньше 0.05

# Проверим процент генотипирования по образцам (call rate)
  
hist(snp_info$CallRate, breaks=50)

# Проведем HWE test
  
catable(snp_info$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)

# Отдельно для контрольной выборки
  
control_ids <- idnames(ge03d2)[phdata(ge03d2)$dm2==0]
catable(summary(gtdata(ge03d2[control_ids,]))$Pexact, c(0.05/nsnps(ge03d2), 0.01, 0.05, 0.1), cum=TRUE)

# В контрольной выборке меньше SNP отклоняются от равновесия


## Теперь посмотрим, как проверить статистику по обазцу

# Статистика по ID
idsummary <- perid.summary(ge03d2)
head(idsummary)

# Описание колонок:
  
# 1. NoMeasured - number of measured SNPs
# 2. NoPoly - number of polymorfic SNPs (SNPs with minor allele frequency > 0)
# 3. Het - heterozigosity rate

# Проверим процент генотипирования по SNP

hist(idsummary$CallPP, breaks=100)

# Проверим гетерозиготность (heterozigosity rate)

hist(idsummary$Het, breaks=100)

# Несколько образцов с высокой гетерозиготностью, что говорит о возможной контаминации образцов

descriptives.marker(ge03d2)
descriptives.marker(ge03d2, ids=control_ids)

## Контроль качества данных генотипирования

# Базовые этапы QC:

# 1. Проверка несоответствий по полу
# a. X-гетерозиготные мужчины
# b. Y-гетерозиготные женщины
 
# 2. Проверка по проценту генотипирования по SNP и поиск аутлаеров по гетерозиготности
# a. ID с низким call rate (<95%, 97% or 99%)
# b. ID с экстримальными значениями гетерозиготности
#  высокие - контаминация ДНК 
#  низкие - инбридинг

# 3. Поиск дубликатов и родственников
# a. Оценка матрицы родства (kinship matrix)
#  IBS = 1 -> монозиготные близнецы или дубликаты
#  IBS > 0.125 -> близкие родственники

# 4. Определение образцов другой этнической принадлежности
# a. PCA анализ с использованием матрицы родства
# b. 1000 Genomes project в качестве референса

# 5. Удаление всех образцов, не прошедших QC
 
# 6. Проверка SNP на процент генотиппирования по образцам
# a. SNP с call rate меньше 97% (95%, 99%)

# 7. Проверка частот аллелей
# MAF > 0.01

# 8. Проверка на равновесие Харди-Вайнберга
 
# 9. Удаление всех SNP, не прошедших QC


# Начнем QC c функцией **check.marker**

qc1 <- check.marker(ge03d2, p.level = 0, maf = 0.01, perid.call = 0.97, callrate = 0.97)

# Определим ID, которые будут проверены на HWE

hweids <- intersect(control_ids, qc1$idok)
qc2 <- check.marker(ge03d2, 
p.level = 1e-6, 
maf = 0.01, 
perid.call = 0.97, 
callrate = 0.97,
hweidsubset = hweids)

# QC окончен, сохраним результат

ge03d2clean <- ge03d2[qc2$idok,qc2$snpok]

## Проверим данные на стартификациюe

# Используем алгоритм, предложенный Price at al.
# Сначала вычислим матрицу родства по аутосомным маркерам:
  
gkin <- ibs(ge03d2clean[, autosomal(ge03d2clean)], weight="freq")
gkin[1:5, 1:5]

# Числа под диагональю отражают оценки геномного родства
# (’genomic kinship’ или ’genome-wide IBS’), сисла по диагонали соответствуют
# 0.5 + гомозиготность генома, наддиагональные элементы говорят о том,
# сколько SNP было генотипировано для обоих объектов (IBD оценивается по этим SNP)

# Переведем матрицу в матрицу расстояний

data.dist <- as.dist(0.5 - gkin)

# Многомерное шкалирование коэффициентов (Classical Multidimensional Scaling)

data.mds <- cmdscale(data.dist)

# Визуализация результата

plot(data.mds)

# Каждая точка на графике - это отдельный индивид, 2D расстояние между ними 
# подстроено так, чтобы быть максимально близким к расстоянию из исходной матрицы расстояний
# Данные явно кластеризуются в две группы 
# Получим id образцов из второго кластера и исключим их, затем повторим QC

clean_ids <- rownames(data.mds)[data.mds[,1]<0.1]
controls_clean <- intersect(clean_ids,idnames(ge03d2clean)[phdata(ge03d2clean)$dm2==0])
qc3 <- check.marker(ge03d2[clean_ids,], 
p.level = 1e-6, 
maf = 0.01, 
perid.call = 0.97, 
callrate = 0.97,
hweidsubset = controls_clean)
ge03d2clean <- ge03d2[qc3$idok,qc3$snpok]


## Сохраним данные для GWAS

save(ge03d2clean, file='.../ge03d2clean.RData')

## Литература: 

# Aulchenko, Yurii S., Karssen, Lennart C., & The GenABEL project developers. (2015). The GenABEL Tutorial. Zenodo. http://doi.org/10.5281/zenodo.19738
 
# Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case-control association studies. Nature Protocols, 5(9), 1564–73. http://doi.org/10.1038/nprot.2010.116

