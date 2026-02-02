####----Project:1.Alpha Diversity----####
####----Author:LiuYueR----####
####----Date:20240131----####
####----Email:yueliu1115@163.com----####


rm(list = ls())

####----load R Package ----####
library(tidyverse)
library(vegan)
library(readxl)
library(Biostrings)
library(ggpubr)
library(ggfun)
library(patchwork)
library(ggpmisc)
library(ggrepel)
library(cowplot)
library(picante)
library(ape)
library(rstatix)

set.seed(1115)

dir.create("1.Alpha_Diversity")

####----O Layer ----####
#####-----load Data-----#####
# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean.txt", header = T, row.names = 1, 
                     sep = "\t", check.names = F,  comment.char = "")
seqtab <- seqtab[which(rowSums(seqtab) != 0),]
dim(seqtab)
colSums(seqtab)
colSums(seqtab) %>% min()

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1)
dim(taxa)

# rare OTU table

otu_rare = as.data.frame(t(rrarefy(t(seqtab), 22000)))
colSums(otu_rare)

####----PD----####
phy_tree <- read.tree(file = "Input_File/otus.tree")
phy_tree

pd_results <- pd(samp = t(otu_rare), tree = phy_tree, include.root = F)

PD <- pd_results

#####-----Alpha diversity-----#####
otutab_rare_t <- t(otu_rare)

# Richness/ Observed species
observed_species <- estimateR(otutab_rare_t)[1, ]

# Chao 1
Chao1 <- estimateR(otutab_rare_t)[2, ]

# ACE
ACE <- estimateR(otutab_rare_t)[4, ]

# Shannon
Shannon_1 <- vegan::diversity(otutab_rare_t, index = 'shannon', base = exp(1))
Shannon_2 <- vegan::diversity(otutab_rare_t, index = 'shannon', base = 2)

# Simpson
# Gini-Simpson 指数代码 我们常用的Simpson 指数
Simpson_Gini <- vegan::diversity(otutab_rare_t, index = "simpson")

# 经典 Simpson 指数（使用频率比较低）
Simpson_Class <- 1- Simpson_Gini

# Inverse Simpson
invsimpson <- 1 / Simpson_Gini

# goods_coverages指数
goods_coverage <- 1 - rowSums(otutab_rare_t == 1) / rowSums(otutab_rare_t)
goods_coverage

# alpha 多样性指数
alpha_result <- data.frame(
  observed_species,
  Chao1,
  ACE,
  Shannon_1,
  Shannon_2,
  Simpson_Gini,
  Simpson_Class,
  invsimpson,
  goods_coverage,
  PD
)

row.names_vector <- row.names(alpha_result)

alpha_result$Sample <- row.names(alpha_result)

alpha_result$Group <- str_remove(alpha_result$Sample, pattern = "_\\d+")
alpha_result$Sample_month <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,1]
alpha_result$Treatment <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,2]
alpha_result$Sample_year <- str_split(alpha_result$Group, pattern = "\\.", simplify = T)[,1]

alpha_result2 <- alpha_result %>%
  dplyr::filter(!str_detect(string = Sample_month, pattern = "2019.7|2020.11"))

write.csv(alpha_result2, file = "./1.Alpha_Diversity/alpha_result_O_Layer_for_Analysis.csv")


####----0-10cm Layer ----####
#####-----load Data-----#####
# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean_0_10cm.txt", header = T, row.names = 1,
                     sep = "\t", check.names = F,  comment.char = "")
seqtab <- seqtab[which(rowSums(seqtab) != 0),]
dim(seqtab)
colSums(seqtab) %>% min()

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1)
dim(taxa)

# rare OTU table
set.seed(1115)
# otu_rare = as.data.frame(t(rrarefy(t(seqtab), round(min(colSums(seqtab)), -4))))
otu_rare = as.data.frame(t(rrarefy(t(seqtab), 25000)))
colSums(otu_rare)

#####-----Alpha diversity-----#####
otutab_rare_t <- t(otu_rare)

# Richness/ Observed species
observed_species <- estimateR(otutab_rare_t)[1, ]

# Chao 1
Chao1 <- estimateR(otutab_rare_t)[2, ]

# ACE
ACE <- estimateR(otutab_rare_t)[4, ]

# Shannon
Shannon_1 <- vegan::diversity(otutab_rare_t, index = 'shannon', base = exp(1))
Shannon_2 <- vegan::diversity(otutab_rare_t, index = 'shannon', base = 2)

# Simpson
# Gini-Simpson 指数代码
Simpson_Gini <- vegan::diversity(otutab_rare_t, index = "simpson")

Simpson_Class <- 1- Simpson_Gini

# goods_coverages指数
goods_coverage <- 1 - rowSums(otutab_rare_t == 1) / rowSums(otutab_rare_t)
goods_coverage

# alpha 多样性指数
alpha_result <- data.frame(
  observed_species,
  Chao1,
  ACE,
  Shannon_1,
  Shannon_2,
  Simpson_Gini,
  Simpson_Class,
  goods_coverage
)

row.names_vector <- row.names(alpha_result)
alpha_result$Sample <- row.names(alpha_result)
alpha_result$Group <- str_remove(alpha_result$Sample, pattern = "_\\d+")
alpha_result$Sample_month <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,1]
alpha_result$Treatment <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,2]
alpha_result$Sample_year <- str_split(alpha_result$Group, pattern = "\\.", simplify = T)[,1]

alpha_result2 <- alpha_result %>%
  dplyr::filter(!str_detect(string = Sample_month, pattern = "2019.7|2020.11"))

write.csv(alpha_result2, file = "./1.Alpha_Diversity/alpha_result_0_10cm_Layer_for_Analysis.csv")