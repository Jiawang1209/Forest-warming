####----Project:4.Lefse_Biomaker_DES----####
####----Author:LiuYueR----####
####----Date:20240131----####
####----Email:yueliu1115@163.com----####


rm(list = ls())

####----load R Package ----####
library(tidyverse)
library(patchwork)
library(pheatmap)
library(indicspecies)
library(DESeq2)
library(magrittr)
library(microeco)
library(vegan)
library(ggfun)
library(aplot)
require(scales)
library(ggstar)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)
library(edgeR)
library(phyloseq)
library(ggvenn)
source("R/5_1.DESeq2_function.R")
source("R/5_2.edgeR_function.R")
source("R/5_3.LefSe_function.R")
source("R/5_4.STAMP_function.R")



set.seed(1115)
dir.create("5.Biomaker")

####----O Layer ----####
#####----load Data ----#####
# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean.txt", header = T, row.names = 1, 
                     sep = "\t", check.names = F,  comment.char = "")
seqtab <- seqtab[which(rowSums(seqtab) != 0),]
dim(seqtab)

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1)
dim(taxa)

# rare OTU table
set.seed(1115)
# otu_rare = as.data.frame(t(rrarefy(t(seqtab), round(min(colSums(seqtab)), -4))))
otu_rare = as.data.frame(t(rrarefy(t(seqtab), 22000)))

colSums(otu_rare)
colSums(seqtab)


otu_rare_relative <- otu_rare / colSums(otu_rare)
colSums(otu_rare_relative)
dim(otu_rare_relative)


# 
Group_info <- data.frame(
  sample = colnames(otu_rare_relative),
  group = str_remove(string = colnames(otu_rare_relative),  pattern = "_\\d+"),
  times = str_split(string = colnames(otu_rare_relative), pattern = "_", simplify = T)[,1]
)

Group_info

table(Group_info$times) %>% names()


#####-----edgeR-----#####
for (time in table(Group_info$times) %>% names() %>%.[-1]) {
  # time = '2019.6'
  print(time)
  sample <- paste0(paste0(rep(time, 12), rep(c("_Control_", "_Warming_"), each=6)), 1:6)
  sample2 <- paste0(time, "_Warming_vs_Control")
  edgeR_function(data = "seqtab",
                 sample1 = sample[1], sample2 = sample[2], sample3 = sample[3], 
                 sample4 = sample[4], sample5 = sample[5], sample6 = sample[6],
                 sample7 = sample[7], sample8 = sample[8], sample9 = sample[9], 
                 sample10 = sample[10], sample11 = sample[11],sample12 = sample[12], 
                 outname = sample2, path = "5.Biomaker")
}

#####-----DESeq2-----#####
for (time in table(Group_info$times) %>% names() %>%.[-1]) {
  # time = '2019.6'
  print(time)
  sample <- paste0(paste0(rep(time, 12), rep(c("_Control_", "_Warming_"), each=6)), 1:6)
  sample2 <- paste0(time, "_Warming_vs_Control")
  DESeq2_function(data = "seqtab",
                  sample1 = sample[1], sample2 = sample[2], sample3 = sample[3], 
                  sample4 = sample[4], sample5 = sample[5], sample6 = sample[6],
                  sample7 = sample[7], sample8 = sample[8], sample9 = sample[9], 
                  sample10 = sample[10], sample11 = sample[11],sample12 = sample[12], 
                  outname = sample2, path = "5.Biomaker")
}

#####-----STAMP-----#####
for (time in table(Group_info$times) %>% names() %>%.[-1]) {
  # time = '2019.6'
  print(time)
  sample <- paste0(paste0(rep(time, 12), rep(c("_Control_", "_Warming_"), each=6)), 1:6)
  sample2 <- paste0(time, "_Warming_vs_Control")
  STAMP_function(data = "seqtab",
                 sample1 = sample[1], sample2 = sample[2], sample3 = sample[3], sample4 = sample[4], sample5 = sample[5], sample6 = sample[6],
                 sample7 = sample[7], sample8 = sample[8], sample9 = sample[9], sample10 = sample[10], sample11 = sample[11],sample12 = sample[12], 
                 groupname  = sample2, path = "5.Biomaker")
}


#####-----Summary-----#####
# 移出两个不要的月份的数据

# list.files(path = "./5.biomaker/", pattern = "2019.7|2020.11", full.names = T) %>% file.remove()

# edgeR 
edgeR_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "EdgeR.DEG.csv", full.names = T),
                   function(x){read_delim(x)} %>% dplyr::filter(change != "Normal") %>% dplyr::pull(ASV)) 
Reduce(union, edgeR_DEASV) %>% length()


# DESeq2
DESeq2_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "DESeq2.DEG.csv", full.names = T),
                    function(x){read_delim(x)} %>% dplyr::filter(change != "Normal") %>% dplyr::pull(ASV)) 
Reduce(union, DESeq2_DEASV) %>% length() 


# STAMP
STAMP_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "STAMP_Genus_Level_stat_out.csv", full.names = T),
                   function(x){read_delim(x) %>% dplyr::filter(p.value < 0.05) %>% dplyr::pull(Variable)})
Reduce(union, STAMP_DEASV) %>% length() 


# 首先获取序列，然后构建进化树
Reduce(union, list(Reduce(union, STAMP_DEASV),
                       Reduce(union, DESeq2_DEASV),
                       Reduce(union, edgeR_DEASV))) -> STAMP_DESeq2_EdgeR.union.ID

length(STAMP_DESeq2_EdgeR.union.ID) # 5810
