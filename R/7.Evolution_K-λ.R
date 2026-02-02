rm(list = ls())

####----load R Package----####
library(evobiR)
library(phytools)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(ggfun)
library(aplot)
library(ggstar)
library(ggnewscale)
library(evobiR)
library(phytools)
library(ape)
library(ggpmisc)
library(corrplot)
library(Biostrings)
library(tidyverse)
library(readxl)
library(writexl)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(lme4)
library(car)
library(ggrepel)
library(ggh4x)
library(broom)
library(sjPlot)
library(emmeans)
library(lmerTest)
library(grid)
library(vegan)
library(scatterpie)
library(readxl)
set.seed(1115)

dir.create("5.Evolution_Analysis")

####----load Data----####
tree <- ggtree::read.tree(file = "Input_File/5810.ASV.seq.mafft.fasttree")
tree<-multi2di(tree)

LMM <- read_xlsx(path = "Input_File/LMM_out_5810_20251103.xlsx")

####----warming----####
trait <- LMM %>%
  dplyr::select(group, term, estimate) %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  tibble::column_to_rownames(var = "group") %>%
  dplyr::select(-term) %>%
  purrr::set_names("Estimate")

trait_ar<-setNames(trait$Estimate,rownames(trait))

signal_ar_K <- phylosig(tree, trait_ar, method = "K", test = T, nsim = 999)
signal_ar_λ <- phylosig(tree, trait_ar, method = "lambda", test = T)

# saveRDS(signal_ar_K, file = "5.Evolution_Analysis/signal_ar_K.rds")
# saveRDS(signal_ar_λ, file = "5.Evolution_Analysis/signal_ar_λ.rds")

####----warming x time----####
trait2 <- LMM %>%
  dplyr::select(group, term, estimate) %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  tibble::column_to_rownames(var = "group") %>%
  dplyr::select(-term) %>%
  purrr::set_names("Estimate")

trait_warming<-setNames(trait2$Estimate,rownames(trait))

signal_warming_K <- phylosig(tree, trait_warming, method = "K", test = T, nsim = 999)
signal_warming_λ <- phylosig(tree, trait_warming, method = "lambda", test = T)

# saveRDS(signal_warming_K, file = "5.Evolution_Analysis/signal_warming_K.rds")
# saveRDS(signal_warming_λ, file = "5.Evolution_Analysis/signal_warming_λ.rds")
