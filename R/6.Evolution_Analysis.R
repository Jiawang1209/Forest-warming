####----Project:5.Evolution Analysis----####
####----Author:LiuYueR----####
####----Date:20240131----####
####----Email:yueliu1115@163.com----####


rm(list = ls())

####----get Evolution tree----####

# perl ~/script/get_fa_by_id.pl STAMP_union_id_308.id.txt otus.fa > STAMP_union_id_308.id.fa
# muscle -in STAMP_union_id_308.id.fa -out STAMP_union_id_308.id.muscle
# iqtree -s STAMP_union_id_308.id.muscle
# 拿到结果 STAMP_union_id_308.id.muscle.treefile

####----load R Package ----####
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
set.seed(1115)

dir.create("5.Evolution_Analysis")

####----load Data ----####
####----O Layer ----####
# tree file
# 换成5810的树才可以
# edgeR 
edgeR_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "EdgeR.DEG.csv", full.names = T),
                          function(x){
                            read_delim(x) %>% 
                              dplyr::filter(change != "Normal") %>% 
                              dplyr::pull(ASV)
                          }
) 
# Reduce(union, edgeR_DEASV) %>% length()


# DESeq2
DESeq2_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "DESeq2.DEG.csv", full.names = T),
                           function(x){
                             read_delim(x) %>% 
                               dplyr::filter(change != "Normal") %>% 
                               dplyr::pull(ASV)
                           }
) 
# Reduce(union, DESeq2_DEASV) %>% length() 


# STAMP
STAMP_DEASV <- purrr::map(list.files(path = "5.Biomaker/", pattern = "STAMP_Genus_Level_stat_out.csv", full.names = T),
                          function(x){
                            read_delim(x) %>% 
                              dplyr::filter(p.value < 0.05) %>% 
                              dplyr::pull(Variable)
                          }
)
# Reduce(union, STAMP_DEASV) %>% length() 

Reduce(union, list(Reduce(union, STAMP_DEASV),
                   Reduce(union, DESeq2_DEASV),
                   Reduce(union, edgeR_DEASV))) -> STAMP_DESeq2_EdgeR.union.ID

# 目前含有5810个ASV的树，然后我们对这个5810的结果构建进化树

# 提取序列文件
ASV_seq <- Biostrings::readDNAStringSet(filepath = "Input_File/otus.fa") %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ASV") %>%
  purrr::set_names(c("ASV", "Seq"))

ASV_seq_sub <- ASV_seq %>% 
  dplyr::filter(ASV %in% STAMP_DESeq2_EdgeR.union.ID)

# ASV_seq_sub %>%
#   dplyr::mutate(out = str_c(">", ASV, "\n", Seq)) %>%
#   dplyr::select(out) %>%
#   write.table(file = "5.Evolution_Analysis/5810.ASV.seq.fa",
#               quote = F,
#               row.names = F,
#               col.names = F)

# 然后用于用这个结果去构建进化树

# nohup muscle -align 5810.ASV.seq.fa -output 5810.ASV.seq.muscle &
# nohup  mafft --auto  5810.ASV.seq.fa >  5810.ASV.seq.mafft 2> 5810.ASV.seq.mafft.log &

# nohup fasttree 5810.ASV.seq.mafft > 5810.ASV.seq.mafft.fasttree 2> 5810.ASV.seq.mafft.fasttree.log &
# fasttree Biomaker.afa > Biomaker.afa.fasttree

# 批量读取差异物种


tree <- ggtree::read.tree(file = "Input_File/5810.ASV.seq.mafft.fasttree")
tree2 <- treeio::read.newick(file = 'Input_File/5810.ASV.seq.mafft.fasttree',
                             node.label = 'support')

as_tibble(tree2)
as_tibble(tree)



# DEASV_union_id
STAMP_DESeq2_EdgeR.union.ID

# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean.txt", header = T, row.names = 1, 
                     sep = "\t", check.names = F,  comment.char = "")
seqtab <- seqtab[which(rowSums(seqtab) != 0),] %>%
  dplyr::filter(rownames(.) %in% STAMP_DESeq2_EdgeR.union.ID)
dim(seqtab)

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1) %>%
  dplyr::filter(rownames(.) %in% STAMP_DESeq2_EdgeR.union.ID)
dim(taxa)

taxa2 <- taxa %>%
  dplyr::mutate(Phylum = ifelse(Phylum == "Thaumarchaeota", "Unassigned", Phylum))

taxa2$Phylum %>% table() %>% sort(., decreasing = T)

# rare OTU table
set.seed(1115)
otu_rare = as.data.frame(t(rrarefy(t(seqtab), colSums(seqtab) %>% min())))
otu_rare_relative <- otu_rare / colSums(otu_rare)

colSums(otu_rare)
colSums(seqtab)

# otu_rare_relative <- otu_rare / colSums(otu_rare)
# colSums(otu_rare_relative)
# dim(otu_rare_relative)

taxa_stamp <- taxa[rownames(taxa) %in% STAMP_DESeq2_EdgeR.union.ID,]


tax_df <- taxa_stamp %>%
  rownames_to_column(var = "label")

tree2 <- groupClade(tree, c(6043,
                            11083,
                            11102,
                            6044,
                            6045,
                            6047,
                            8143,
                            8283,
                            7582,
                            8069,
                            7463,
                            8354,
                            10459,
                            10574,
                            8355,
                            8716,
                            8357,
                            8772,
                            8773,
                            8886,
                            9448
                            ))


# 将进化树的结果和分类结果合并在一起

tree_label_tax_info <- as_tibble(tree2) %>% 
  dplyr::left_join(tax_df, by = c("label" = "label")) 


ggtree(tree2, branch.length = "none", mapping = aes(color = group), size = 2) %<+% tax_df +
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  new_scale_color() + 
  geom_tiplab(aes(color = Phylum), offset = 0.5) +
  scale_color_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                               "Proteobacteria" = "#f16913", 
                               "Verrucomicrobia" = "#ff7f00",
                               "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                               "Firmicutes" = "#ffed6f",
                               "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                               "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                               "candidate_division_WPS-1" = "#fa9fb5", "candidate_division_WPS-1_Uassigned" = "#fa9fb5",
                               "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                               "Actinobacteria" = "#980043")) +
  new_scale_color() +
  # geom_tree(aes(color = Phylum)) +
  geom_nodelab(aes(label = node)) -> tmp.p2

# reroot
ggtree::rotate(tmp.p2, node = 6043) -> tmp.p2.rotate



####----root一下看看效果----####
tree_reroot <- root(tree2 %>% as.phylo(), node=6043)

tree_reroot2 <- groupClade(tree_reroot, c(6043,
                            11083,
                            11102,
                            6044,
                            6045,
                            6047,
                            8143,
                            8283,
                            7582,
                            8069,
                            7463,
                            8354,
                            10459,
                            10574,
                            8355,
                            8716,
                            8357,
                            8772,
                            8773,
                            8886,
                            9448
))

ggtree(tree_reroot2, branch.length = "none", 
       mapping = aes(color = group),
       size = 2) %<+% tax_df +
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  new_scale_color() + 
  geom_tiplab(aes(color = Phylum), offset = 0.5) +
  scale_color_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#fa9fb5", "candidate_division_WPS-1_Uassigned" = "#fa9fb5",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) +
  new_scale_color() +
  # geom_tree(aes(color = Phylum)) +
  geom_nodelab(aes(label = node)) -> tmp.p2_reroot


ggsave(filename = "5.Evolution_Analysis/20250311.tmptree2.reroot.pdf",
       plot = tmp.p2_reroot,
       height = 800,
       width = 100,
       limitsize = F)



# 然后从进化树的结果中
get_taxa_name(tmp.p2) %>%
  as.data.frame() %>%
  purrr::set_names("ASV") %>%
  tidyr::fill(Phylum) %>%
  dplyr::mutate(Phylum1 = str_remove(string = Phylum, pattern = "_Uassigned"),
                Phylum2 = str_remove(string = Phylum, pattern = ".*_")) -> Tree_phylum_annotation_from_tree
  

Tree_phylum_annotation_from_tree


as_tibble(tree2) %>% 
  dplyr::left_join(Tree_phylum_annotation_from_tree, by = c("label" = "ASV")) %>%
  dplyr::filter(str_detect(label, pattern = "ASV")) %>%
  dplyr::mutate(group = as.character(group)) -> tax_df2


####----20251020 绘制进化树 + 线性混合模型----####

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

seatab_sub <- otu_rare %>%
  dplyr::filter(rownames(.) %in% as_tibble(tree2)[["label"]]) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(across(where(is.numeric), ~ scale2(.x, na.rm = TRUE))) %>%
  dplyr::mutate(Time = c(
    rep(0, times = 12),
    rep(0.87, times = 12),
    rep(1.04, times = 12),
    rep(1.85, times = 12),
    rep(2.04, times = 12),
    rep(2.79, times = 12),
    rep(3.06, times = 12),
    rep(3.93, times = 12),
    rep(4.06, times = 12),
    rep(4.83, times = 12),
    rep(5.04, times = 12)
  )) %>%
  dplyr::mutate(Treatment = rep(c(rep("Control",6),
                                  rep("Warming", 6)),
                                times = 11)) %>%
  dplyr::mutate(across(Time, ~ scale2(.x, na.rm = T)))  %>%
  dplyr::mutate(block = rep(c(2,2,3,3,6,6, 1,1,4,4,5,5), times = 11)) %>% # 这是不去处2017的时候
  dplyr::mutate(block = factor(block)) %>%
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))


# 线性混合模型
LMM_result <- list()
LMM_everytime_result <- list()
ASV_id <- colnames(seatab_sub)[1:5810]

for (g in ASV_id) {
  # g = "ASV_1006"
  print(g)
  
  data_tmp <- seatab_sub %>% 
    dplyr::select(all_of(c(g, "Time", "Treatment","block"))) %>%
    dplyr::mutate(Value = .data[[g]])
  
  # LMM
  fm1 <- lmerTest::lmer(Value~Treatment * Time +(1|block),data=data_tmp)
  
  out1 <- broom.mixed::tidy(fm1, effects = "fixed", conf.int = TRUE) %>%
    dplyr::mutate(group = g)
  
  emm_time <- emmeans(fm1, pairwise ~ Treatment | Time, 
                      at = list(Time = c(-1.67325487237833, -1.13062359048539, -1.02459219057527,
                                         -0.519383755709437, -0.400878073456956, 0.0669075143817825,
                                         0.235310326003728, 0.777941607896665, 0.859024443122047,
                                         1.33928431330315, 1.470264277898))
  )
  
  ctr <- contrast(emm_time, method = "revpairwise")
  
  ctr %>% broom.mixed::tidy()
  
  ctr_tc <- contrast(emm_time, method = "trt.vs.ctrl", ref = "Control")
  
  out2 <- ctr_tc %>% broom.mixed::tidy(conf.int = TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(group = g)
  
  LMM_result[[g]] <- out1
  LMM_everytime_result[[g]] <- out2
  
}

LMM_result_df <- do.call(rbind, LMM_result) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) %>%
  dplyr::filter(term %in% c("TreatmentWarming", "TreatmentWarming:Time")) %>%
  dplyr::mutate(
    Change = case_when(
      estimate > 0 & p.value < 0.05 ~ "Up",
      estimate < 0 & p.value < 0.05  ~ "Down",
      .default = "Normal"
    )
  )

dim(LMM_result_df)

write_xlsx(LMM_result_df,
           path = "5.Evolution_Analysis/LMM_out_5810_20251103.xlsx")

LMM_result_df$Change %>% table()


LMM_everytime_result_df <- do.call(rbind, LMM_everytime_result) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) 

# 然后统计一下有多少上调和多少下调，然后绘制扇形图

# 这个是 warming x time的结果
LMM_result_df %>%
  dplyr::select(group, term, Change, estimate) %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::left_join(Tree_phylum_annotation_from_tree, by = c("group" = "ASV")) %>%
  dplyr::rename(ASV = group) -> LMM_result_df_LMM_out

# 这个是 Warmig的结果
  LMM_result_df %>%
    dplyr::select(group, term, Change, estimate) %>%
    dplyr::filter(term == "TreatmentWarming") %>%
    dplyr::left_join(Tree_phylum_annotation_from_tree, by = c("group" = "ASV")) %>%
    dplyr::rename(ASV = group)  -> LMM_result_df_LMM_out2


scale3 <- function(x, na.rm = FALSE) {x/( max(x) - min(x))}

# warming x time phylum stat
table(LMM_result_df_LMM_out$Change,
      LMM_result_df_LMM_out$Phylum1) %>%
  as.data.frame() %>%
  purrr::set_names(c("Change", "Phylum", "Number")) %>%
  tidyr::pivot_wider(names_from = Change, values_from = Number) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Total = sum(c(Down, Up))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(x = rep(seq(1,18, 2), times = 2)[-length(rep(seq(1,18, 2), times = 2))],
                y = rep(c(1, 3), each = 9)[-length(rep(c(1, 3), each = 8))],
                group = 1:17,
                radius = scale3(Total, na.rm = T) 
  )-> LMM_stat

LMM_stat

ggplot() + 
  geom_scatterpie(data = LMM_stat, 
                  mapping = aes(x=x, y = y, group = group
                                # , r = radius*1.25
                  ),
                  cols = colnames(LMM_stat)[c(2,4)],
                  alpha = 0.85) + 
  # geom_scatterpie_legend(LMM_stat$radius*1.25, x=-5, y=2) + 
  geom_text(data = LMM_stat,
            mapping = aes(x = x, y = y + 1, label = Phylum)) + 
  scale_fill_manual(values = c("#6577bd", "#f7ec16")) +
  coord_equal() + 
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) -> LMM_Phylum_stat

LMM_Phylum_stat


ggsave(filename = "5.Evolution_Analysis/LMM_Phylum_stat_warmingxtime.pdf",
       plot = LMM_Phylum_stat,
       height = 10,
       width = 35)


# warming phylum stat

table(LMM_result_df_LMM_out2$Change,
      LMM_result_df_LMM_out2$Phylum1) %>%
  as.data.frame() %>%
  purrr::set_names(c("Change", "Phylum", "Number")) %>%
  tidyr::pivot_wider(names_from = Change, values_from = Number) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Total = sum(c(Down, Up))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(x = rep(seq(1,18, 2), times = 2)[-length(rep(seq(1,18, 2), times = 2))],
                y = rep(c(1, 3), each = 9)[-length(rep(c(1, 3), each = 8))],
                group = 1:17,
                radius = scale3(Total, na.rm = T) #,
                # radius2 = case_when(
                #   radius < 0.3  & 0.1 < radius ~ radius * 5,
                #   radius < 0.1  ~ radius * 20,
                #   .default = as.numeric(radius)),
                # transformed_scaled = scales::rescale(radius*5, to = c(0, 2)),
                # log_trans = log10(radius + 1) * 10
  )-> LMM_stat2

LMM_stat2

ggplot() + 
  geom_scatterpie(data = LMM_stat2, 
                  mapping = aes(x=x, y = y, group = group
                                # , r = radius*1.25
                  ),
                  cols = colnames(LMM_stat)[c(2,4)],
                  alpha = 0.85) + 
  # geom_scatterpie_legend(LMM_stat$radius*1.25, x=-5, y=2) + 
  geom_text(data = LMM_stat2,
            mapping = aes(x = x, y = y + 1, label = Phylum)) + 
  scale_fill_manual(values = c("#6577bd", "#f7ec16")) +
  coord_equal() + 
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) -> LMM_Phylum_stat2


ggsave(filename = "5.Evolution_Analysis/LMM_Phylum_stat_warming.pdf",
       plot = LMM_Phylum_stat2,
       height = 10,
       width = 35)




background_df <- data.frame(
  ASV = LMM_result_df$group,
  Up = rep(0.7, 5810),
  Down = rep(-0.7, 5810)
) %>%
  tidyr::pivot_longer(cols = -ASV, names_to = "Type", values_to = "Value")

# 先画里面的进化树
p1 <- ggtree(tree2, layout='slanted', branch.length='none',
             aes(color = group)) +
  layout_circular() + 
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  xlim(NA, 10000)

p1

# 然后依次画外面的修饰的图
p1 +
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = background_df,
    geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", 
    stat="identity", 
    offset = 0.2, 
    size = 0.02, 
    alpha = 0.5,
    show.legend = F) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16")) + 
  ggnewscale::new_scale_fill() + 
  ggnewscale::new_scale_color() +
  # 先添加增温的效应
  geom_fruit(
    data = LMM_result_df %>% 
      dplyr::filter(term == "TreatmentWarming") %>% 
      dplyr::rename(label = group) %>%
      dplyr::left_join(Tree_phylum_annotation_from_tree, by = c("label" = "ASV")), 
    geom = geom_bar,
    mapping = aes(y = label, x = estimate, fill = Phylum),
    orientation="y", stat="identity", offset = -0.2, size = 0.25,
    axis.params=list(
      axis = "x", 
      hjust = 1, 
      vjust = 0,
      nbreak = 3, 
      color = "#000000",
      linetype = 2)) + 
  scale_color_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) +
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043"))+ 
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() + 
  geom_fruit(
    data = background_df, 
    geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", 
    stat="identity", 
    offset = 1.0, 
    size = 0.02, 
    alpha = 0.5,
    show.legend = F) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16"))
  geom_fruit(
    data = LMM_result_df, geom = geom_tile,
    mapping = aes(y = ASV, x = 1, fill = Change),
    height = 1.5, width = 400, offset = -0.1, size = 0.02) +
  scale_fill_manual(values = c("#6577bd", "#f7ec16")) +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df_LMM_out, geom = geom_tile,
    mapping = aes(y = ASV, x = 3, fill = Phylum1 ),
    height = 1.5,  width = 100, offset = 0.00, size = 0.02) + 
  scale_color_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043"))  -> p2
  
# p2

ggsave(filename = "./5.Evolution_Analysis/20250313.tmptree5.pdf",
       plot = p2,
       height = 10,
       width = 12)

####----reroot之后结果再画图----####
# 先画里面的进化树



p1_reroot2 <- ggtree(tree_reroot2, layout='slanted', branch.length='none',
             aes(color = group)) +
  layout_circular() + 
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  xlim(NA, 15000)

p1_reroot2

background_df1 <- data.frame(
  ASV = LMM_result_df$group,
  Up = rep(0.4, 5810),
  Down = rep(-0.4, 5810)
) %>%
  tidyr::pivot_longer(cols = -ASV, names_to = "Type", values_to = "Value")

background_df2 <- data.frame(
  ASV = LMM_result_df$group,
  Up = rep(0.5, 5810),
  Down = rep(-0.5, 5810)
) %>%
  tidyr::pivot_longer(cols = -ASV, names_to = "Type", values_to = "Value")


# 然后依次画外面的修饰的图
p1_reroot2 +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = background_df1, geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", stat="identity", offset = 0.25, size = 0.02, alpha = 0.25) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16")) + 
  ggnewscale::new_scale_fill() + 
  # 添加 warming 的结果
  geom_fruit(
    data = LMM_result_df_LMM_out2, geom = geom_bar,
    mapping = aes(y = ASV, x = estimate, fill = Phylum),
    orientation="y", stat="identity", offset = -0.2, size = 0.25
    # axis.params=list(
    #   axis = "x", hjust = 1, vjust = 0,
    #   nbreak = 3, color = "#000000", linetype = 2)
    ) + 
  # 不同的物种配色
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) +
  # 添加背景
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = background_df2, geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", stat="identity", offset = 0.26, size = 0.02, alpha = 0.25) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16")) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df_LMM_out, geom = geom_bar,
    mapping = aes(y = ASV, x = estimate, fill = Phylum),
    orientation="y", stat="identity", offset = -0.2, size = 0.25
    # axis.params=list(
    #   axis = "x", hjust = 1, vjust = 0,
    #   nbreak = 3, color = "#000000", linetype = 2)
    ) + 
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043"))  -> p2_reroot2

# p2
ggsave(filename = "5.Evolution_Analysis/20251020.reroot_LMM.pdf",
       plot = p2_reroot2,
       height = 12.5,
       width = 15)



#####-----reroot add warming effect-----#####
p1_reroot2_1 <- ggtree(tree_reroot2, layout='slanted', branch.length='none',
                     aes(color = group)) +
  layout_circular() + 
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  xlim(NA, 15000)

# 添加ID 然后进行区分
p1_reroot2_1_add_label <- ggtree(tree_reroot2, layout='fan', branch.length='none',
                       aes(color = group)) +
  # layout_circular() + 
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  geom_tiplab() +
  geom_nodelab(aes(label = node), size = 5, color = "red") + 
  xlim(NA, 80)

ggsave(filename = "./5.Evolution_Analysis/p1_reroot2_1_add_label.pdf",
       plot = p1_reroot2_1_add_label,
       height = 400,
       width = 430,
       limitsize = F)


p1_reroot2_1_add_label2 <- ggtree(tree_reroot2, layout='slanted', branch.length='none',
                                 aes(color = group)) +
  layout_circular() +
  scale_color_manual(values = c(
    # Acidobacteria / Acidobacteria_Uassigned. Spirochaetes 12 Thaumarchaeota 13 
    "0" = "#1f78b4","5" = "#1f78b4","12" = "#1f78b4", "13" = "#1f78b4", 
    # Armatimonadetes / Armatimonadetes_Uassigned -- Chloroflexi / Chloroflexi_Uassigned # Uassigned 
    "14" = "#969696",
    # Bacteroidetes Bacteroidetes_Uassigned
    "17" = "#df65b0",
    # candidate_division_WPS-1 / 	candidate_division_WPS-1_Uassigned
    "19" = "#e7298a",
    # Planctomycetes / Planctomycetes_Uassigned
    "20" = "#ce1256", 
    # Actinobacteria / Actinobacteria_Uassigned
    "21" = "#980043",
    # Candidatus_Saccharibacteria / Candidatus_Saccharibacteria_Uassigned
    "16" = "#80b1d3",
    # Chlamydiae / Chlamydiae_Uassigned Firmicutes
    "10" = "#fd8d3c", 
    # Gemmatimonadetes / Gemmatimonadetes_Uassigned
    "7" = "#ef3b2c", "8" = "#ef3b2c",
    # Latescibacteria
    "3" = "#33a02c",
    # Proteobacteria / Proteobacteria_Uassigned
    "11" = "#f16913", "4" = "#f16913", "9" = "#f16913",
    # Verrucomicrobia
    "6" = "#ff7f00",
    # 特殊的分支
    "2" = "#bdbdbd", 
    "1" = "#95d6c7","15" = "#1f78b4","18" = "#fccde5")) + 
  geom_tiplab() +
  geom_nodelab(aes(label = node), size = 5, color = "red") + 
  xlim(NA, 80)

ggsave(filename = "./5.Evolution_Analysis/p1_reroot2_1_add_label2.pdf",
       plot = p1_reroot2_1_add_label2,
       height = 400,
       width = 430,
       limitsize = F)


# 然后依次画外面的修饰的图
p1_reroot2_1 +
  # 首先是增温的效应
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = background_df, 
    geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", stat="identity", offset = 0.18, size = 0.02, alpha = 0.5) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16")) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df_LMM_out2, geom = geom_bar,
    mapping = aes(y = ASV, x = `TreatmentWarming.mean`, fill = Phylum),
    orientation="y", stat="identity", offset = -0.2, size = 0.25,
    axis.params=list(
      axix = "x", hjust = 1, vjust = 0,
      nbreak = 3, color = "#000000", linetype = 2)) + 
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) +
  new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df2, geom = geom_tile,
    mapping = aes(y = ASV, x = other, fill = signif_Treatment),
    height = 1.5,  width = 400, offset = -0.18, size = 0.02) +
  scale_fill_manual(values = c("#7fcdbb", "#fec44f")) + 
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = LMM_result_df2, geom = geom_tile,
    mapping = aes(y = ASV, x = 1, fill = Change),
    height = 1.5, width = 400, offset = -0.1, size = 0.02) +
  scale_fill_manual(values = c("#6577bd", "#f7ec16")) + 
  # 其次是增温和时间的共同的效应
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = background_df, geom = geom_bar,
    mapping = aes(y = ASV, x = Value, fill = Type),
    orientation="y", stat="identity", offset = 0.22, size = 0.02, alpha = 0.5) + 
  scale_fill_manual(values = c("#66c2a5", "#f7ec16")) + 
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df_LMM_out, geom = geom_bar,
    mapping = aes(y = ASV, x = `TreatmentWarming:Time.mean`, fill = Phylum),
    orientation="y", stat="identity", offset = -0.2, size = 0.25,
    axis.params=list(
      axix = "x", hjust = 1, vjust = 0,
      nbreak = 3, color = "#000000", linetype = 2)) + 
  # 不同的物种配色
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) +
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df, geom = geom_tile,
    mapping = aes(y = ASV, x = other, fill = signif_Treatment),
    height = 1.5,  width = 400, offset = -0.2, size = 0.02) +
  scale_fill_manual(values = c("#7fcdbb", "#fec44f")) + 
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = LMM_result_df, geom = geom_tile,
    mapping = aes(y = ASV, x = 1, fill = Change),
    height = 1.5, width = 400, offset = -0.1, size = 0.02) +
  scale_fill_manual(values = c("#6577bd", "#f7ec16")) +
  # 最后添加的是物种注释
  ggnewscale::new_scale_fill() + 
  geom_fruit(
    data = LMM_result_df_LMM_out, geom = geom_tile,
    mapping = aes(y = ASV, x = 3, fill = Phylum1 ),
    height = 1.5,  width = 100, offset = -0.01, size = 0.02) + 
  scale_fill_manual(values = c("Acidobacteria" = "#1f78b4", "Uassigned" = "#969696", "Latescibacteria" = "#33a02c",
                                "Proteobacteria" = "#f16913", 
                                "Verrucomicrobia" = "#ff7f00",
                                "Gemmatimonadetes_Uassigned" = "#ef3b2c", "Gemmatimonadetes" = "#ef3b2c", 
                                "Firmicutes" = "#ffed6f",
                                "Candidatus_Saccharibacteria_Uassigned" = "#80b1d3", "Candidatus_Saccharibacteria" = "#80b1d3",
                                "Bacteroidetes" = "#df65b0", "Bacteroidetes_Uassigned" = "#df65b0",
                                "candidate_division_WPS-1" = "#e7298a", "candidate_division_WPS-1_Uassigned" = "#e7298a",
                                "Planctomycetes" = "#ce1256",  "Planctomycetes_Uassigned" = "#ce1256",
                                "Actinobacteria" = "#980043")) -> p2_reroot2_add_warming_effects

# p2_reroot2_add_warming_effects

# p2

ggsave(filename = "./5.Evolution_Analysis/20250313.tmptree5_reroot2_add_warming_effects.pdf",
       plot = p2_reroot2_add_warming_effects,
       height = 15,
       width = 20)



table(tax_df2$Phylum,
      tax_df2$group)


get_taxa_name(p2) %>%
  as.data.frame() %>%
  purrr::set_names("ASV") %>%
  dplyr::left_join(tax_df2, by = c("ASV" = "label")) %>%
  View()
  



# 
ggtree(tree2, layout = "fan", open.angle=90) %<+% tax_df +
  geom_tiplab(aes(color = Phylum), offset = 0.5) +
  geom_tree(aes(color = Phylum)) +
  # xlim(NA,30) +
  geom_tree2()
# 
p_tree <- ggtree(tree, branch.length = "none", layout = "fan", open.angle=90) %<+% tax_df +
  geom_tiplab(aes(color = Phylum), offset = 0.5) +
  geom_tree(aes(color = Phylum)) +
  # xlim(NA,30) +
  geom_tree2()
# 



# 添加上之前的热图 探索一下
gene_name <- ggtree::get_taxa_name(p_tree)



taxa_sub <- taxa %>% 
  dplyr::filter(rownames(.) %in% Biomaker$ASV) %>%
  tibble::rownames_to_column(var = "label")


nodedf <- data.frame(
  node = c(1047, 1087, 1411, 1485)
)

p2 <- ggtree(tree2, 
             aes(color = group),
             layout="fan", 
             size=0.25, 
             linetype=1,
             open.angle = 5,
             branch.length = "none"
)  %<+% tax_df +
  # geom_nodelab(aes(label = node), size = 2) + 
  # geom_tiplab(aes(color = Phylum), offset = 0.5) + 
  geom_hilight(data = nodedf, mapping = aes(node=node),
               extendto=1, alpha=0.3, 
               fill="grey", color="grey50",
               size=0.05) + 
  geom_treescale() + 
  geom_star(aes(fill=Phylum, starshape=Phylum),
            position="identity",starstroke=0.5) + 
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970","#7FFF00"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE) + 
  scale_starshape_manual(values=c(1:15),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE) + 
  scale_color_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                              "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                              "#EE6A50","#8DEEEE", "#006400","#800000",
                              "#B0171F","#191970","#7FFF00"),
                     guide=guide_legend(keywidth = 0.5, 
                                        keyheight = 0.5, order=1,
                                        override.aes=list(starshape=15)),
                     na.translate=FALSE)

# ggtree::rotate(p2,775) 

# p2

p2 + 
  new_scale_fill() + 
  geom_fruit(data=edgeR_DEASV_file_change.df2, geom=geom_tile,
             mapping=aes(y=label, x=Group, fill=Value2),
             color = "grey50", offset = 0.01, size = 0.1, alpha = 0.85,
             pwidth = 0.4)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) + 
  scale_fill_manual(values=c("#3288bd","#66c2a5","#abdda4",
                             "#fdae61", "#f46d43","#d53e4f","#696969"),
                    guide=guide_legend(keywidth = 0.3,
                                       keyheight = 0.3, order=4)) +
  # scale_fill_gradient2(low = "#66c2a5", mid = "#ffffbf", high = "#f46d43", midpoint = 0) + 
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.spacing.y = unit(0.02, "cm"),
  ) -> p3

p3

ggsave(filename = "5.Evolution_Analysis/Biomaker_logfoldchange_evolution_2.pdf",
       plot = p3,
       height = 15,
       width = 15,
       limitsize = F)
