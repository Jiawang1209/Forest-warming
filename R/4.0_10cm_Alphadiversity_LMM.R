rm(list = ls())

####----load R Package----####
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
library(MuMIn)
library(GGally)
library(lme4)
library(car)
library(ggh4x)
library(ggnewscale)
library(aplot)
library(broom)
library(sjPlot)
library(emmeans)
library(lmerTest)
library(grid)
library(writexl)

dir.create("4.0_10cm_Alphadiversity_LMM")

####----load Data----####
#####-----0-10cm-----#####
# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean_0_10cm.txt", header = T, row.names = 1, 
                     sep = "\t", check.names = F,  comment.char = "")
dim(seqtab)
seqtab <- seqtab[which(rowSums(seqtab) != 0),]
dim(seqtab)
colSums(seqtab)
colSums(seqtab) %>% min()

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1)
dim(taxa)

# rare
set.seed(1115)
otu_rare = as.data.frame(t(rrarefy(t(seqtab), 26000)))
colSums(otu_rare)

#####-----PD-----#####
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

alpha_result

row.names_vector <- row.names(alpha_result)

alpha_result$Sample <- row.names(alpha_result)

alpha_result$Group <- str_remove(alpha_result$Sample, pattern = "_\\d+")
alpha_result$Sample_month <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,1]
alpha_result$Treatment <- str_split(alpha_result$Group, pattern = "_", simplify = T)[,2]
alpha_result$Sample_year <- str_split(alpha_result$Group, pattern = "\\.", simplify = T)[,1]


write.csv(alpha_result, file = "./38.0_10cm_Alphadiversity_LMM/alpha_result_0_10cm_Layer.csv")

alpha_result2 <- alpha_result %>%
  dplyr::filter(!str_detect(string = Sample_month, pattern = "2019.7|2020.11"))

write.csv(alpha_result2, file = "./38.0_10cm_Alphadiversity_LMM/alpha_result_0_10cm_Layer_for_Analysis.csv")

####----Analysis----####
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


alpha_result2 <- read_delim(file = "./4.0_10cm_Alphadiversity_LMM/alpha_result_0_10cm_Layer_for_Analysis.csv",
                            delim = ",",
                            col_names = T) %>%
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
  dplyr::mutate(across(observed_species:PD, ~ scale2(.x, na.rm = T))) %>%
  dplyr::mutate(across(Time, ~ scale2(.x, na.rm = T)))  %>%
  dplyr::mutate(block = rep(c(2,2,3,3,6,6, 1,1,4,4,5,5), times = 11)) %>% # 这是不去处2017的时候
  dplyr::mutate(block = factor(block)) %>%
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))


alpha_result2


####----LMM new----####
emm_all_list <- list()
emm_time_list <- list()

gene_id <- colnames(alpha_result2)[2:11]

for (g in gene_id) {
  # g = "observed_species"
  print(g)
  
  # sub data
  data_tmp <- alpha_result2 %>% 
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
  
  emm_all_list[[g]] <- out1
  emm_time_list[[g]] <- out2
  
  
}

emm_all_df <- do.call(rbind, emm_all_list) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) %>%
  dplyr::filter(term %in% c("TreatmentWarming", "TreatmentWarming:Time")) %>%
  dplyr::mutate(group = case_when(
    group == "observed_species" ~ "Richness",
    .default = group
  ))


emm_time_df <- do.call(rbind, emm_time_list) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) %>%
  dplyr::mutate(group = case_when(
    group == "observed_species" ~ "Richness",
    .default = group
  ))


####----LMM new plot----####
group_order <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(estimate) %>%
  dplyr::filter(group %in% c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
  dplyr::mutate(group = case_when(
    group == "Shannon_1" ~ "Shannon",
    group == "Simpson_Gini" ~ "Simpson",
    .default = group
  )) %>%
  dplyr::pull(group)

group_order

group_order <- c("Richness", "PD", "Shannon", "Simpson")


p_warming <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::filter(group %in% c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
  dplyr::mutate(group = case_when(
    group == "Shannon_1" ~ "Shannon",
    group == "Simpson_Gini" ~ "Simpson",
    .default = group
  )) %>%
  dplyr::mutate(group = factor(group, levels = rev(group_order), ordered = T)) %>%
  ggplot() + 
  geom_bar(aes(x = estimate, y = group, fill = group), stat = "identity", width = 0.55) + 
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low, xmax = 0, color = group),
                width = 0.1, linewidth = 1, show.legend = F) + 
  geom_text(aes(x = conf.low - 0.2, y = group, label = signif, color = group),
            size = 10, vjust = 0.75, show.legend = F) + 
  scale_x_continuous(expand = expansion(mult = c(0.1, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  scale_fill_manual(values = c(
    "Richness" = "#74a9cf",
    "Shannon" = "#fb9a99",
    "Simpson" = "#ff7f00",
    "PD" = "#66c2a5"
  ),
  name = "") + 
  scale_color_manual(values = c(
    "Richness" = "#74a9cf",
    "Shannon" = "#fb9a99",
    "Simpson" = "#ff7f00",
    "PD" = "#66c2a5"
  ))  + 
  labs(x = "Effect size", y = "Alpha diversity") + 
  ggtitle(label = "Warming") + 
  # geom_hline(yintercept = 0) + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 1.5),
    panel.background = element_rect(fill='#fff7ec',colour = 'black'),
    axis.text = element_text(color = "#000000", size = 15),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,1,1,0,"cm"),
    # aspect.ratio = 1.5
  )


p_warming


p_warmingTime <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::filter(group %in% c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
  dplyr::mutate(group = case_when(
    group == "Shannon_1" ~ "Shannon",
    group == "Simpson_Gini" ~ "Simpson",
    .default = group
  )) %>%
  dplyr::mutate(group = factor(group, levels = rev(group_order), ordered = T)) %>%
  ggplot() + 
  geom_bar(aes(x = estimate, y = group, fill = group), stat = "identity", width = 0.55) + 
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low, xmax = 0, color = group),
                width = 0.1, linewidth = 1, show.legend = F) + 
  geom_text(aes(x = conf.low - 0.15, y = group, label = signif, color = group),
            size = 10, vjust = 0.75, show.legend = F) + 
  scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  scale_fill_manual(values = c(
    "Richness" = "#74a9cf",
    "Shannon" = "#fb9a99",
    "Simpson" = "#ff7f00",
    "PD" = "#66c2a5"
  ),
  name = "") + 
  scale_color_manual(values = c(
    "Richness" = "#74a9cf",
    "Shannon" = "#fb9a99",
    "Simpson" = "#ff7f00",
    "PD" = "#66c2a5"
  ))  + 
  labs(x = "Effect size", y = "") + 
  ggtitle(label = "Warming x Time") + 
  # geom_hline(yintercept = 0) + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 1.5),
    panel.background = element_rect(fill='#fff7ec',colour = 'black'),
    axis.text = element_text(color = "#000000", size = 15),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,0,1,0,"cm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
    # aspect.ratio = 1.5
  )

p_warmingTime


p_time_df <- emm_time_df %>%
  dplyr::filter(group %in%  c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
  dplyr::mutate(group = case_when(
    group == "Shannon_1" ~ "Shannon",
    group == "Simpson_Gini" ~ "Simpson",
    .default = group
  )) %>%
  dplyr::mutate(Time = rep(c("2017/08",
                             "2019/07", "2019/09", 
                             "2020/07", "2020/09",
                             "2021/06", "2021/09",
                             "2022/08", "2022/10",
                             "2023/07", "2023/09"),  times = 4)) %>%
  dplyr::mutate(Time = factor(Time, levels = c("2017/08",
                                               "2019/07", "2019/09", 
                                               "2020/07", "2020/09",
                                               "2021/06", "2021/09",
                                               "2022/08", "2022/10",
                                               "2023/07", "2023/09"), ordered = T)) %>%
  dplyr::select(Time, group, estimate, conf.low, conf.high, signif) %>%
  dplyr::mutate(group = factor(group, levels =rev(group_order), ordered = T))

p_time <- ggplot(data = p_time_df, mapping = aes(x = Time, y = group)) + 
  geom_tile(data = p_time_df %>% dplyr::filter(group == "Simpson"), mapping = aes(fill = estimate), color = "#000000") + 
  scale_fill_gradient(low = "#fff7bc", high = "#ec7014", name = "Simpson",
                      guide = guide_colorbar(order = 4),
                      limits = c(-0.8, 0.8),
                      breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
                      labels = c(-0.8, -0.4, 0, 0.4, 0.8)
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(data = p_time_df %>% dplyr::filter(group == "Shannon"), mapping = aes(fill = estimate), color = "#000000") + 
  scale_fill_gradient(low = "#fde0dd", high = "#dd3497", name = "Shannon",
                      guide = guide_colorbar(order = 3),
                      limits = c(-0.8, 0.8),
                      breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
                      labels = c(-0.8, -0.4, 0, 0.4, 0.8)
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(data = p_time_df %>% dplyr::filter(group == "PD"), mapping = aes(fill = estimate), color = "#000000") + 
  scale_fill_gradient(low = "#c7eae5", high = "#35978f", name = "PD",
                      guide = guide_colorbar(order = 2),
                      limits = c(-0.8, 0.8),
                      breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
                      labels = c(-0.8, -0.4, 0, 0.4, 0.8)
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(data = p_time_df %>% dplyr::filter(group == "Richness"), mapping = aes(fill = estimate), color = "#000000") + 
  scale_fill_gradient(low = "#ece7f2", high = "#2b8cbe", name = "Richness",
                      guide = guide_colorbar(order = 1),
                      limits = c(-0.8, 0.8),
                      breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
                      labels = c(-0.8, -0.4, 0, 0.4, 0.8)
  ) +
  geom_text(data = p_time_df, mapping = aes(label = signif), size = 10, vjust = 0.75) + 
  labs(x = "Date", y = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 1.5),
    panel.background = element_rect(fill='#fff7ec',colour = 'black'),
    axis.text = element_text(color = "#000000", size = 15),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(color = "#000000", size = 17.5),
    legend.position = "top",
    legend.key = element_rect(color = "#000000"),
    legend.ticks = element_line(color = "#000000"),
    legend.frame = element_rect(color = "#000000"),
    legend.text = element_text(color = "#000000", size = 11),
    legend.title = element_text(color = "#000000", size = 15),
    plot.margin = margin(1,1,1,0,"cm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


p_time


p_combine <- p_warming +
  plot_spacer() + 
  p_warmingTime + 
  plot_spacer() + 
  p_time + 
  plot_layout(# guides = "collect",
    widths = c(1, 0.1, 1, 0.1, 3))

p_combine


ggsave(filename = "4.0_10cm_Alphadiversity_LMM/Figure1_LMM_0_10cm_Alpha_diversity.pdf",
       plot = p_combine,
       height = 5.5,
       width = 20)


####----save----####
p_time_df

write_xlsx(p_time_df, path = "4.0_10cm_Alphadiversity_LMM/LMM_timepoint_out.xlsx")

p_Warming_WT <- dplyr::bind_rows(
  emm_all_df %>%
    dplyr::filter(term == "TreatmentWarming") %>%
    dplyr::filter(group %in% c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
    dplyr::mutate(group = case_when(
      group == "Shannon_1" ~ "Shannon",
      group == "Simpson_Gini" ~ "Simpson",
      .default = group
    )),
  emm_all_df %>%
    dplyr::filter(term == "TreatmentWarming:Time") %>%
    dplyr::filter(group %in% c("Richness", "PD", "Shannon_1", "Simpson_Gini")) %>%
    dplyr::mutate(group = case_when(
      group == "Shannon_1" ~ "Shannon",
      group == "Simpson_Gini" ~ "Simpson",
      .default = group
    ))
)

write_xlsx(p_Warming_WT, path = "4.0_10cm_Alphadiversity_LMM/LMM_Warming_WT_out.xlsx")
