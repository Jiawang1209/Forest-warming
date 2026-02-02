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
library(rstatix)
library(GGally)
library(lme4)
library(grid)
library(emmeans)
library(lmerTest)
library(ggnewscale)
library(ggh4x)

dir.create("13.functional_gene")

####----load Data----####

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

Sample_id <- read_xlsx(path = "Sample_ID_File.xlsx", col_names = T) %>%
  dplyr::mutate(Sample2 = str_c(Year, Date, Treatmen, sep = "_")) %>%
  dplyr::select(Sample, Sample2)

df_clean <- read_xlsx(path = "Relative_Expression_Functional_Gene_clean.xlsx",
                      sheet = 2)

functional_gene_estimate_add_group <- read_xlsx(path = "functional_gene_estimate_add_group_20251104.xlsx", sheet = 1)
functional_gene_estimate_add_group2 <- read_xlsx(path = "functional_gene_estimate_add_group_20251104.xlsx", sheet = 2)

# 这里是先做了标准化
df_clean_group_scale <- df_clean %>%
  dplyr::mutate(across(`iso-plu`:`dsrB`, ~ scale2(.x, na.rm = TRUE))) %>%
  tidyr::pivot_longer(cols = -c(Sample, Sample2, Time, Treatment),
                      names_to = "Gene",
                      values_to = "relative_abudance") %>%
  dplyr::left_join(functional_gene_estimate_add_group2, by = c("Gene" = "group")) %>%
  dplyr::select(-Group) %>%
  dplyr::group_by(Sample, Sample2, Time, Treatment, Functional) %>%
  dplyr::summarise(mean = mean(relative_abudance)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = Functional,
                     values_from = mean) %>%
  dplyr::arrange(Time, Treatment) %>%
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
  dplyr::mutate(across(Time, ~ scale2(.x, na.rm = T)))  %>%
  dplyr::mutate(block = rep(c(2,2,3,3,6,6, 1,1,4,4,5,5), times = 11)) %>%
  dplyr::mutate(block = factor(block)) %>%
  dplyr::mutate(Treatment = case_when(
    Treatment == "control" ~ "Control",
    Treatment == "warming" ~ "Warming"
  )) %>%
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))

df_clean_group_scale


emm_all_list <- list()
emm_time_list <- list()

gene_id <- colnames(df_clean_group_scale)[5:24]

for (g in gene_id) {
  # g = "Ammonification"
  print(g)
  
  data_tmp <- df_clean_group_scale %>% 
    dplyr::select(all_of(c(g, "Time", "Treatment","block"))) %>%
    dplyr::mutate(Value = .data[[g]])
  data_tmp
  
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
  dplyr::left_join(functional_gene_estimate_add_group2 %>%
                     dplyr::select(-1) %>%
                     dplyr::distinct(Functional, .keep_all = T),
                   by = c("group" = "Functional")
  ) %>%
  dplyr::mutate(change = case_when(
    p.value < 0.05 & estimate > 0 ~ "Positive",
    p.value < 0.05 & estimate < 0 ~ "Negative",
    .default = "Nochange"
  )) 


emm_time_df <- do.call(rbind, emm_time_list) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) %>%
  dplyr::left_join(functional_gene_estimate_add_group2 %>%
                     dplyr::select(-1) %>%
                     dplyr::distinct(Functional, .keep_all = T),
                   by = c("group" = "Functional")
  ) %>%
  dplyr::mutate(change = case_when(
    p.value < 0.05 & estimate > 0 ~ "Positive",
    p.value < 0.05 & estimate < 0 ~ "Negative",
    .default = "Nochange"
  )) 

####----Plot----####
# Warming
y_order <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(Group, desc(estimate)) %>%
  dplyr::pull(group)

y_order

y_order <- c(y_order[1:3],y_order[5:7], y_order[4], 
             y_order[8:20]
             # y_order[10:11],
             # y_order[13:15],
             # y_order[12],
             # y_order[16:20]
             )


p_warming <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::mutate(group = factor(group, levels = rev(y_order), ordered = T)) %>%
  ggplot() + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = Inf, fill = "#feb24c", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 14.5, fill = "#fa9fb5", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 11.5, ymax = 13.5, fill = "#7fcdbb", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 5.5, ymax = 11.5, fill = "#74c476", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 5.5, fill = "#9e9ac8", alpha = 0.3) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2.5, fill = "#fd8d3c", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(aes(x = estimate, y = group, color = change), size = 2.5) + 
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low , xmax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = emm_all_df %>%
              dplyr::filter(term == "TreatmentWarming") %>% 
              dplyr::filter(change == "Positive"),
            mapping = aes(x = conf.high + 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) +
  geom_text(data = emm_all_df %>%
              dplyr::filter(term == "TreatmentWarming") %>% 
              dplyr::filter(change == "Negative"),
            mapping = aes(x = conf.low - 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) + 
  scale_color_manual(values = c("Positive" = "#ef6548",
                                "Negative" = "#3690c0",
                                "Nochange" = "#878787")) + 
  scale_x_continuous(
    limits = c(-2.5, 1.5),
    expand = expansion(mult = c(0.02, 0.02)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming") + 
  labs(x = "Effect size", y = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 1),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3,
                               color = rev(c(rep("#fd8d3c", 6),
                                             rep("#f768a1", 1),
                                             rep("#41b6c4", 2),
                                             rep("#41ab5d", 6),
                                             rep("#807dba", 3),
                                             rep("#fc4e2a", 2))),
                               size = 15),
    axis.title = element_text(color = "#000000", size = 15),
    # plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 2.5
  )  + 
  coord_cartesian(clip = "off") +
  coord_cartesian(clip = "off") + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#fd8d3c", color = "#fd8d3c", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 0.5,
                    ymax = 2.5) + 
  annotation_custom(grob = textGrob(label = "S cycling", gp = gpar(col = "#fc4e2a", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 1.5,
                    ymax = 1.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#9e9ac8", color = "#9e9ac8", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 2.5,
                    ymax = 5.5) + 
  annotation_custom(grob = textGrob(label = "P cycling", gp = gpar(col = "#807dba", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 4,
                    ymax = 4) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#74c476", color = "#74c476", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 5.5,
                    ymax = 11.5) + 
  annotation_custom(grob = textGrob(label = "N cycling", gp = gpar(col = "#41ab5d", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 8.5,
                    ymax = 8.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#7fcdbb", color = "#7fcdbb", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 11.5,
                    ymax = 13.5) + 
  annotation_custom(grob = textGrob(label = "Methane\nmetabolism", gp = gpar(col = "#41b6c4", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 12.5,
                    ymax = 12.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#fa9fb5", color = "#fa9fb5", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 13.5,
                    ymax = 14.5) + 
  # annotation_custom(grob = textGrob(label = "C fixation", gp = gpar(col = "#f768a1", cex = 1.5), rot = 90),
  #                   xmin = -6,
  #                   xmax = -6,
  #                   ymin = 35,
  #                   ymax = 35) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#feb24c", color = "#feb24c", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 14.5,
                    ymax = 20.65) + 
  annotation_custom(grob = textGrob(label = "C cycling", gp = gpar(col = "#fd8d3c", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 17,
                    ymax = 17) 
  # annotation_custom(grob = textGrob(label = "C degradation", gp = gpar(col = "#fd8d3c", cex = 1.25), rot = 90),
  #                   xmin = -4.75,
  #                   xmax = -4.75,
  #                   ymin = 18,
  #                   ymax = 18)
  
p_warming

# Warming x Time
p_warmingTime <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::mutate(group = factor(group, levels = rev(y_order), ordered = T)) %>%
  ggplot() + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = Inf, fill = "#feb24c", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 14.5, fill = "#fa9fb5", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 11.5, ymax = 13.5, fill = "#7fcdbb", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 5.5, ymax = 11.5, fill = "#74c476", alpha = 0.3) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 5.5, fill = "#9e9ac8", alpha = 0.3) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2.5, fill = "#fd8d3c", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(aes(x = estimate, y = group, color = change), size = 2.5) + 
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low , xmax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = emm_all_df %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>% 
              dplyr::filter(change == "Positive"),
            mapping = aes(x = conf.high + 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) +
  geom_text(data = emm_all_df %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>% 
              dplyr::filter(change == "Negative"),
            mapping = aes(x = conf.low - 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) + 
  scale_color_manual(values = c("Positive" = "#ef6548",
                                "Negative" = "#3690c0",
                                "Nochange" = "#878787")) + 
  scale_x_continuous(
    limits = c(-2.5, 1.5),
    expand = expansion(mult = c(0.02, 0.02)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming x Time") + 
  labs(x = "Effect size", y = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 1),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3,
                               color = rev(c(rep("#fd8d3c", 6),
                                             rep("#f768a1", 1),
                                             rep("#41b6c4", 2),
                                             rep("#41ab5d", 6),
                                             rep("#807dba", 3),
                                             rep("#fc4e2a", 2))),
                               size = 15),
    axis.title = element_text(color = "#000000", size = 15),
    # plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 2.5
  )  + 
  coord_cartesian(clip = "off") +
  coord_cartesian(clip = "off") + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#fd8d3c", color = "#fd8d3c", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 0.5,
                    ymax = 2.5) + 
  annotation_custom(grob = textGrob(label = "S cycling", gp = gpar(col = "#fc4e2a", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 1.5,
                    ymax = 1.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#9e9ac8", color = "#9e9ac8", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 2.5,
                    ymax = 5.5) + 
  annotation_custom(grob = textGrob(label = "P cycling", gp = gpar(col = "#807dba", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 4,
                    ymax = 4) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#74c476", color = "#74c476", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 5.5,
                    ymax = 11.5) + 
  annotation_custom(grob = textGrob(label = "N cycling", gp = gpar(col = "#41ab5d", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 8.5,
                    ymax = 8.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#7fcdbb", color = "#7fcdbb", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 11.5,
                    ymax = 13.5) + 
  annotation_custom(grob = textGrob(label = "Methane\nmetabolism", gp = gpar(col = "#41b6c4", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 12.5,
                    ymax = 12.5) + 
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#fa9fb5", color = "#fa9fb5", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 13.5,
                    ymax = 14.5) + 
  # annotation_custom(grob = textGrob(label = "C fixation", gp = gpar(col = "#f768a1", cex = 1.5), rot = 90),
  #                   xmin = -6,
  #                   xmax = -6,
  #                   ymin = 35,
  #                   ymax = 35) +
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#feb24c", color = "#feb24c", lwd = 0, alpha = 0.2)),
                    xmin = -2.6,
                    xmax = -6.5,
                    ymin = 14.5,
                    ymax = 20.65) + 
  annotation_custom(grob = textGrob(label = "C cycling", gp = gpar(col = "#fd8d3c", cex = 1.5), rot = 90),
                    xmin = -6,
                    xmax = -6,
                    ymin = 17,
                    ymax = 17) 
  # annotation_custom(grob = textGrob(label = "C degradation", gp = gpar(col = "#fd8d3c", cex = 1.25), rot = 90),
  #                   xmin = -4.75,
  #                   xmax = -4.75,
  #                   ymin = 18,
  #                   ymax = 18)

# p_warmingTime


# every time
# everytime
emm_time_df_2 <- emm_time_df %>%
  dplyr::mutate(Time = rep(c("2017/08",
                             "2019/07", "2019/09", 
                             "2020/07", "2020/09",
                             "2021/06", "2021/09",
                             "2022/08", "2022/10",
                             "2023/07", "2023/09"),  times = 20)) %>%
  dplyr::mutate(Time = factor(Time, levels = c("2017/08",
                                               "2019/07", "2019/09", 
                                               "2020/07", "2020/09",
                                               "2021/06", "2021/09",
                                               "2022/08", "2022/10",
                                               "2023/07", "2023/09"), ordered = T)) %>%
  dplyr::select(Time, group, estimate, conf.low, conf.high, signif, p.value, Group) %>%
  dplyr::mutate(group = factor(group, levels = rev(y_order), ordered = T)) %>%
  dplyr::mutate(change = case_when(
    p.value < 0.05 & estimate > 0 ~ "Positive",
    p.value < 0.05 & estimate < 0 ~ "Negative",
    .default = "Nochange"
  )) %>%
  dplyr::mutate(Year = str_split(Time, pattern = "/", simplify = T)[,1])

# emm_time_df_2


p_emm_time_df_2 <- ggplot(data = emm_time_df_2, mapping = aes(x = interaction(Time, Year), y = group)) + 
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(Group == c("S cycling")),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#fff7bc", high = "#fc8d59",
                      guide = guide_legend(order = 6),
                      name = "S cycling") + 
  new_scale_fill() +
  # P Cycling  
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(Group == "P cycling"),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#dadaeb", high = "#8c96c6",
                      guide = guide_legend(order = 5),
                      name = "P cycling") + 
  new_scale_fill() +
  # N Cycle
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(Group == "N cycling"),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#d9f0a3", high = "#78c679",
                      guide = guide_legend(order = 4),
                      name = "N cycling") + 
  new_scale_fill() +
  # Methane metabolism
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(Group == "Methane metabolism"),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#e0f3db", high = "#7bccc4",
                      guide = guide_legend(order = 3),
                      name = "Methane metabolism") + 
  new_scale_fill() +
  # C fixation
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(group == "C fixation"),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#fde0dd", high = "#f768a1",
                      guide = guide_legend(order = 2),
                      name = "C fixation") + 
  new_scale_fill() +
  # C degradation
  geom_tile(data = emm_time_df_2 %>% dplyr::filter(group %in% c("Pectin", "Cellulose", "Chitin", 
                                                                "Lignin", "Hemicellulose", "Starch")),
            aes(x = interaction(Time, Year, sep = ":"), y = group, fill = estimate), color = "#969696") + 
  scale_fill_gradient(low = "#fff7bc", high = "#ec7014",
                      guide = guide_legend(order = 1),
                      name = "C degradation") + 
  geom_text(aes(x = interaction(Time, Year, sep = ":"), y = group, label = signif), size = 5.5) + 
  scale_x_discrete(guide = guide_axis_nested(delim = ":")) +
  guides(x.sec = guide_axis_manual(breaks = 1:11,
                                   labels = c("17.8.29", 
                                              "19.7.22", "19.9.22",
                                              "20.7.17", "20.9.21",
                                              "21.6.24", "21.9.19",
                                              "22.8.22", "22.10.1",
                                              "23.7.7", "23.9.23"))) + 
  labs(x = "Date", y = "") + 
  # coord_fixed(ratio = 0.35) +
  theme_bw()  + 
  theme(
    axis.title.x=element_text(colour = 'black',size=15),
    axis.title.y = element_text(colour = 'black',size=15),
    axis.text = element_text(colour = "black", size = 11),
    axis.text.x = element_text(colour = "black", size = 15, angle = 45, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    # legend.background = element_roundrect(color = "#808080",linetype = 1),
    # legend.position = "top",
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10),
    # panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2, color = "#000000"),
    aspect.ratio = 2.5
    # panel.background = element_blank()
  ) 

# p_emm_time_df_2


p_combine <- p_warming + 
  plot_spacer() + 
  p_warmingTime + 
  # plot_spacer() + 
  p_emm_time_df_2 + 
  plot_layout(widths = c(1, 0.001, 1, 1))


p_combine


ggsave(filename = "13.functional_gene/New_functional_gene_estimate_add_group_20251104.pdf",
       plot = p_combine,
       height = 13,
       width = 30)


####----20251214----####
emm_all_df_C_group <- emm_all_df %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::filter(Group == "C cycling") %>%
  dplyr::filter(group != "C fixation") %>%
  dplyr::arrange(estimate) %>%
  dplyr::mutate(group = factor(group, levels = group, ordered = T))

emm_all_df_C_group %>%
  ggplot(aes(x = estimate, y = group)) + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = Inf, ymax = -Inf, fill = "#feb24c", alpha = 0.3) + 
  geom_bar(aes(x = estimate, y = group, fill = change), stat = "identity", width = 0.5) + 
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low , xmax = conf.high, color = change), width = 0.2, linewidth = 1) +
  geom_text(data = emm_all_df_C_group %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>% 
              dplyr::filter(change == "Positive"),
            mapping = aes(x = conf.high + 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) +
  geom_text(data = emm_all_df_C_group %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>% 
              dplyr::filter(change == "Negative"),
            mapping = aes(x = conf.low - 0.3, y = group, label = signif, color = change), 
            size = 7.5,
            vjust = 0.75,
            show.legend = F) + 
  scale_fill_manual(values = c("Positive" = "#ef6548",
                                "Negative" = "#3690c0",
                                "Nochange" = "#878787")) + 
  scale_color_manual(values = c("Positive" = "#ef6548",
                               "Negative" = "#3690c0",
                               "Nochange" = "#878787")) + 
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  labs(x = "Effect size", y = "") + 
  coord_flip() + 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 1),
    # panel.background = element_rect(fill = "#fd8d3c"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(color = "#000000", size = 15),
    axis.title = element_text(color = "#000000", size = 15),
    # plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 0.75
  ) -> p_emm_all_df_C_group_bar

ggsave(filename = "13.functional_gene/emm_all_df_C_group_bar.pdf",
       plot = p_emm_all_df_C_group_bar,
       height = 6,
       width = 8
       )
