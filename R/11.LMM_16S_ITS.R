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

dir.create("11.LMM_16S_ITS")

####----16S Phylum----####
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
Phylum_16S <- read_delim(file = "Tax_OTU_Abundance/O_Layer_16S_Phylum_abundance.csv", 
                         col_names = T, delim = ",") %>%
  dplyr::mutate(across(where(is.numeric), ~ scale2(.x, na.rm = TRUE)))  %>%
  # dplyr::slice(-c(1:12)) %>%
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
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))

Phylum_16S_relative <- read_delim(file = "Tax_OTU_Abundance/O_Layer_16S_Phylum_relative_abundance.csv", 
                         col_names = T, delim = ",") %>%
  dplyr::mutate(across(where(is.numeric), ~ scale2(.x, na.rm = TRUE)))  %>%
  # dplyr::slice(-c(1:12)) %>%
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
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))

####----LMM New 20251017----####
emm_all_list <- list()
emm_time_list <- list()


Phylum_16S_id <- colnames(Phylum_16S)[6:26]

for (g in Phylum_16S_id) {
  print(g)
  
  data_tmp <- Phylum_16S %>% 
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
  dplyr::filter(!group %in% c("Unassigned","Thaumarchaeota","Chlamydiae", 
                              "BRC1", "candidate_division_WPS-2", "Cyanobacteria", "Spirochaetes"))


emm_time_df <- do.call(rbind, emm_time_list) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) %>%
  dplyr::filter(!group %in% c("Unassigned","Thaumarchaeota","Chlamydiae", 
                              "BRC1", "candidate_division_WPS-2", "Cyanobacteria", "Spirochaetes"))


####----LMM new plot 20251017----####
group_order <- emm_all_df %>%
  dplyr::mutate(group = case_when(
    group == "candidate_division_WPS-1" ~ "Candidate_division_WPS-1",
    group == "candidate_division_WPS-2" ~ "Candidate_division_WPS-2",
    .default = group
  )) %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(estimate) %>%
  dplyr::pull(group)

# group_order

warming_df <- emm_all_df %>%
  dplyr::mutate(group = case_when(
    group == "candidate_division_WPS-1" ~ "Candidate_division_WPS-1",
    group == "candidate_division_WPS-2" ~ "Candidate_division_WPS-2",
    .default = group
  )) %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(estimate) %>%
  dplyr::mutate(change = case_when(
    p.value < 0.05 & estimate > 0 ~ "Positive",
    p.value < 0.05 & estimate < 0 ~ "Negative",
    .default = "Nochange"
  )) %>%
  dplyr::mutate(group = factor(group, levels = group_order, ordered = T))

p_warming_df <- ggplot(data = warming_df) + 
  # annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(x = estimate, y = group, color = change), size = 2.5) +
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low , xmax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = warming_df %>% dplyr::filter(change == "Positive"),
            mapping = aes(x = conf.high + 0.75, y = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  geom_text(data = warming_df %>% dplyr::filter(change == "Negative"),
            mapping = aes(x = conf.low - 0.5, y = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  scale_color_manual(values = c("Positive" = "#d7301f",
                                "Negative" = "#0570b0",
                                "Nochange" = "#525252")) + 
  scale_x_continuous(
    limits = c(-2, 2),
    expand = expansion(mult = c(0.1, 0.1)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming") + 
  labs(x = "Effect size", y = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 2),
    panel.background = element_rect(fill='#fff7ec',colour = 'black'),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3, color = c(rep("#0570b0", 6), 
                                                   rep("#525252",5), 
                                                   rep("#d7301f", 2), 
                                                   rep("#525252",2), 
                                                   rep("#d7301f", 5)), 
                               size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,1,1,1,"cm"),
    # aspect.ratio = 1.5
  )
  

# p_warming_df



p_warming_df2 <- ggplot(data = warming_df %>%
                          dplyr::mutate(group = factor(group, levels = rev(levels(group)), ordered = T))) + 
  # annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(y = estimate, x = group, color = change), size = 2.5) +
  geom_errorbar(aes(y = estimate, x = group, ymin = conf.low , ymax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = warming_df %>% dplyr::filter(change == "Positive"),
            mapping = aes(y = conf.high + 0.75, x = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  geom_text(data = warming_df %>% dplyr::filter(change == "Negative"),
            mapping = aes(y = conf.low - 0.5, x = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  scale_color_manual(values = c("Positive" = "#d7301f",
                                "Negative" = "#0570b0",
                                "Nochange" = "#525252")) + 
  scale_y_continuous(
    limits = c(-2, 2),
    # expand = expansion(mult = c(0.1, 0.1)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming") +
  labs(y = "Effect size", x = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 2),
    panel.background = element_rect(fill='#fff7ec',colour = 'black'),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1,
                               color = rev(c(rep("#0570b0", 5), 
                                             rep("#525252",2), 
                                             rep("#d7301f", 1), 
                                             rep("#525252",2), 
                                             rep("#d7301f", 4))),
                               face = 3),
    axis.text.y = element_text(color = "#000000", size = 15),
    # axis.text.x = element_text(face = 3, , 
    #                            size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,1,1,1,"cm"),
    # aspect.ratio = 1.5
  )


# p_warming_df2




warmingTime_df <- emm_all_df %>%
  dplyr::mutate(group = case_when(
    group == "candidate_division_WPS-1" ~ "Candidate_division_WPS-1",
    group == "candidate_division_WPS-2" ~ "Candidate_division_WPS-2",
    .default = group
  )) %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::arrange(estimate) %>%
  dplyr::mutate(change = case_when(
    p.value < 0.05 & estimate > 0 ~ "Positive",
    p.value < 0.05 & estimate < 0 ~ "Negative",
    .default = "Nochange"
  )) %>%
  dplyr::mutate(group = factor(group, levels = group_order, ordered = T))

p_warmingTime_df <- ggplot(data = warmingTime_df) +
  # annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(x = estimate, y = group, color = change), size = 2.5) +
  geom_errorbar(aes(x = estimate, y = group, xmin = conf.low , xmax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = warmingTime_df %>% dplyr::filter(change == "Positive"),
            mapping = aes(x = conf.high + 0.75, y = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  geom_text(data = warmingTime_df %>% dplyr::filter(change == "Negative"),
            mapping = aes(x = conf.low - 0.75, y = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  scale_color_manual(values = c("Positive" = "#d7301f",
                                "Negative" = "#0570b0",
                                "Nochange" = "#525252")) + 
  scale_x_continuous(
    limits = c(-2, 2),
    expand = expansion(mult = c(0.1, 0.1)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming x Time") + 
  labs(x = "Effect size", y = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 2),
    panel.background = element_rect(fill="#fff7ec",colour = 'black'),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3, color = c(rep("#0570b0", 4), 
                                                   rep("#525252",8), 
                                                   rep("#d7301f", 3), 
                                                   rep("#525252",3), 
                                                   rep("#d7301f", 2)), 
                               size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,1,1,1,"cm"),
    # aspect.ratio = 1.5
  )


# p_warmingTime_df


p_warmingTime_df2 <- ggplot(data = warmingTime_df %>%
                              dplyr::mutate(group = factor(group, levels = rev(levels(group)), ordered = T))) + 
  # annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(y = estimate, x = group, color = change), size = 2.5) +
  geom_errorbar(aes(y = estimate, x = group, ymin = conf.low , ymax = conf.high, color = change), width = 0.2, linewidth = 1) + 
  geom_text(data = warmingTime_df %>% dplyr::filter(change == "Positive"),
            mapping = aes(y = conf.high + 0.75, x = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  geom_text(data = warmingTime_df %>% dplyr::filter(change == "Negative"),
            mapping = aes(y = conf.low - 0.75, x = group, label = signif, color = change), size = 7.5, vjust = 0.75) + 
  scale_color_manual(values = c("Positive" = "#d7301f",
                                "Negative" = "#0570b0",
                                "Nochange" = "#525252")) + 
  scale_y_continuous(
    limits = c(-2, 2),
    # expand = expansion(mult = c(0.1, 0.1)),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2)
  ) + 
  ggtitle(label = "Warming x Time") +
  labs(y = "Effect size", x = "") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 2),
    panel.background = element_rect(fill="#fff7ec",colour = 'black'),
    axis.text.x = element_text(size = 15, color = rev(c(rep("#0570b0", 4), 
                                                        rep("#525252",4), 
                                                        rep("#d7301f", 2), 
                                                        rep("#525252",3), 
                                                        rep("#d7301f", 1))),
                               angle = 45, hjust = 1, face = 3),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(1,1,1,1,"cm"),
    # aspect.ratio = 1.5
  )


# p_warmingTime_df2


p_combine <- p_warming_df + p_warmingTime_df




p_combine2 <- p_warming_df2 +  p_warmingTime_df2

p_combine2


ggsave(filename = "11.LMM_16S_ITS/New_LMM_out_16S_Phylum_20251214.pdf",
       plot = p_combine2,
       height = 6,
       width = 17)


####----old LMM----####
colnames(Phylum_16S)
glimpse(Phylum_16S)

Phylum_16S_LMM_result <- list()
LMM_everytime_result_16S <- list()

Phylum_16S_LMM_result_relative <- list()
LMM_everytime_result_16S_relative <- list()

Phylum_16S_LMM_result_new <- list()
Phylum_16S_LMM_result_relative_new <- list()

Phylum_16S_id <- colnames(Phylum_16S)[6:29]

for (g in Phylum_16S_id) {
  print(g)
  # g = "Firmicutes"
  data_tmp <- Phylum_16S %>% 
    dplyr::select(all_of(c(g, "Time", "Treatment","Replicate"))) %>%
    dplyr::mutate(Value = .data[[g]])
  data_tmp
  
  data_tmp_relative <- Phylum_16S_relative %>% 
    dplyr::select(all_of(c(g, "Time", "Treatment","Replicate"))) %>%
    dplyr::mutate(Value = .data[[g]])
  data_tmp_relative
  
  if (sum(data_tmp$Value) == 0 | is.nan(sum(data_tmp$Value))) {
    next
  }
  
  fm1 <- lmer(Value~Treatment*Time+(1|Replicate),data=data_tmp)
  # fm1 <- lmerTest::lmer(Value~Treatment*Time+(1|Replicate),data=data_tmp)
  # fm1 <- lmer(Value~Treatment*block+(1|Time),data=data_tmp)
  # fm1 <- lmer(Value~Treatment*Time,data=data_tmp)
  presult<-car::Anova(fm1,type=2)

  
  coefs<-coef(summary(fm1))[ , "Estimate"] 
  names(coefs)<-paste0(names(coefs),".mean")
  
  SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
  names(SEvalues)<-paste0(names(SEvalues),".se")
  
  tvalues<-coef(summary(fm1))[ , "t value"] ##t values
  names(tvalues)<-paste0(names(tvalues),".t")
  
  chisqP<-c(presult[,1],presult[,3])
  names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
  
  result<-c(coefs,tvalues,SEvalues,chisqP)
  
  Phylum_16S_LMM_result[[g]] <- result
  
  data_tmp2 <- Phylum_16S %>% 
    dplyr::select(all_of(c(g, "Time", "Treatment","Replicate"))) %>%
    dplyr::mutate(Value = .data[[g]]) %>%
    dplyr::mutate(Time = as.factor(Time)) %>%
    dplyr::mutate(Treatment = factor(Treatment, levels = c("Warming","Control"), ordered = T))
  data_tmp2
  
  fm2 <- lme4::lmer(Value~Treatment*Time+(1|Replicate),data=data_tmp2) 
  
  
  # # method1
  # # 提取模型中每个时间点的Treatment的效应
  em_results <- emmeans(fm2, ~ Treatment * Time, by = "Time")
  # em_results
  summary(em_results, by = "Time")
  # 
  em_results_time <- emmeans(fm2, ~ Treatment | Time, by = "Time")
  summary(em_results_time)
  # 
  # # 
  # # # 进行比较
  contrast_results <- contrast(em_results,  method = "pairwise")
  contrast_results
  # # 
  contrast_results_df <- broom::tidy(contrast_results)
  # # 
  contrast_results_df <- contrast_results_df %>% 
    dplyr::mutate(Gene = g)
  
  LMM_everytime_result_16S[[g]] <- contrast_results_df
  
  # new result abosolute
  fm1 <- lmerTest::lmer(Value~Treatment*Time+(1|Replicate),data=data_tmp)
  
  # 系数表：固定效应参数的点估计、SE、t、df、p
  smry <- summary(fm1)$coefficients
  
  ci <- suppressMessages(confint(fm1, method = "Wald", parm = "beta_", level = 0.95))
  
  ci_df <- data.frame(term = rownames(ci), ci.lower = ci[,1], ci.upper = ci[,2], row.names = NULL)
  
  coef_df <- data.frame(
    term      = rownames(smry),
    estimate  = smry[, "Estimate"],
    std.error = smry[, "Std. Error"],
    statistic = smry[, "t value"],
    df        = smry[, "df"],              # 由 lmerTest 计算
    p.value   = smry[, "Pr(>|t|)"],
    row.names = NULL
  ) %>%
    left_join(ci_df, by = "term") %>%
    mutate(
      sig_by_p  = p.value < 0.05,
      sig_by_ci = !(ci.lower <= 0 & ci.upper >= 0)  # TRUE 表示 95%CI 不跨 0（与显著性一致）
    ) %>%
    dplyr::mutate(g = g)
  
  coef_df
  
  Phylum_16S_LMM_result_new[[g]] <- coef_df
  
  # new result relative
  fm1 <- lmerTest::lmer(Value~Treatment*Time+(1|Replicate),data=data_tmp_relative)
  
  # 系数表：固定效应参数的点估计、SE、t、df、p
  smry <- summary(fm1)$coefficients
  
  ci <- suppressMessages(confint(fm1, method = "Wald", parm = "beta_", level = 0.95))
  
  ci_df <- data.frame(term = rownames(ci), ci.lower = ci[,1], ci.upper = ci[,2], row.names = NULL)
  
  coef_df <- data.frame(
    term      = rownames(smry),
    estimate  = smry[, "Estimate"],
    std.error = smry[, "Std. Error"],
    statistic = smry[, "t value"],
    df        = smry[, "df"],              # 由 lmerTest 计算
    p.value   = smry[, "Pr(>|t|)"],
    row.names = NULL
  ) %>%
    left_join(ci_df, by = "term") %>%
    mutate(
      sig_by_p  = p.value < 0.05,
      sig_by_ci = !(ci.lower <= 0 & ci.upper >= 0)  # TRUE 表示 95%CI 不跨 0（与显著性一致）
    ) %>%
    dplyr::mutate(g = g)
  
  coef_df
  
  Phylum_16S_LMM_result_relative_new[[g]] <- coef_df
  
}

Phylum_16S_LMM_result_new_df <- do.call(rbind, Phylum_16S_LMM_result_new)
Phylum_16S_LMM_result_relative_new_df <- do.call(rbind, Phylum_16S_LMM_result_relative_new)

####----16S Phylum new plot 202510----####
Phylum_16S_LMM_result_new_df2 <- Phylum_16S_LMM_result_new_df %>%
  dplyr::filter(term %in% c("TreatmentWarming", "TreatmentWarming:Time")) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) 

# relative
Phylum_16S_LMM_result_relative_new_df2 <- Phylum_16S_LMM_result_relative_new_df %>%
  dplyr::filter(term %in% c("TreatmentWarming", "TreatmentWarming:Time", "Time")) %>%
  dplyr::mutate(signif = case_when(
    `p.value` < 0.05 & `p.value` > 0.01 ~ "*",
    `p.value` < 0.01 & `p.value` > 0.001 ~ "**",
    `p.value` < 0.001 ~ "***",
    .default = ""
  )) 

g_order <- Phylum_16S_LMM_result_relative_new_df2 %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(estimate) %>%
  dplyr::pull(g)

p1 <- Phylum_16S_LMM_result_relative_new_df2 %>%
  dplyr::filter(term == "TreatmentWarming:Time") %>%
  dplyr::arrange(desc(estimate)) %>%
  dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
  dplyr::mutate(change = case_when(
    estimate > 0 ~ "Up",
    estimate <= 0 ~ "Down"
  )) %>%
  ggplot() + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_point(aes(x = estimate, y = g, color = change), size = 2.5) + 
  geom_errorbar(aes(x = estimate, y = g, xmin = ci.lower, xmax = ci.upper, color = change), width = 0.2, linewidth = 1) +
  geom_text(data = Phylum_16S_LMM_result_relative_new_df2 %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>%
              dplyr::arrange(desc(estimate)) %>%
              dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
              dplyr::mutate(change = case_when(
                estimate > 0 ~ "Up",
                estimate <= 0 ~ "Down"
              )) %>%
              dplyr::filter(estimate > 0),
            mapping = aes(x = ci.upper + 0.2, y = g, label = signif, color = change), 
            size = 7.5,
            vjust = 0.7) +
  geom_text(data = Phylum_16S_LMM_result_relative_new_df2 %>%
              dplyr::filter(term == "TreatmentWarming:Time") %>%
              dplyr::arrange(desc(estimate)) %>%
              dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
              dplyr::mutate(change = case_when(
                estimate > 0 ~ "Up",
                estimate <= 0 ~ "Down"
              )) %>%
              dplyr::filter(estimate <= 0),
            mapping = aes(x = ci.lower - 0.2, y = g, label = signif, color = change), 
            size = 7.5,
            vjust = 0.7) +
  geom_vline(xintercept = 0) + 
  ggtitle(label = "Warming * Time") + 
  labs(x = "Effect Size", y = "") + 
  scale_color_manual(values = c("Up" = "#ef6548",
                                "Down" = "#3690c0"))+
  scale_x_continuous(
    limits = c(-1.0, 1.0),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c(-1, -0.5, 0, 0.5, 1)
  ) + 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 0.5),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3, color = c(rep("#41b6c4", 9),
                                                   rep("#41b6c4", 12)),
                               size = 15),
    axis.title = element_text(color = "#000000", size = 15),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 2.5
  ) + 
  coord_cartesian(clip = "off") +
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = "#7fcdbb", color = "#7fcdbb", lwd = 0, alpha = 0.2)),
                    xmin = -1.02,
                    xmax = -3,
                    ymin = 0.35,
                    ymax = 21.6) + 
  annotation_custom(grob = textGrob(label = "Bacterial\nphyla", gp = gpar(col = "#4292c6", cex = 1.5), rot = 90),
                    xmin = -2.8,
                    xmax = -2.8,
                    ymin = 12,
                    ymax = 12) 

p1

p2 <- Phylum_16S_LMM_result_relative_new_df2 %>%
  dplyr::filter(term == "TreatmentWarming") %>%
  dplyr::arrange(desc(estimate)) %>%
  dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
  dplyr::mutate(change = case_when(
    estimate > 0 ~ "Up",
    estimate <= 0 ~ "Down"
  )) %>%
  ggplot() + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#7fcdbb", alpha = 0.2) +
  geom_point(aes(x = estimate, y = g, color = change), size = 2.5) + 
  geom_errorbar(aes(x = estimate, y = g, xmin = ci.lower, xmax = ci.upper, color = change), width = 0.2, linewidth = 1) +
  geom_text(data = Phylum_16S_LMM_result_relative_new_df2 %>%
              dplyr::filter(term == "TreatmentWarming") %>%
              dplyr::arrange(desc(estimate)) %>%
              dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
              dplyr::mutate(change = case_when(
                estimate > 0 ~ "Up",
                estimate <= 0 ~ "Down"
              )) %>%
              dplyr::filter(estimate > 0),
            mapping = aes(x = ci.upper + 0.4, y = g, label = signif, color = change), 
            size = 7.5,
            vjust = 0.7) +
  geom_text(data = Phylum_16S_LMM_result_relative_new_df2 %>%
              dplyr::filter(term == "TreatmentWarming") %>%
              dplyr::arrange(desc(estimate)) %>%
              dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
              dplyr::mutate(change = case_when(
                estimate > 0 ~ "Up",
                estimate <= 0 ~ "Down"
              )) %>%
              dplyr::filter(estimate <= 0),
            mapping = aes(x = ci.lower - 0.4, y = g, label = signif, color = change), 
            size = 7.5,
            vjust = 0.7) +
  geom_vline(xintercept = 0) + 
  ggtitle(label = "Warming") + 
  labs(x = "Effect Size", y = "") + 
  scale_color_manual(values = c("Up" = "#ef6548",
                                "Down" = "#3690c0")) +
  scale_x_continuous(
    limits = c(-3, 3),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c(-3, -2, -1, 0, 1, 2, 3)
  ) + 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 0.5),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(color = "#000000", size = 15),
    axis.text.y = element_text(face = 3, color = c(rep("#41b6c4", 9),
                                                   rep("#41b6c4", 12)),
                               size = 15),
    axis.title = element_text(color = "#000000", size = 15),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 2.5
  ) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(grob = grid::rectGrob(gp = gpar(fill = "#7fcdbb", color = "#7fcdbb", lwd = 0, alpha = 0.2)),
                    xmin = -3.05,
                    xmax = -9,
                    ymin = 0.35,
                    ymax = 21.6) + 
  annotation_custom(grob = textGrob(label = "Bacterial\nphyla", gp = gpar(col = "#4292c6", cex = 1.5), rot = 90),
                    xmin = -8.35,
                    xmax = -8.35,
                    ymin = 12,
                    ymax = 12) 

p2

p3 <- Phylum_16S_LMM_result_relative_new_df2 %>%
  dplyr::filter(term == "Time") %>%
  dplyr::arrange(desc(estimate)) %>%
  dplyr::mutate(g = factor(g, levels = g_order, ordered = T)) %>%
  ggplot() + 
  geom_point(aes(x = estimate, y = g)) + 
  geom_errorbar(aes(x = estimate, y = g, xmin = ci.lower, xmax = ci.upper), width = 0.15) +
  geom_text(aes(x = ci.upper + 0.1, y = g, label = signif), size = 5) +
  geom_vline(xintercept = 0) + 
  ggtitle(label = "Time") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(color = "#000000", size = 12.5)
  )

p_combine <- p2 + 
  plot_spacer() + 
  p1 +  
  plot_layout(widths = c(1, 0.01, 1)) + 
  plot_annotation(
    tag_levels = c('A'),
    theme = theme(plot.tag = element_text(size = 10))
  )

p_combine


ggsave(filename = "11.LMM_16S_ITS/LMM_Phylum.pdf",
       plot = p_combine,
       height = 12.5,
       width = 20)

