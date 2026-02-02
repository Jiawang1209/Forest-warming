rm(list = ls())

####----load R Package----####
library(tidyverse)
library(ggfun)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(lme4)
library(car)
library(ggrepel)
library(ggh4x)
library(ggnewscale)
library(aplot)
library(broom)
library(sjPlot)
library(emmeans)
library(lmerTest)
library(linkET)
library(grid)
library(vegan)
library(psych)
library(linkET)
library(readxl)
library(psych)
library(ggnewscale)
library(broom)
library(sjPlot)
library(emmeans)
library(lmerTest)
library(corrplot)
library(grid)
source("R/cor_test.R")
set.seed(2024)

dir.create("9.Env_Mirco_Correlation")

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


####----load Data----####
#####-----Env scale data-----#####
Envdf <- read_delim(file = "Input_File/data_export_new_20250714.csv", col_names = T, delim = ",") %>%
  dplyr::filter(!str_detect(string = Name, pattern = "0-10")) %>%
  dplyr::filter(!Date1 %in% c("190705", "201105")) %>%
  dplyr::select(Year, Duration1, Treatment, Plot, 
                Water,
                DOC, TDN, DON, 
                NH4, NO3, DIN,
                TC, TN, TC_TN 
                ) %>%
  dplyr::select(Year, Duration1, Treatment, Plot, 
                Water,
                DOC, TC,
                TN, TDN, DON, DIN, NH4, NO3,
                TC_TN
                ) %>%
  dplyr::rename(Time = Duration1)

head(Envdf)

colnames(Envdf)

table(Envdf$Year)

Envdf



Envdf2 <- Envdf %>%
  dplyr::mutate(across(5:dim(Envdf)[2], as.numeric)) %>%
  dplyr::mutate(across(5:dim(Envdf)[2], ~scale2(.x, na.rm = TRUE)))

Envdf2

# 在这里，我也不做线性混合模型，我就是单纯的使用原始数据去做分析

####----alpha Diversity----####

Bactarial_alpha <- read_delim(file = "./1.Alpha_Diversity/alpha_result_O_Layer.csv", col_names = T, delim = ",") %>%
  dplyr::filter(!str_detect(string = Sample_month, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120)

Bactarial_alpha2 <- Bactarial_alpha %>%
  dplyr::mutate(across(2:9, as.numeric)) %>%
  dplyr::mutate(across(2:9, ~scale2(.x, na.rm = TRUE))) %>%
  dplyr::mutate(Plot = Envdf$Plot,
                Time = Envdf$Time,
                Year = Envdf$Year)

Bactarial_beta <- read_delim(file = "2.Beta_Diversity/Beta_Diversity_result_O_Layer.csv", col_names = T, delim = ",") %>%
  dplyr::select(-1) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1)


Bactarial_beta_pcoa <- read_delim(file = "2.Beta_Diversity/Beta_Diversity_result_O_Layer_pcoa1to10.csv", col_names = T, delim = ",") %>%
  dplyr::rename(Sample = `...1`) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1) %>%
  purrr::set_names(str_c("B_", str_c("PCOA", 1:10)))

Bactarial_beta_nmds <- read_delim(file = "2.Beta_Diversity/Beta_Diversity_result_O_Layer_NMDS.csv", col_names = T, delim = ",") %>%
  dplyr::rename(Sample = `...1`) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1) %>%
  purrr::set_names(str_c("B_", str_c("NMDS", 1:2)))



Fungle_alpha <- read_delim(file = "../ITS/1.Alpha_Diversity/alpha_result_O_Layer.csv", col_names = T, delim = ",")%>%
  dplyr::filter(!str_detect(string = Sample_month, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120)

Fungle_beta <- read_delim(file = "../ITS/2.Beta_Diversity/Beta_Diversity_result_O_Layer.csv", col_names = T, delim = ",") %>%
  dplyr::select(-1) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1)


Fungle_beta_pcoa <- read_delim(file = "../ITS/2.Beta_Diversity/Beta_Diversity_result_O_Layer_pcoa1to10.csv", col_names = T, delim = ",") %>%
  dplyr::rename(Sample = `...1`) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1) %>%
  purrr::set_names(str_c("F_", str_c("PCOA", 1:10)))

Fungle_beta_nmds <- read_delim(file = "../ITS/2.Beta_Diversity/Beta_Diversity_result_O_Layer_NMDS.csv", col_names = T, delim = ",") %>%
  dplyr::rename(Sample = `...1`) %>%
  dplyr::filter(!str_detect(string = Sample, pattern = "2019.7|2020.11")) %>%
  dplyr::slice(13:120) %>%
  dplyr::select(-1) %>%
  purrr::set_names(str_c("F_", str_c("NMDS", 1:2)))

# beta diversity
Envdf_clean <- Envdf %>% 
  dplyr::select(5:dim(Envdf)[2])

Bactarial_alpha_clean <- Bactarial_alpha %>% 
  dplyr::select(2:9) %>%
  purrr::set_names(str_c("B", colnames(.), sep = "_"))

Fungle_alpha_clean <- Fungle_alpha %>%
  dplyr::select(2:9) %>%
  purrr::set_names(str_c("F", colnames(.), sep = "_"))

BF_alpha_clean <- cbind(Bactarial_alpha_clean, Fungle_alpha_clean) %>%
  dplyr::select(1,4,6,2,3,5,7,8,
                9,12,14,10,11,13,15,16)


Bactarial_beta_clean <- Bactarial_beta %>%
  purrr::set_names(str_c("B", colnames(.), sep = "_"))
# 
Fungle_beta_clean <- Fungle_beta %>%
  purrr::set_names(str_c("F", colnames(.), sep = "_"))

BF_beta_clean <- cbind(Bactarial_beta_clean, Fungle_beta_clean)
BF_beta_clean_pcoa <- cbind(Bactarial_beta_pcoa, Fungle_beta_pcoa)
BF_beta_clean_nmds <- cbind(Bactarial_beta_nmds, Fungle_beta_nmds)

BF <- cbind(BF_alpha_clean, BF_beta_clean, BF_beta_clean_pcoa, BF_beta_clean_nmds)

BF_scale <- BF %>%
  dplyr::mutate(across(1:304, ~scale2(.x, na.rm = TRUE)))

Envdf_clean_scale <- Envdf_clean %>%
  dplyr::mutate(across(1:dim(Envdf_clean)[2], ~scale2(.x, na.rm = TRUE)))


####----warming----####
mantel <- mantel_test(BF %>% dplyr::slice(which(Envdf$Treatment == "warming")),
                      Envdf_clean %>% dplyr::slice(which(Envdf$Treatment == "warming")),
                      mantel_fun = "mantel.randtest",
                      spec_select = list(
                        `Bacterial diversity`= 2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

mantel

sum(mantel$p < 0.05)

mantel %>%
  dplyr::filter(p < 0.05)

tmp_cor_warming <- corr.test(BF %>% dplyr::slice(which(Envdf$Treatment == "warming")) %>%
                               dplyr::select(B_Shannon_1),
                             Envdf_clean %>% dplyr::slice(which(Envdf$Treatment == "warming"))
)

corrtest <- tmp_cor_warming$r %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "spec") %>%
  dplyr::mutate(spec = "Bacterial diversity") %>%
  tidyr::pivot_longer(cols = -spec, names_to = "env", values_to = "r") %>%
  dplyr::left_join(tmp_cor_warming$p %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "spec") %>%
                     dplyr::mutate(spec = "Bacterial diversity") %>%
                     tidyr::pivot_longer(cols = -spec, names_to = "env", values_to = "p"),
                   by = c("spec", "env")) %>%
  dplyr::mutate(rd = cut(r, 
                         breaks = c(-Inf,-0.2, 0, 0.4, Inf),
                         labels = c("Cor <= -0.2", 
                                    "-0.2 < Cor <= 0", 
                                    "0 < Cor <= 0.4", 
                                    "Cor >= 0.4")),
                pd = cut(p, breaks = c(-Inf, 0.05, Inf),
                         labels = c("< 0.05", ">= 0.05")))

corrtest


p <- qcorrplot(correlate(Envdf_clean%>% dplyr::slice(which(Envdf$Treatment == "warming"))),
               type = "upper", diag = F, grid_col = NA) +
  geom_tile(color = "#000000", fill = NA, linewidth = 0.85) +
  geom_point(aes(size = abs(r), fill = r), shape = 21) +
  # geom_square()
  geom_mark(sep = '', size = 10, only_mark = T, vjust = 0.75) + 
  scale_fill_gradient2(low = "#4d9221", high = "#c51b7d", name = "Pearson's r \n (Env-Env)") +
  scale_size(range = c(8, 25), name = "Pearson's r \n (Env-Env)") + 
  guides(size = "none") + 
  new_scale("size") +
  geom_couple(aes(colour = r, 
                  size = abs(r),
                  linetype = pd), 
              data = corrtest,
              label.size = 5,
              nudge_x = 0.1
              # curvature = nice_curvature()
  ) +
  scale_size(range = c(2,10),
             name = "Pearson's abs(r) \n (Diversity-Env)",
             breaks = c(0.1, 0.25, 0.5)) +
  scale_linetype_manual(values = c("< 0.05" = 1,
                                   ">= 0.05" = 2),
                        name = "Pearson's p \n (Diversity-Env)") +
  scale_colour_gradient2(low = "#2166ac",
                         mid = "white",
                         high = "#b2182b",
                         midpoint = 0,
                         name = "Pearson's r \n (Diversity-Env)",
                         limit = c(-0.8, 0.8)) +
  guides(size = guide_legend(title = "Pearson's r \n (Diversity-Env)",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_colorbar(title = "Pearson's r \n (Diversity-Env)", 
                               override.aes = list(size = 3), 
                               order = 1),
         linetype = guide_legend(order = 3),
         fill = guide_colorbar(title = "Pearson's r \n (Env-Env)", order = 4)
         ) + 
  theme(
    plot.margin = margin(10,10,10,20,unit = "pt"),
    legend.frame = element_rect(color = "#000000"),
    legend.ticks = element_line(color = "#000000")
  ) + 
  coord_fixed(clip = "off")

# p


p

ggsave(filename = "9.Env_Mirco_Correlation/New_mental_test_warming_20260120.pdf",
       plot = p,
       height = 12,
       width = 12)


####----control----####
mantel2 <- mantel_test(BF %>% dplyr::slice(which(Envdf$Treatment == "control")),
                      Envdf_clean %>% dplyr::slice(which(Envdf$Treatment == "control")),
                      mantel_fun = "mantel.randtest",
                      spec_select = list(
                        `Bacterial diversity`= 2
                      )) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

mantel2

sum(mantel2$p < 0.05)

mantel2 %>%
  dplyr::filter(p < 0.05)


tmp_cor_control <- corr.test(BF %>% dplyr::slice(which(Envdf$Treatment == "control")) %>%
                               dplyr::select(B_Shannon_1),
                             Envdf_clean %>% dplyr::slice(which(Envdf$Treatment == "control"))
)

corrtest2 <- tmp_cor_control$r %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "spec") %>%
  dplyr::mutate(spec = "Bacterial diversity") %>%
  tidyr::pivot_longer(cols = -spec, names_to = "env", values_to = "r") %>%
  dplyr::left_join(tmp_cor_control$p %>%
                     as.data.frame() %>%
                     tibble::rownames_to_column(var = "spec") %>%
                     dplyr::mutate(spec = "Bacterial diversity") %>%
                     tidyr::pivot_longer(cols = -spec, names_to = "env", values_to = "p"),
                   by = c("spec", "env")) %>%
  dplyr::mutate(rd = cut(r, 
                         breaks = c(-Inf,-0.2, 0, 0.4, Inf),
                         labels = c("Cor <= -0.2", 
                                    "-0.2 < Cor <= 0", 
                                    "0 < Cor <= 0.4", 
                                    "Cor >= 0.4")),
                pd = cut(p, breaks = c(-Inf, 0.05, Inf),
                         labels = c("< 0.05", ">= 0.05"))) %>%
  dplyr::mutate(pd = factor(pd, levels = c("< 0.05", ">= 0.05"), ordered = T))

corrtest2

p2 <- qcorrplot(correlate(Envdf_clean%>% dplyr::slice(which(Envdf$Treatment == "control"))),
               type = "upper", diag = F, grid_col = NA) +
  geom_tile(color = "#000000", fill = NA, linewidth = 0.85) +
  geom_point(aes(size = abs(r), fill = r), shape = 21) +
  # geom_square()
  geom_mark(sep = '', size = 10, only_mark = T, vjust = 0.75) + 
  scale_fill_gradient2(low = "#4d9221", high = "#c51b7d", name = "Pearson's r \n (Env-Env)") +
  scale_size(range = c(8, 25), name = "Pearson's r \n (Env-Env)") + 
  guides(size = "none") + 
  new_scale("size") +
  geom_couple(aes(colour = r, 
                  size = abs(r),
                  linetype = pd), 
              data = corrtest2,
              label.size = 5,
              nudge_x = 0.1
              # curvature = nice_curvature()
  ) +
  scale_size(range = c(1,6),
             name = "Pearson's abs(r) \n (Diversity-Env)",
             breaks = c(0.005,  0.01, 0.02)) +
  scale_linetype_manual(values = c("< 0.05" = 1,
                                   ">= 0.05" = 2),
                        drop = FALSE,
                        name = "Pearson's p \n (Diversity-Env)") +
  scale_colour_gradient2(low = "#2166ac",
                         mid = "white",
                         high = "#b2182b",
                         midpoint = 0,
                         name = "Pearson's r \n (Diversity-Env)",
                         limit = c(-0.6, 0.6)) +
  guides(size = guide_legend(title = "Pearson's r \n (Diversity-Env)",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_colorbar(title = "Pearson's r \n (Diversity-Env)", 
                                 override.aes = list(size = 3), 
                                 order = 1),
         linetype = guide_legend(order = 3),
         fill = guide_colorbar(title = "Pearson's r \n (Env-Env)", order = 4)
  ) + 
  theme(
    plot.margin = margin(10,10,10,20,unit = "pt"),
    legend.frame = element_rect(color = "#000000"),
    legend.ticks = element_line(color = "#000000")
  ) + 
  coord_fixed(clip = "off")

p2


ggsave(filename = "9.Env_Mirco_Correlation/New_mental_test_control_20260120.pdf",
       plot = p2,
       height = 12,
       width = 12)
