rm(list = ls())

####----load R Package----####
library(tidyverse)
library(ggpmisc)
library(ggrepel)
library(cowplot)
library(rstatix)
library(ggfun)
library(patchwork)
library(readxl)
library(lubridate)
library(lme4)
library(grid)
library(emmeans)
library(lmerTest)

dir.create("10.Diversity_respiration_relationship")

####----load Data----####
soil_respiration <- read_delim(file = "Input_File/20250623 soil CO2_dmean_mgCm-2h-1.csv",
                               delim = ",",
                               col_names = T) %>%
  dplyr::mutate(date = as_date(date)) %>%
  dplyr::mutate(
    Treatment = case_when(
      plot %in% c(1,4,5) & date >= as_date("2018-09-10") & date <= as_date("2018-11-15") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2019-06-08") & date <= as_date("2019-11-30") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-06-16") & date <= as_date("2020-08-26") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-09-12") & date <= as_date("2020-12-06") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2021-03-30") & date <= as_date("2021-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2022-03-26") & date <= as_date("2022-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2023-03-20") & date <= as_date("2023-12-07") ~ "Warming",
      .default = "Control"
    )
  ) %>%
  dplyr::filter(
    date >= as_date("2018-09-10") & date <= as_date("2018-11-15") |
      date >= as_date("2019-06-08") & date <= as_date("2019-11-30") |
      date >= as_date("2020-06-16") & date <= as_date("2020-08-26") |
      date >= as_date("2020-09-12") & date <= as_date("2020-12-06") |
      date >= as_date("2021-03-30") & date <= as_date("2021-12-07") |
      date >= as_date("2022-03-26") & date <= as_date("2022-12-07") |
      date >= as_date("2023-03-20") & date <= as_date("2023-12-07") 
  ) %>%
  dplyr::filter(
    year %in% 2019:2023,
    month(date) %in% c(6:9)
  ) %>%
  dplyr::mutate(`value_mgC m-2day-1` = `value_mgC m-2day-1` / (1000*1000*1000) * 10000) %>%
  dplyr::group_by(year, plot) %>%
  dplyr::summarise(sum = sum(`value_mgC m-2day-1`, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Treatment = case_when(
    plot %in% c(2,3,6) ~ "Control" ,
    plot %in% c(1,4,5) ~ "Warming"
  )) %>%
  dplyr::arrange(plot, year, Treatment) %>%
  dplyr::select(plot, year, Treatment, sum) %>%
  dplyr::rename(value = sum)


table(soil_respiration$year,
      soil_respiration$Treatment)


####----LMM CO2----####
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

soil_respiration_4_LMM <- soil_respiration %>%
  dplyr::arrange(year, Treatment, plot) %>%
  dplyr::mutate(Time = case_when(
    year == 2019 ~ 1,
    year == 2020 ~ 2,
    year == 2021 ~ 3,
    year == 2022 ~ 4,
    year == 2023 ~ 5
  )) %>%
  dplyr::mutate(across(Time, ~ scale2(.x, na.rm = T))) %>%
  dplyr::mutate(across(value, ~ scale2(.x, na.rm = T))) %>%
  dplyr::mutate(plot =factor(plot)) %>%
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))

fm1 <- lmerTest::lmer(value~Treatment * Time +(1|plot),data=soil_respiration_4_LMM)

out1 <- broom.mixed::tidy(fm1, effects = "fixed", conf.int = TRUE)

out1

# bar plot
soil_respiration %>%
  dplyr::group_by(year, Treatment) %>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value),
                   se = sd/sqrt(n())
                   ) %>%
  dplyr::ungroup() %>%
  ggplot() + 
  geom_bar(aes(x = year, y = mean, group = Treatment, fill = Treatment),
           stat = "identity", position = "dodge",
           color = "#000000",
           linewidth = 0.5) + 
  geom_errorbar(aes(x = year, y = mean, ymin = mean-se, ymax = mean+se,
                    group = Treatment),
                position = position_dodge(0.85),
                width = 0.25) + 
  scale_fill_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) + 
  labs(x = "Year", y = "Respiration") + 
  theme_bw() + 
  theme(panel.background = element_rect(colour = 'black', linewidth = 0.75),
        axis.title.x=element_text(colour = 'black',size=15),
        axis.title.y = element_text(colour = 'black',size=15),
        axis.text = element_text(colour = "black", size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.background = element_roundrect(color = "#808080",linetype = 1),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 20),
        panel.grid = element_blank(),
        # aspect.ratio = 1
  ) + 
  annotation_custom(
    grob = grid::textGrob(label = expression("Warming: " * italic(P) == 0.781),
                          gp = grid::gpar(col = "#000000", fontsize = 15)),
    xmin = 2020,
    xmax = 2020,
    ymin = 5.5,
    ymax = 5.5

  ) +
  annotation_custom(
    grob = grid::textGrob(label = expression("Time: " * italic(P) == 0.028),
                          gp = grid::gpar(col = "#000000", fontsize = 15)),
    xmin = 2019.9,
    xmax = 2019.9,
    ymin = 5.3,
    ymax = 5.3
  ) + 
  annotation_custom(
    grob = grid::textGrob(label = expression("Warming x Time: " * italic(P) == 0.988),
                          gp = grid::gpar(col = "#000000", fontsize = 15)),
    xmin = 2020.2,
    xmax = 2020.2,
    ymin = 5.1,
    ymax = 5.1
  )

ggsave(filename = "10.Diversity_respiration_relationship/soil_respiration_year_6-9_2.pdf",
       height = 7.5,
       width = 10)

####----寻找关键时期，增温增加呼吸，实质上就是增温导致CO2增加----####

soil_respiration2 <- read_delim(file = "Input_File/20250623 soil CO2_dmean_mgCm-2h-1.csv",
                               delim = ",",
                               col_names = T) %>%
  dplyr::mutate(date = as_date(date)) %>%
  dplyr::mutate(
    Treatment = case_when(
      plot %in% c(1,4,5) & date >= as_date("2018-09-10") & date <= as_date("2018-11-15") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2019-06-08") & date <= as_date("2019-11-30") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-06-16") & date <= as_date("2020-08-26") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-09-12") & date <= as_date("2020-12-06") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2021-03-30") & date <= as_date("2021-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2022-03-26") & date <= as_date("2022-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2023-03-20") & date <= as_date("2023-12-07") ~ "Warming",
      .default = "Control"
    )
  ) %>%
  dplyr::filter(
    year %in% 2019:2023,
    month(date) %in% c(7:8)
  ) %>%
  dplyr::mutate(`value_mgC m-2day-1` = `value_mgC m-2day-1` / (1000*1000*1000) * 10000) %>%
  dplyr::group_by(year, plot) %>%
  dplyr::summarise(sum = sum(`value_mgC m-2day-1`, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Treatment = case_when(
    plot %in% c(2,3,6) ~ "Control" ,
    plot %in% c(1,4,5) ~ "Warming"
  )) %>%
  dplyr::arrange(plot, year, Treatment) %>%
  dplyr::select(plot, year, Treatment, sum) %>%
  dplyr::rename(value = sum)

soil_respiration2 %>%
  dplyr::group_by(year, Treatment) %>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value),
                   se = sd/sqrt(n())
  ) %>%
  dplyr::ungroup() %>%
  ggplot() + 
  geom_bar(aes(x = year, y = mean, group = Treatment, fill = Treatment),
           stat = "identity", position = "dodge",
           color = "#000000",
           linewidth = 0.5) + 
  geom_errorbar(aes(x = year, y = mean, ymin = mean-se, ymax = mean+se,
                    group = Treatment),
                position = position_dodge(0.85),
                width = 0.25) + 
  scale_fill_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) + 
  labs(x = "Year", y = "Respiration") + 
  theme_bw() + 
  theme(panel.background = element_rect(colour = 'black', linewidth = 0.75),
        axis.title.x=element_text(colour = 'black',size=15),
        axis.title.y = element_text(colour = 'black',size=15),
        axis.text = element_text(colour = "black", size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.background = element_roundrect(color = "#808080",linetype = 1),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 20),
        panel.grid = element_blank(),
        # aspect.ratio = 1
  )

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

soil_respiration_4_LMM <- soil_respiration2 %>%
  dplyr::arrange(year, Treatment, plot) %>%
  dplyr::mutate(Time = case_when(
    year == 2019 ~ 1,
    year == 2020 ~ 2,
    year == 2021 ~ 3,
    year == 2022 ~ 4,
    year == 2023 ~ 5
  )) %>%
  dplyr::mutate(across(Time, ~ scale2(.x, na.rm = T))) %>%
  dplyr::mutate(across(value, ~ scale2(.x, na.rm = T))) %>%
  dplyr::mutate(plot =factor(plot)) %>%
  dplyr::mutate(Treatment = factor(Treatment)) %>%
  dplyr::mutate(Treatment = relevel(Treatment, ref = "Control"))

fm1 <- lmerTest::lmer(value~Treatment * Time +(1|plot),data=soil_respiration_4_LMM)

out1 <- broom.mixed::tidy(fm1, effects = "fixed", conf.int = TRUE)

out1

####----2024年5月， 2024年9月份----####
soil_respiration_2024 <- read_delim(file = "Input_File/20250623 soil CO2_dmean_mgCm-2h-1.csv",
                               delim = ",",
                               col_names = T) %>%
  dplyr::mutate(date = as_date(date)) %>%
  dplyr::mutate(
    Treatment = case_when(
      plot %in% c(1,4,5) & date >= as_date("2018-09-10") & date <= as_date("2018-11-15") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2019-06-08") & date <= as_date("2019-11-30") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-06-16") & date <= as_date("2020-08-26") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2020-09-12") & date <= as_date("2020-12-06") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2021-03-30") & date <= as_date("2021-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2022-03-26") & date <= as_date("2022-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2023-03-20") & date <= as_date("2023-12-07") ~ "Warming",
      plot %in% c(1,4,5) & date >= as_date("2024-03-26") & date <= as_date("2024-12-07") ~ "Warming",
      .default = "Control"
    )
  ) %>%
  dplyr::filter(
    date >= as_date("2018-09-10") & date <= as_date("2018-11-15") |
      date >= as_date("2019-06-08") & date <= as_date("2019-11-30") |
      date >= as_date("2020-06-16") & date <= as_date("2020-08-26") |
      date >= as_date("2020-09-12") & date <= as_date("2020-12-06") |
      date >= as_date("2021-03-30") & date <= as_date("2021-12-07") |
      date >= as_date("2022-03-26") & date <= as_date("2022-12-07") |
      date >= as_date("2023-03-20") & date <= as_date("2023-12-07") |
      date >= as_date("2024-03-26") & date <= as_date("2024-12-07") 
  ) %>%
  dplyr::filter(
    year %in% 2024,
    month(date) %in% c(5:9)
  ) %>%
  dplyr::mutate(`value_mgC m-2day-1` = `value_mgC m-2day-1` / (1000*1000*1000) * 10000) %>%
  dplyr::mutate(month = month(date))
soil_respiration_2024

ggplot(data = soil_respiration_2024) + 
  geom_line(aes(x = date, y = `value_mgC m-2day-1`,
                color = Treatment)) + 
  facet_wrap(~month, scale = "free") + 
  scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a"))


soil_respiration_2024 %>%
  dplyr::filter(month %in% c(5, 9))  %>%
  dplyr::mutate(month = case_when(
    month == 5 ~ "May",
    month == 9 ~ "September"
  )) %>%
  ggplot() + 
  geom_line(aes(x = date, y = `value_mgC m-2day-1`,
                color = Treatment)) + 
  geom_point(aes(x = date, y = `value_mgC m-2day-1`,
                 color = Treatment)) + 
  facet_wrap(~month, scale = "free", nrow = 2) + 
  scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
  ylab(expression("CO"[2]~"flux"~"(t C ha"^{-1}*")")) + 
  xlab("Month") + 
  theme_bw() + 
  theme(panel.background = element_rect(colour = 'black', linewidth = 0.75),
        axis.title.x=element_text(colour = 'black',size=25),
        axis.title.y = element_text(colour = 'black',size=25),
        axis.text = element_text(colour = "black", size = 15),
        plot.title = element_text(hjust = 0.5, size = 25),
        legend.background = element_roundrect(color = "#808080",linetype = 1),
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 20),
        panel.grid = element_blank(),
        aspect.ratio = 0.85,
        legend.position = "none"
  )

ggsave(filename = "10.Diversity_respiration_relationship/CO2_2024_5_9.pdf",
       height = 10,
       width = 5)

####----sessionInfo----####
sessionInfo()
