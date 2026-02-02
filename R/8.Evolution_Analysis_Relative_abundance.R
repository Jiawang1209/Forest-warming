rm(list = ls())

####----load R Package----####
library(tidyverse)
library(readxl)
library(vegan)
library(ggalt)
library(patchwork)
library(ggpmisc)

dir.create("8.Tax_OTU_Abundance")

####----load Data----####
####----load Phylum relative abundance----####
Phylum_relative <- read_delim(file = "Input_File/O_Layer_16S_Phylum_relative_abundance.csv",
                              col_names = T, 
                              delim = ",") %>%
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
  dplyr::select(-1) %>%
  purrr::set_names(str_to_sentence(colnames(.))) %>%
  dplyr::select(-Thaumarchaeota)


####----ratio----####
relative_abundance_line_ratio_list <- list()

formula <- y ~ x

for (p in colnames(Phylum_relative)[5:24]) {
  print(p)
  # relative abundance line
  plot <- Phylum_relative %>%
    dplyr::select(Time, Treatment, all_of(p)) %>%
    dplyr::group_split(Time) %>%
    purrr::map(., function(x){
      dplyr::bind_cols(
        x %>% dplyr::slice(1:6),
        x %>% dplyr::slice(7:12)
      ) %>%
        purrr::set_names(c("Time1","Treatment1", "Value1",
                           "Time2","Treatment2", "Value2"))
    }) %>%
    do.call(rbind, .) %>%
    dplyr::select(Time1, Value1, Value2) %>%
    purrr::set_names(c("Time", "Control", "Warming")) %>%
    dplyr::mutate(Ratio = Warming/Control) %>%
    dplyr::mutate(Phylum = p) %>%
    ggplot(aes(x = Time, y = Ratio)) + 
    geom_point(aes(x = Time, y = Ratio, colour = Phylum), size = 5, alpha = 0.75) +
    geom_smooth(
      method = "lm", formula = formula, se = F, linewidth = 3) + 
    stat_poly_eq(
      use_label(c("R2", "p")),formula = formula, size = 5) + 
    scale_color_manual(values = c("Acidobacteria" = "#266d9f",
                                 "Actinobacteria" = "#8a1c41",
                                 "Armatimonadetes" = "grey",
                                 "Brc1" = "grey",
                                 "Bacteroidetes" = "#c25f97",
                                 "Candidatus_saccharibacteria" = "#7ca5c2",
                                 "Chlamydiae" = "grey",
                                 "Chloroflexi" = "grey",
                                 "Cyanobacteria" = "grey",
                                 "Firmicutes" = "#df9228",
                                 "Gemmatimonadetes" = "#d43e2f",
                                 "Latescibacteria" = "#3a9239",
                                 "Nitrospirae" = "grey",
                                 "Planctomycetes" = "#b51e4f",
                                 "Proteobacteria" = "#d96521",
                                 "Spirochaetes" = "#266d9f",
                                 "Unassigned" = "grey",
                                 "Verrucomicrobia" = "#e37920",
                                 "Candidate_division_wps-1" = "#c73379",
                                 "Candidate_division_wps-2" = "grey"
                                 )) + 
    ggtitle(label = p) + 
    labs(x = "Year",
         y = "Relative abundance ratio (Warming/Control)") + 
    theme_classic() + 
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(color = "#000000", size = 15),
      axis.title = element_text(color = "#000000", size = 17.5),
      plot.title = element_text(color = "#000000", size = 17.5, hjust = 0.5, face = "italic"),
      plot.margin = margin(5,5,5,5,unit = "pt"),
      axis.line = element_line(color = "#000000", linewidth = 1),
      aspect.ratio = 0.85
    )
  
  relative_abundance_line_ratio_list[[p]] <- plot
}


colnames(Phylum_relative)

Phylum_relative %>%
  dplyr::select(Time, Treatment, all_of("Candidate_division_wps-2")) %>%
  dplyr::mutate(`Candidate_division_wps-2` = `Candidate_division_wps-2` + 1e-6) %>%
  dplyr::group_split(Time) %>%
  purrr::map(., function(x){
    dplyr::bind_cols(
      x %>% dplyr::slice(1:6),
      x %>% dplyr::slice(7:12)
    ) %>%
      purrr::set_names(c("Time1","Treatment1", "Value1",
                         "Time2","Treatment2", "Value2"))
  }) %>%
  do.call(rbind, .) %>%
  dplyr::select(Time1, Value1, Value2) %>%
  purrr::set_names(c("Time", "Control", "Warming")) %>%
  dplyr::mutate(Ratio = Warming/Control) %>%
  dplyr::mutate(Phylum = "Candidate_division_wps-2") %>%
  ggplot(aes(x = Time, y = Ratio)) + 
  geom_point(aes(x = Time, y = Ratio, colour = Phylum), size = 5, alpha = 0.75) +
  geom_smooth(
    method = "lm", formula = formula, se = F, linewidth = 3) + 
  stat_poly_eq(
    use_label(c("R2", "p")),formula = formula, size = 5) + 
  scale_color_manual(values = c("Acidobacteria" = "#266d9f",
                                "Actinobacteria" = "#8a1c41",
                                "Armatimonadetes" = "grey",
                                "Brc1" = "grey",
                                "Bacteroidetes" = "#c25f97",
                                "Candidatus_saccharibacteria" = "#7ca5c2",
                                "Chlamydiae" = "grey",
                                "Chloroflexi" = "grey",
                                "Cyanobacteria" = "grey",
                                "Firmicutes" = "#df9228",
                                "Gemmatimonadetes" = "#d43e2f",
                                "Latescibacteria" = "#3a9239",
                                "Nitrospirae" = "grey",
                                "Planctomycetes" = "#b51e4f",
                                "Proteobacteria" = "#d96521",
                                "Spirochaetes" = "#266d9f",
                                "Unassigned" = "grey",
                                "Verrucomicrobia" = "#e37920",
                                "Candidate_division_wps-1" = "#c73379",
                                "Candidate_division_wps-2" = "grey"
  )) + 
  ggtitle(label = "Candidate_division_wps-2") + 
  labs(x = "Year",
       y = "Relative abundance ratio (Warming/Control)") + 
  theme_classic() + 
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "#000000", size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 17.5, hjust = 0.5, face = "italic"),
    plot.margin = margin(5,5,5,5,unit = "pt"),
    axis.line = element_line(color = "#000000", linewidth = 1),
    aspect.ratio = 0.85
  ) -> tmp1


Phylum_relative %>%
  dplyr::select(Time, Treatment, all_of("Spirochaetes")) %>%
  dplyr::mutate(`Spirochaetes` = `Spirochaetes` + 1e-6) %>%
  dplyr::group_split(Time) %>%
  purrr::map(., function(x){
    dplyr::bind_cols(
      x %>% dplyr::slice(1:6),
      x %>% dplyr::slice(7:12)
    ) %>%
      purrr::set_names(c("Time1","Treatment1", "Value1",
                         "Time2","Treatment2", "Value2"))
  }) %>%
  do.call(rbind, .) %>%
  dplyr::select(Time1, Value1, Value2) %>%
  purrr::set_names(c("Time", "Control", "Warming")) %>%
  dplyr::mutate(Ratio = Warming/Control) %>%
  dplyr::mutate(Phylum = "Spirochaetes") %>%
  ggplot(aes(x = Time, y = Ratio)) + 
  geom_point(aes(x = Time, y = Ratio, colour = Phylum), size = 5, alpha = 0.75) +
  geom_smooth(
    method = "lm", formula = formula, se = F, linewidth = 3) + 
  stat_poly_eq(
    use_label(c("R2", "p")),formula = formula, size = 5) + 
  scale_color_manual(values = c("Acidobacteria" = "#266d9f",
                                "Actinobacteria" = "#8a1c41",
                                "Armatimonadetes" = "grey",
                                "Brc1" = "grey",
                                "Bacteroidetes" = "#c25f97",
                                "Candidatus_saccharibacteria" = "#7ca5c2",
                                "Chlamydiae" = "grey",
                                "Chloroflexi" = "grey",
                                "Cyanobacteria" = "grey",
                                "Firmicutes" = "#df9228",
                                "Gemmatimonadetes" = "#d43e2f",
                                "Latescibacteria" = "#3a9239",
                                "Nitrospirae" = "grey",
                                "Planctomycetes" = "#b51e4f",
                                "Proteobacteria" = "#d96521",
                                "Spirochaetes" = "#266d9f",
                                "Unassigned" = "grey",
                                "Verrucomicrobia" = "#e37920",
                                "Candidate_division_wps-1" = "#c73379",
                                "Candidate_division_wps-2" = "grey"
  )) + 
  ggtitle(label = "Spirochaetes") + 
  labs(x = "Year",
       y = "Relative abundance ratio (Warming/Control)") + 
  theme_classic() + 
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "#000000", size = 15),
    axis.title = element_text(color = "#000000", size = 17.5),
    plot.title = element_text(color = "#000000", size = 17.5, hjust = 0.5, face = "italic"),
    plot.margin = margin(5,5,5,5,unit = "pt"),
    axis.line = element_line(color = "#000000", linewidth = 1),
    aspect.ratio = 0.85
  ) -> tmp2

tmp2

p_combine_tmp <- tmp1 + tmp2 + plot_layout(guides = "collect")


relative_abundance_line_ratio_list[["Spirochaetes"]] <- tmp2
relative_abundance_line_ratio_list[["Candidate_division_wps-2"]] <- tmp1

patchwork::wrap_plots(relative_abundance_line_ratio_list,
                      nrow = 5,
                      guides = "collect") -> p_combine2

ggsave(filename = "8.Tax_OTU_Abundance/Phylum_year_ratio_line_add_point_line.pdf",
       plot = p_combine2,
       height = 20,
       width = 22,
       device = cairo_pdf)

####----ratio smooth----####
relative_abundance_line_ratio_list_smooth <- list()

# formula <- y ~ x

for (p in colnames(Phylum_relative)[5:24]) {
  print(p)
  # relative abundance line
  plot <- Phylum_relative %>%
    dplyr::select(Time, Treatment, all_of(p)) %>%
    dplyr::mutate(dplyr::across(all_of(p), ~ .x + 1e-6)) %>% 
    dplyr::group_split(Time) %>%
    purrr::map(., function(x){
      dplyr::bind_cols(
        x %>% dplyr::slice(1:6),
        x %>% dplyr::slice(7:12)
      ) %>%
        purrr::set_names(c("Time1","Treatment1", "Value1",
                           "Time2","Treatment2", "Value2"))
    }) %>%
    do.call(rbind, .) %>%
    dplyr::select(Time1, Value1, Value2) %>%
    purrr::set_names(c("Time", "Control", "Warming")) %>%
    dplyr::mutate(Ratio = Warming/Control) %>%
    dplyr::mutate(Phylum = p) %>%
    ggplot(aes(x = Time, y = Ratio)) + 
    geom_point(aes(x = Time, y = Ratio, colour = Phylum), size = 5, alpha = 0.75) +
    geom_smooth(se = F, linewidth = 3) + 
    # stat_poly_eq(
    #   use_label(c("R2", "p")),formula = formula, size = 5) + 
    scale_color_manual(values = c("Acidobacteria" = "#266d9f",
                                  "Actinobacteria" = "#8a1c41",
                                  "Armatimonadetes" = "#349189",
                                  "Brc1" = "#b81f75",
                                  "Bacteroidetes" = "#c25f97",
                                  "Candidatus_saccharibacteria" = "#7ca5c2",
                                  "Chlamydiae" = "#b81f75",
                                  "Chloroflexi" = "#54b332",
                                  "Cyanobacteria" = "#54b332",
                                  "Firmicutes" = "#df9228",
                                  "Gemmatimonadetes" = "#d43e2f",
                                  "Latescibacteria" = "#3a9239",
                                  "Nitrospirae" = "#4f2580",
                                  "Planctomycetes" = "#b51e4f",
                                  "Proteobacteria" = "#d96521",
                                  "Spirochaetes" = "#266d9f",
                                  # "Unassigned" = "grey",
                                  "Verrucomicrobia" = "#e37920",
                                  "Candidate_division_wps-1" = "#c73379",
                                  "Candidate_division_wps-2" = "#b81c25"
    )) + 
    ggtitle(label = p) + 
    labs(x = "Year",
         y = "Relative abundance ratio (Warming/Control)") + 
    theme_classic() + 
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(color = "#000000", size = 15),
      axis.title = element_text(color = "#000000", size = 17.5),
      plot.title = element_text(color = "#000000", size = 17.5, hjust = 0.5, face = "italic"),
      plot.margin = margin(5,5,5,5,unit = "pt"),
      axis.line = element_line(color = "#000000", linewidth = 1),
      aspect.ratio = 0.85
    )
  
  relative_abundance_line_ratio_list_smooth[[p]] <- plot
}


patchwork::wrap_plots(relative_abundance_line_ratio_list_smooth,
                      nrow = 5,
                      guides = "collect") -> p_combine_smooth

ggsave(filename = "8.Tax_OTU_Abundance/Phylum_year_ratio_line_add_point_smooth.pdf",
       plot = p_combine_smooth,
       height = 20,
       width = 22,
       device = cairo_pdf)


####----sessionInfo----####
sessionInfo()
