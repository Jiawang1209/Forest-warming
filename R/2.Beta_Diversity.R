####----Project:2.Beta Diversity----####
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
library(rstatix)
library(GGally)

dir.create("2.Beta_Diversity")

set.seed(1115)

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
otu_rare = as.data.frame(t(rrarefy(t(seqtab), 22000)))

colSums(otu_rare)
colSums(seqtab)

#####-----PCoA month Levels-----#####
pcoa_out_list <- list()
# pcoa
for (month in unique(str_split(colnames(otu_rare), pattern = "_", simplify = T)[,1])) {
  print(month)
  # month = "2020.6"
  otu_rare_sub <- seqtab %>%
    dplyr::select(starts_with(month))
  
  otu_rare_sub <- otu_rare_sub[which(rowSums(otu_rare_sub) != 0), ]
  
  group <- data.frame(
    sample = colnames(otu_rare_sub),
    group = str_remove(string = colnames(otu_rare_sub), pattern = "_\\d+") %>%
      str_remove(string = ., pattern = ".*_")
                      )
  
  # PERMANOVA统计检验
  adonis_result <- adonis2(t(otu_rare_sub) ~ group, data = group,permutations = 999,method="bray")
  adonis_result
  adonis_result$R2
  adonis_result$`Pr(>F)`
  
  dune_adonis <- paste0("adonis R2: ",round(adonis_result$R2[1],2), "; P-value: ", adonis_result$`Pr(>F)`[1])
  dune_adonis
  
  distance<-vegdist(t(otu_rare_sub), method='bray')
  
  as.matrix(distance) %>% as.data.frame() %>% write.csv(file = paste0("./2.Beta_Diversity/",month, "beta_bray.csv"))
  
  pcoa<- cmdscale(distance,k=(nrow(t(otu_rare_sub))-1),eig=TRUE)
  
  plot_data<-data.frame({pcoa$point})[1:2]
  
  names(plot_data)[1:2]<-c('PCoA1','PCoA2')
  
  eig=pcoa$eig
  
  plot_data$group = str_remove(rownames(plot_data), pattern = "_\\d+")
  plot_data$sample <- rownames(plot_data)
  plot_data$Treatment = str_remove(string = plot_data$group, pattern = ".*_")
  plot_data$Data = str_remove(string = plot_data$group, pattern = "_.*")
  
  
  sig = paste0("P = ", adonis_result$`Pr(>F)`[1])
  
  sig
  
  p_pcoa <- ggplot(plot_data, aes(x=PCoA1,y=PCoA2, color=Treatment))+
    # geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    stat_ellipse(level=0.75,size=1)+
    # ggforce::geom_mark_ellipse(aes(x = PCoA1, y = PCoA2,color = Treatment, fill = Treatment), 
    #                            expand = unit(1, "mm"),
    #                            alpha = 0.5) + 
    geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    # geom_text_repel(aes(x = PCoA1, y = PCoA2, label = sample)) +
    # labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),
    #      y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep="")) +
    labs(x = "", y = "") + 
    # ggtitle(label = paste("Group:", month, "\n", dune_adonis)) + 
    # facet_wrap(~Data) + 
    xlim(c(-0.45, 0.45)) +
    ylim(c(-0.4, 0.4)) +
    scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
    scale_fill_manual(values = c("Control"= "#a6bddb", "Warming" = "#fdbb84")) + 
    geom_vline(aes(xintercept=0),linetype="dotted",linewidth = 1.25)+
    geom_hline(aes(yintercept=0),linetype="dotted",linewidth = 1.25)+
    coord_fixed() + 
    theme_bw() + 
    theme(panel.background = element_rect(fill='#fff7ec',colour = 'black'),
          panel.grid = element_blank(),
          axis.title.x=element_text(colour = 'black',size=10),
          axis.title.y = element_text(colour = 'black',size=10),
          plot.title = element_text(hjust = 0.5, size = 20),
          # axis.text = element_text(colour = "black", size = 10),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.background = element_roundrect(color = "#808080",linetype = 1),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 12),
          # strip.text = element_text(size = 25),
          # strip.background = element_rect(color = "black"),
          legend.position = "none",
          panel.border = element_rect(linewidth = 1.5)
          ) + 
    # annotate(geom = "text", x = -0.25, y = 0.35, label = sig, size = 3) 
    annotate(geom = "text", x = -0.22, y = 0.35, 
             label = paste0("italic(P) == ", adonis_result$`Pr(>F)`[1]), 
             size = 4,
             parse = TRUE)
  
  p_pcoa
  
  pcoa_out_list[[month]] <- p_pcoa
}

pcoa_out_list[[1]]
pcoa_out_list[[2]]
pcoa_out_list[[3]]
pcoa_out_list[[4]]
pcoa_out_list[[5]]
pcoa_out_list[[6]]
pcoa_out_list[[7]]
pcoa_out_list[[8]]
pcoa_out_list[[9]]
pcoa_out_list[[10]]
pcoa_out_list[[11]]


# 拼图
p_pcoa_combine <- pcoa_out_list[[1]] + pcoa_out_list[[2]] +  pcoa_out_list[[3]] + pcoa_out_list[[4]] + pcoa_out_list[[5]] +
  pcoa_out_list[[6]] + pcoa_out_list[[7]] + pcoa_out_list[[8]] + pcoa_out_list[[9]] +
  pcoa_out_list[[10]] + pcoa_out_list[[11]] + plot_layout(guides = "collect", ncol = 3) 
p_pcoa_combine

ggsave(filename = "./2.Beta_Diversity/pcoa_month_col1.pdf",
       plot = p_pcoa_combine,
       height = 10,
       width = 10)


#####-----NMDS month Levels-----#####
nmds_out_list <- list()

# nmds
for (month in unique(str_split(colnames(otu_rare), pattern = "_", simplify = T)[,1])) {
  print(month)
  # month = "2019.6"
  otu_rare_sub <- seqtab %>%
    dplyr::select(starts_with(month))
  
  otu_rare_sub <- otu_rare_sub[which(rowSums(otu_rare_sub) != 0), ]
  
  group <- data.frame(
    sample = colnames(otu_rare_sub),
    group = str_remove(string = colnames(otu_rare_sub), pattern = "_\\d+") %>%
      str_remove(string = ., pattern = ".*_")
  )
  
  # 计算bray_curtis距离
  otu.distance <- vegdist(t(otu_rare_sub), method = 'bray')
  # NMDS排序分析 metaMDS 
  df_nmds <- metaMDS(otu.distance, k = 2)
  
  #结果查看——关注stress、points及species三个指标
  summary(df_nmds)

  #适用性检验-基于stress结果
  #应力函数值（<=0.2合理）
  df_nmds_stress <- df_nmds$stress
  df_nmds_stress
  
  #检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
  stressplot(df_nmds)
  
  #提取作图数据
  plot_data <- as.data.frame(df_nmds$points)
  colnames(plot_data) <- c('NMDS1', 'NMDS2')
  
  plot_data$group = str_remove(rownames(plot_data), pattern = "_\\d+")
  plot_data$sample <- rownames(plot_data)
  plot_data$Treatment = str_remove(string = plot_data$group, pattern = ".*_")
  plot_data$Data = str_remove(string = plot_data$group, pattern = "_.*")
  
  p_nmds <- ggplot(plot_data, aes(x=NMDS1,y=NMDS2, color=Treatment))+
    # geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    stat_ellipse(level=0.75,size=1)+
    # ggforce::geom_mark_ellipse(aes(x = NMDS1, y = NMDS2,color = Treatment, fill = Treatment), 
    #                            expand = unit(1, "mm"),
    #                            alpha = 0.5) + 
    geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    # geom_text_repel(aes(x = PCoA1, y = PCoA2, label = sample)) +
    # labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),
    #      y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
    # ggtitle(label = paste("Group:", month, "\n", "Stress = ", round(df_nmds_stress, 3))) + 
    facet_wrap(~Data) + 
    xlim(c(-0.65, 0.65)) +
    ylim(c(-0.55, 0.55)) +
    scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
    scale_fill_manual(values = c("Control"= "#a6bddb", "Warming" = "#fdbb84")) + 
    geom_vline(aes(xintercept=0),linetype="dotted")+
    geom_hline(aes(yintercept=0),linetype="dotted")+
    annotate(geom = "text", x = 0, y = 0.45, label = paste("Stress = ", round(df_nmds_stress, 3)), size=6) +
    theme(panel.background = element_rect(fill='white',colour = 'black'),
          axis.title.x=element_text(colour = 'black',size=20),
          axis.title.y = element_text(colour = 'black',size=20),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(colour = "black", size = 15),
          legend.background = element_roundrect(color = "#808080",linetype = 1),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 25),
          strip.background = element_rect(color = "black")
          )
  
  p_nmds
  
  nmds_out_list[[month]] <- p_nmds
}


nmds_out_list[[1]]
nmds_out_list[[2]]
nmds_out_list[[3]]
nmds_out_list[[4]]
nmds_out_list[[5]]
nmds_out_list[[6]]
nmds_out_list[[7]]
nmds_out_list[[8]]
nmds_out_list[[9]]
nmds_out_list[[10]]
nmds_out_list[[11]]


# 拼图
p_nmds_combine <- nmds_out_list[[1]] + nmds_out_list[[2]] +  nmds_out_list[[3]] + nmds_out_list[[4]] + nmds_out_list[[5]] +
  nmds_out_list[[6]] + nmds_out_list[[7]] + nmds_out_list[[8]] + nmds_out_list[[9]] +
  nmds_out_list[[10]] + nmds_out_list[[11]] + plot_layout(guides = "collect", nrow = 3) 


ggsave(filename = "./2.Beta_Diversity/nmds_month.pdf",
       plot = p_nmds_combine,
       height = 10,
       width = 15)


####----0-10cm ----####
#####----load Data ----#####
# OTU table
seqtab <- read.table(file = "./Input_File/16S_otutab_clean_0_10cm.txt", header = T, row.names = 1, 
                     sep = "\t", check.names = F,  comment.char = "")
seqtab <- seqtab[which(rowSums(seqtab) != 0),]
dim(seqtab)

# taxa file
taxa <- read.table(file = "./Input_File/taxonomy.txt", header = T, row.names = 1)
dim(taxa)

# rare OTU table
otu_rare = as.data.frame(t(rrarefy(t(seqtab), 25000)))

colSums(otu_rare)
colSums(seqtab)

#####-----PCoA month Levels-----#####
pcoa_out_list <- list()
# pcoa
for (month in unique(str_split(colnames(otu_rare), pattern = "_", simplify = T)[,1])) {
  print(month)
  # month = "2019.6"
  otu_rare_sub <- seqtab %>%
    dplyr::select(starts_with(month))
  
  otu_rare_sub <- otu_rare_sub[which(rowSums(otu_rare_sub) != 0), ]
  
  group <- data.frame(
    sample = colnames(otu_rare_sub),
    group = str_remove(string = colnames(otu_rare_sub), pattern = "_\\d+") %>%
      str_remove(string = ., pattern = ".*_")
  )
  
  # PERMANOVA统计检验
  adonis_result <- adonis2(t(otu_rare_sub) ~ group, data = group,permutations = 999,method="bray")
  adonis_result
  adonis_result$R2
  adonis_result$`Pr(>F)`
  
  dune_adonis <- paste0("adonis R2: ",round(adonis_result$R2[1],2), "; P-value: ", adonis_result$`Pr(>F)`[1])
  dune_adonis
  
  distance<-vegdist(t(otu_rare_sub), method='bray')
  
  as.matrix(distance) %>% as.data.frame() %>% write.csv(file = paste0("./2.Beta_Diversity/",month, "beta_bray.csv"))
  
  pcoa<- cmdscale(distance,k=(nrow(t(otu_rare_sub))-1),eig=TRUE)
  
  plot_data<-data.frame({pcoa$point})[1:2]
  
  names(plot_data)[1:2]<-c('PCoA1','PCoA2')
  
  eig=pcoa$eig
  
  plot_data$group = str_remove(rownames(plot_data), pattern = "_\\d+")
  plot_data$sample <- rownames(plot_data)
  plot_data$Treatment = str_remove(string = plot_data$group, pattern = ".*_")
  plot_data$Data = str_remove(string = plot_data$group, pattern = "_.*")
  
  p_pcoa <- ggplot(plot_data, aes(x=PCoA1,y=PCoA2, color=Treatment))+
    # geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    stat_ellipse(level=0.75,size=1)+
    # ggforce::geom_mark_ellipse(aes(x = PCoA1, y = PCoA2,color = Treatment, fill = Treatment), 
    #                            expand = unit(1, "mm"),
    #                            alpha = 0.5) + 
    geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    # geom_text_repel(aes(x = PCoA1, y = PCoA2, label = sample)) +
    labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),
         y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
    # ggtitle(label = paste("Group:", month, "\n", dune_adonis)) + 
    facet_wrap(~Data) + 
    xlim(c(-0.5, 0.7)) +
    ylim(c(-0.35, 0.35)) +
    scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
    scale_fill_manual(values = c("Control"= "#a6bddb", "Warming" = "#fdbb84")) + 
    geom_vline(aes(xintercept=0),linetype="dotted")+
    geom_hline(aes(yintercept=0),linetype="dotted")+
    annotate(geom = "text", x = 0, y = 0.35, label = dune_adonis, size=6) +
    theme(panel.background = element_rect(fill='white',colour = 'black'),
          axis.title.x=element_text(colour = 'black',size=20),
          axis.title.y = element_text(colour = 'black',size=20),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(colour = "black", size = 15),
          legend.background = element_roundrect(color = "#808080",linetype = 1),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 25),
          strip.background = element_rect(color = "black")
    )
  
  p_pcoa
  
  pcoa_out_list[[month]] <- p_pcoa
}

pcoa_out_list[[1]]
pcoa_out_list[[2]]
pcoa_out_list[[3]]
pcoa_out_list[[4]]
pcoa_out_list[[5]]
pcoa_out_list[[6]]
pcoa_out_list[[7]]
pcoa_out_list[[8]]
pcoa_out_list[[9]]
pcoa_out_list[[10]]
pcoa_out_list[[11]]
# pcoa_out_list[[12]]
# pcoa_out_list[[13]]

# 拼图
p_pcoa_combine <- pcoa_out_list[[1]] + pcoa_out_list[[2]] +  pcoa_out_list[[3]] + pcoa_out_list[[4]] + pcoa_out_list[[5]] +
  pcoa_out_list[[6]] + pcoa_out_list[[7]] + pcoa_out_list[[8]] + pcoa_out_list[[9]] +
  pcoa_out_list[[10]] + pcoa_out_list[[11]] + plot_layout(guides = "collect") 
p_pcoa_combine

ggsave(filename = "./2.Beta_Diversity/pcoa_month_0_10cm.pdf",
       plot = p_pcoa_combine,
       height = 13,
       width = 20)

#####-----NMDS month Levels-----#####
nmds_out_list <- list()

# nmds
for (month in unique(str_split(colnames(otu_rare), pattern = "_", simplify = T)[,1])) {
  print(month)
  # month = "2019.6"
  otu_rare_sub <- seqtab %>%
    dplyr::select(starts_with(month))
  
  otu_rare_sub <- otu_rare_sub[which(rowSums(otu_rare_sub) != 0), ]
  
  group <- data.frame(
    sample = colnames(otu_rare_sub),
    group = str_remove(string = colnames(otu_rare_sub), pattern = "_\\d+") %>%
      str_remove(string = ., pattern = ".*_")
  )
  
  # 计算bray_curtis距离
  otu.distance <- vegdist(t(otu_rare_sub), method = 'bray')
  # NMDS排序分析 metaMDS 
  df_nmds <- metaMDS(otu.distance, k = 2)
  
  #结果查看——关注stress、points及species三个指标
  summary(df_nmds)
  
  #适用性检验-基于stress结果
  #应力函数值（<=0.2合理）
  df_nmds_stress <- df_nmds$stress
  df_nmds_stress
  
  #检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
  stressplot(df_nmds)
  
  #提取作图数据
  plot_data <- as.data.frame(df_nmds$points)
  colnames(plot_data) <- c('NMDS1', 'NMDS2')
  
  plot_data$group = str_remove(rownames(plot_data), pattern = "_\\d+")
  plot_data$sample <- rownames(plot_data)
  plot_data$Treatment = str_remove(string = plot_data$group, pattern = ".*_")
  plot_data$Data = str_remove(string = plot_data$group, pattern = "_.*")
  
  p_nmds <- ggplot(plot_data, aes(x=NMDS1,y=NMDS2, color=Treatment))+
    # geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    stat_ellipse(level=0.75,size=1)+
    # ggforce::geom_mark_ellipse(aes(x = NMDS1, y = NMDS2,color = Treatment, fill = Treatment), 
    #                            expand = unit(1, "mm"),
    #                            alpha = 0.5) + 
    geom_point(aes(shape=Treatment), alpha=0.6, size=5)+
    # geom_text_repel(aes(x = PCoA1, y = PCoA2, label = sample)) +
    # labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 4),"%)",sep=""),
    #      y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 4),"%)",sep=""))+
    # ggtitle(label = paste("Group:", month, "\n", "Stress = ", round(df_nmds_stress, 3))) + 
    facet_wrap(~Data) + 
    xlim(c(-0.65, 0.65)) +
    ylim(c(-0.55, 0.55)) +
    scale_color_manual(values = c("Control" = "#74a9cf", "Warming" = "#fc4e2a")) + 
    scale_fill_manual(values = c("Control"= "#a6bddb", "Warming" = "#fdbb84")) + 
    geom_vline(aes(xintercept=0),linetype="dotted")+
    geom_hline(aes(yintercept=0),linetype="dotted")+
    annotate(geom = "text", x = 0, y = 0.45, label = paste("Stress = ", round(df_nmds_stress, 3)), size=6) +
    theme(panel.background = element_rect(fill='white',colour = 'black'),
          axis.title.x=element_text(colour = 'black',size=20),
          axis.title.y = element_text(colour = 'black',size=20),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(colour = "black", size = 15),
          legend.background = element_roundrect(color = "#808080",linetype = 1),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 25),
          strip.background = element_rect(color = "black")
    )
  
  p_nmds
  
  nmds_out_list[[month]] <- p_nmds
}


nmds_out_list[[1]]
nmds_out_list[[2]]
nmds_out_list[[3]]
nmds_out_list[[4]]
nmds_out_list[[5]]
nmds_out_list[[6]]
nmds_out_list[[7]]
nmds_out_list[[8]]
nmds_out_list[[9]]
nmds_out_list[[10]]
nmds_out_list[[11]]


# 拼图
p_nmds_combine <- nmds_out_list[[1]] + nmds_out_list[[2]] +  nmds_out_list[[3]] + nmds_out_list[[4]] + nmds_out_list[[5]] +
  nmds_out_list[[6]] + nmds_out_list[[7]] + nmds_out_list[[8]] + nmds_out_list[[9]] +
  nmds_out_list[[10]] + nmds_out_list[[11]] + plot_layout(guides = "collect") 


ggsave(filename = "./2.Beta_Diversity/nmds_month_0_10cm.pdf",
       plot = p_nmds_combine,
       height = 13,
       width = 20)
