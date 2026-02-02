STAMP_function <- function(data,
                           sample1, sample2, sample3,sample4, sample5, sample6,
                           sample7, sample8, sample9,sample10, sample11, sample12,
                           groupname, path){

  seqtab <- get(data)
  
  # get Group Table
  group_table <- data.frame(
    sample = c(sample1, sample2, sample3,sample4, sample5, sample6,
               sample7, sample8, sample9,sample10, sample11, sample12),
    group = rep(c("control","treatment"), each = 6)
  )
  
  # trans to relative abundance
  seqtab_rare <- seqtab / colSums(seqtab)
  
  # get sub seqtab
  seqtab_rare_sub <- seqtab_rare[, c(sample1, sample2, sample3,sample4, sample5, sample6,
                                     sample7, sample8, sample9,sample10, sample11, sample12)]
  
  # trans to 100% format
  seqtab_rare_sub_100 <- seqtab_rare_sub * 100
  dim(seqtab_rare_sub_100)
  
  # # filter > 0.5%
  seqtab_rare_sub_100_filter <- seqtab_rare_sub_100 %>% dplyr::filter(apply(seqtab_rare_sub_100, 1, mean) > 0.01)
  dim(seqtab_rare_sub_100_filter)
  seqtab_rare_sub_100_filter.t <- t(seqtab_rare_sub_100_filter)
  
  data1 = data.frame(seqtab_rare_sub_100_filter.t,
                     group_table$group)
  colnames(data1)[length(colnames(data1))] <- "Group"
  
  #  t-test select signif Genus
  glimpse(data1)
  
  dfg <- data1 %>%
    dplyr::select_if(is.numeric) %>%
    purrr::map_df(~ broom::tidy(t.test(. ~ Group, data = data1)), .id = "var")
  
  dfg$p.value = p.adjust(dfg$p.value, "none")
  
  dfg <-  dfg %>% dplyr::filter(p.value < 0.05)
  dfg
  
  dfg.t.test <- dfg
  
  # plot
  dfg.t.test
  
  abun.bar = data1[,c(dfg.t.test$var, "Group")] %>% 
    tidyr::pivot_longer(cols = !Group,
                        names_to = "Variable",
                        values_to = "Value") %>%
    dplyr::group_by(Variable, Group) %>%
    dplyr::summarise(Mean = mean(Value))
  
  abun.bar
  
  # 右侧散点图
  diff.mean <- dfg.t.test[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0, levels(factor(data1$Group))[1], levels(factor(data1$Group))[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]
  diff.mean
  
  # 左侧条形图
  # set order
  abun.bar$Variable  =  factor(abun.bar$Variable,levels = rev(diff.mean$var))
  abun.bar$Group <- str_remove(string = abun.bar$Group, pattern = ".*_")
  
  p1 <- ggplot(abun.bar, aes(Variable, Mean, fill = Group)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("control" = "#74a9cf", "treatment" = "#fc4e2a")) + 
    coord_flip() + 
    xlab("") +#X轴标签
    ylab("Mean proportion (%)") + 
    labs(title=paste0("Group: ", groupname))  + 
    geom_bar(stat = "identity", position = "dodge") +
    theme(panel.background = element_rect(fill = 'transparent'),#主题设定，
          panel.grid = element_blank(),#背景格子为空
          axis.ticks.length = unit(0.4,"lines"), #坐标轴刻度线长，正数向外
          axis.ticks = element_line(color='black'),#坐标轴刻度线颜色
          axis.line = element_line(colour = "black"),#坐标轴线颜色
          axis.title.x=element_text(colour='black', size=12,face = "bold"),#X轴文本设置
          axis.text=element_text(colour='black',size=10,face = "bold"),#坐标轴文本设定
          legend.title=element_blank(),#图例标题为空
          legend.position = "top",#图例位置，左下
          legend.direction = "horizontal",#图例方向，水平排列
          legend.key.width = unit(0.8,"cm"),#图例方块宽
          legend.key.height = unit(0.5,"cm"),
          legend.background = element_roundrect(color = "#808080",linetype = 1),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5))#图例方块高
  
  p1
  
  for (i in seq(1, nrow(diff.mean)-1, 2)) {
    p1 <- p1 + 
      annotate(geom = "rect", xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = "#bdbdbd", alpha = 0.25)
  }
  
  p1 <- p1 + geom_bar(stat = "identity", position = "dodge", color = "black")
  
  p1
  
  # 右侧散点图
  
  diff.mean$var = factor(diff.mean$var,levels = levels(abun.bar$Variable))
  diff.mean$p.value = as.numeric(diff.mean$p.value)
  diff.mean$p.value = round(diff.mean$p.value,3)#保留3位小数
  diff.mean$p.value = as.character(diff.mean$p.value)
  diff.mean$Group2 <- str_remove(string = diff.mean$Group, pattern = ".*_")
  
  p2 <- ggplot(diff.mean, aes(var, estimate, fill = Group2)) + 
    geom_point(shape = 21, size = 3) + 
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")  + 
    theme(panel.background = element_rect(fill = 'transparent'),#主题设定，
          panel.grid = element_blank(),#背景格子为空
          axis.ticks.length = unit(0.4,"lines"), #坐标轴刻度线长，正数向外
          axis.ticks = element_line(color='black'),#坐标轴刻度线颜色
          axis.line = element_line(colour = "black"),#坐标轴线颜色
          axis.title.x=element_text(colour='black', size=12,face = "bold"),#X轴文本设置
          axis.text=element_text(colour='black',size=10,face = "bold"),#坐标轴文本设定
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title=element_blank(),#图例标题为空
          legend.position = "none",#图例位置，左下
          legend.direction = "horizontal",#图例方向，水平排列
          legend.key.width = unit(0.8,"cm"),#图例方块宽
          legend.key.height = unit(0.5,"cm"),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5))
  
  p2
  
  for (i in seq(1, nrow(diff.mean)-1, 2)) {
    p2 <- p2 + 
      annotate(geom = "rect", xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf, fill = "#bdbdbd", alpha = 0.25)
  }
  
  # p2
  
  p2 <- p2 + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(0.8), width = 0.25, linewidth = 0.5) +#误差线
    geom_point(shape = 21,size = 3) +#散点图参数
    scale_fill_manual(values = c("control" = "#74a9cf", "treatment" = "#fc4e2a")) +#点颜色
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')#加虚线
  
  p2
  
  # 最左侧 p value
  p3 <- ggplot(diff.mean) + 
    geom_text(aes(x = 0, y = var), label = diff.mean$p.value) + 
    labs(title="P-Value")  + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5))
  
  p3
  
  p <- p1 + p2 + p3 + plot_layout(widths = c(6,4,2))
  
  p
  
  ggsave(filename = paste0("./5.Biomaker/" ,groupname, "_STAMP_Genus_Level.pdf"),
         plot = p,
         height = 15,
         width = 18)
  
  diff.mean
  
  abun.bar %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = "Group", values_from = "Mean") %>%
    left_join(diff.mean, by = c("Variable" = "var")) %>%
    write.csv(file = paste0("./5.Biomaker/" ,groupname, "_STAMP_Genus_Level_stat_out.csv"))
}
