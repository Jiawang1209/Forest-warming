LefSe_function <- function(data,
                           sample1, sample2, sample3,sample4, sample5, sample6,
                           sample7, sample8, sample9,sample10, sample11, sample12,
                           groupname, path){

  seqtab <- get(data)
  
  # get sub seqtab
  seqtab_sub <- seqtab[, c(sample1, sample2, sample3,sample4, sample5, sample6,
                           sample7, sample8, sample9,sample10, sample11, sample12)]
  
  # get Group Table
  group_table <- data.frame(
    sample = c(sample1, sample2, sample3,sample4, sample5, sample6,
               sample7, sample8, sample9,sample10, sample11, sample12),
    group = rep(c("control","treatment"), each = 6)
  )
  rownames(group_table) <- group_table$sample
  group_table
  
  # tax_table
  tax_table <- taxa[rownames(taxa) %in% rownames(seqtab_sub), ]
  dim(tax_table)
  tax_table %<>% tidy_taxonomy
  
  
  # create microtable object
  dataset <- microtable$new(sample_table = group_table,
                            otu_table = seqtab_sub, 
                            tax_table = tax_table)
  
  dataset
  
  # LEfSe Analysis
  lefse <- trans_diff$new(dataset = dataset, 
                          method = "lefse", 
                          group = "group", 
                          # alpha = 0.01, 
                          taxa_level = "Species",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")
  # 查看分析结果
  head(lefse$res_diff)
  
  write.csv(lefse$res_diff,
            file = paste0(path, "/", groupname, "_Lefse_Species_levels.table.csv"),
            quote = F)
  
  lefse$plot_diff_bar(use_number = 1:30, 
                      width = 0.8, 
                      group_order = c("control", "treatment"),
                      color_values = c("#74a9cf", "#fc4e2a"))
  
  ggsave(filename = paste0(path, "/", groupname, "_Lefse_Species_levels.pdf"),
         height = 8,
         width = 12)
  
  # LEfSe Analysis
  lefse <- trans_diff$new(dataset = dataset, 
                          method = "lefse", 
                          group = "group", 
                          # alpha = 0.01, 
                          taxa_level = "all",
                          lefse_subgroup = NULL,
                          p_adjust_method = "none")
  
  write.csv(lefse$res_diff,
            file = paste0(path, "/", groupname, "_Lefse_all_levels.table.csv"),
            quote = F)
  
  lefse$plot_diff_cladogram(use_taxa_num = 200, 
                            use_feature_num = 50, 
                            clade_label_level = 5,
                            group_order = c("control", "treatment"),
                            color = c("#74a9cf", "#fc4e2a"))
  
  ggsave(filename = paste0(path, "/", groupname,"_Lefse_all_levels.pdf"),
         height = 12,
         width = 18)
  
}
