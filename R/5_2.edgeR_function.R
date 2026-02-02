edgeR_function <- function(data,
                           sample1, sample2, sample3,sample4, sample5, sample6,
                           sample7, sample8, sample9,sample10, sample11, sample12,
                           outname, path){
  
  seatab <- get(data)
  # get OTU Table
  df <- seqtab[, c(sample7, sample8, sample9,sample10, sample11, sample12,
                   sample1, sample2, sample3, sample4, sample5, sample6)]
  
  # get group
  group = factor(rep(c("treatment","control"), each = 6), levels = c("treatment", "control"), ordered = T)
  
  # DEGList object
  d <- DGEList(counts=df,group=group)
  
  # data control
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  
  # Normal
  d <- calcNormFactors(d)
  d$samples
  
  #构建分组矩阵，进行差异表达分析
  dge=d
  design <- model.matrix(~0+group)
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(group)
  
  # 构建分组矩阵
  design
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast=c(1,-1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG =nrDEG 
  
  edgeR_DEG <- edgeR_DEG %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    dplyr::mutate(change = case_when(
      logFC >= 1 & FDR < 0.05 ~ "Up",
      logFC <= -1 & FDR < 0.05 ~ "Down",
      abs(logFC) < 1 | FDR > 0.05 ~ "Normal"
    ))
  
  table(edgeR_DEG$change)
  
  write.csv(edgeR_DEG, file = paste0(path,'/', outname, ".EdgeR.DEG.csv"))
}