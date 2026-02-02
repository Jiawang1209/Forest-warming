DESeq2_function <- function(data,
                            sample1, sample2, sample3,sample4, sample5, sample6,
                            sample7, sample8, sample9,sample10, sample11, sample12,
                            outname, path){
  
  # get OTU Table
  seatab <- get(data)
  df <- seqtab[, c(sample1, sample2, sample3,sample4, sample5, sample6,
                   sample7, sample8, sample9,sample10, sample11, sample12)]
  
  # get metadata
  metadata <- data.frame(
    sample = c(sample1, sample2, sample3,sample4, sample5, sample6,
               sample7, sample8, sample9,sample10, sample11, sample12),
    group = rep(c("control","treatment"), each = 6)
  )
  
  # create dds object
  dds <-DESeqDataSetFromMatrix(countData=round(df,0),
                               colData=metadata,
                               design=~group)
  
  # check dds
  nrow(dds)
  rownames(dds)
  # filter dds
  #dds <- dds[rowSums(counts(dds))>1,]
  nrow(dds)
  # DESeq program
  dds <- DESeq(dds)
  
  #LogFC adjust
  ## set contrast parameter
  contrast=c("group", "treatment", "control")
  
  # get DEG results
  dd1 <- results(dds, contrast=contrast)
  dd1Ordered <- dd1[order(dd1$pvalue),] # 按照P值排序
  dd1 <- as.data.frame(dd1Ordered)
  head(dd1)
  dd1 <- na.omit(dd1)
  dim(dd1)
  # get DEG results
  res <- dd1 %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    dplyr::mutate(change = case_when(
      log2FoldChange >= 1 & pvalue < 0.05 ~ "Up",
      log2FoldChange <= -1 & pvalue < 0.05 ~ "Down",
      abs(log2FoldChange) < 1 | pvalue > 0.05 ~ "Normal"
    ))
  table(res$change)
  nrDEG <- res
  
  # output
  write.csv(nrDEG, file = paste0(path,'/', outname, ".DESeq2.DEG.csv"))
  
}