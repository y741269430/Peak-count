myGrouplist <- function(x){
  num <- stringr::str_split_fixed(colnames(x), "_", n = 2)[,1]
  num <- factor(num, levels = unique(num)); num
  num <- data.frame(table(num)); num
  
  grouplist <- rep(num[1:nrow(num),1], num[1:nrow(num),2]); grouplist
}
myNormal <- function(x, grouplist){
  colData <- data.frame(row.names = colnames(x), grouplist)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                        colData = colData,
                                        design = ~ grouplist)
  keep <- rowSums(DESeq2::counts(dds)) >= 50
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)  
  nor <- DESeq2::counts(dds, normalized = T)
  nor2 <- as.data.frame(nor)
}
myPCA <- function(x, grouplist, label_option = "ind") {
  factoextra::fviz_pca_ind(prcomp(t(x),scale. = T), 
                           repel = TRUE, pointsize = 5, 
                           palette = "jco", 
                           label = ifelse(label_option == "none", "none", "ind"), 
                           mean.point = F, 
                           col.ind = grouplist) + 
    ggplot2::coord_fixed(1) + 
    ggplot2::labs(title = element_blank())+ 
    ggplot2::theme_bw()+
    ggplot2::theme(text = element_text(family = "serif"),
                   plot.title = element_text(hjust = 0.5),
                   axis.title = element_text(size = 15),
                   axis.text = element_text(size = 15),
                   legend.title = element_blank(), 
                   legend.text = element_text(face = "bold", size = 15),
                   legend.position = "bottom") + 
    xlim(-250, 250) + ylim(-250, 250)
}
myDESeq2 <- function(x, ct, tr, id = 'ENSEMBL'){
  # raw counts
  localgeneid <- read.csv("downloads/Motif_database/rmdup_id.csv")[,-5]
  
  condition <- factor(c(rep("ctrl", ct), rep("treatment", tr)), 
                      levels = c("ctrl", "treatment"))
  
  colData <- data.frame(row.names = colnames(x), condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                        colData = colData,
                                        design = ~condition)
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq2::DESeq(dds)    
  nor <- DESeq2::counts(dds, normalized = T)
  
  DEG_DESeq2 <- merge(data.frame(results(dds, contrast = c("condition", "treatment", "ctrl"))), nor, by = "row.names")
  DEG_DESeq2 <- na.omit(DEG_DESeq2)
  colnames(DEG_DESeq2)[1] <- id
  colnames(DEG_DESeq2)[which(colnames(DEG_DESeq2) == 'pvalue')] = 'pvalue'
  colnames(DEG_DESeq2)[which(colnames(DEG_DESeq2) == 'log2FoldChange')] = 'log2FC'
  DEG <- merge(localgeneid, DEG_DESeq2, id)
  
  return(DEG)
}
myDotplot <- function(x, showCategory = 10, title = NULL, color_scheme = 5, angle = 45, term_key = ''){
  
  if ("compareClusterResult" %in% slotNames(x)) {
    x@compareClusterResult <- x@compareClusterResult[grep(term_key, x@compareClusterResult$Description), ]
  } else if ("result" %in% slotNames(x)) {
    x@result <- x@result[grep(term_key, x@result$Description), ]
  }
  
  # 颜色方案列表
  color_schemes <- list(
    c('red', 'blue'),  # 方案1
    c('#4182A4', '#b3eebe'),  # 方案2
    c('#b22832', '#f6f0bc'),  # 方案3
    c('#5B7E74', '#e5e4da'),  # 方案4
    c('#dc6767', '#fbd6d4', '#b9cbe4','#357cb7')    # 方案4
    # 添加更多方案...
  )
  
  # 获取指定的颜色方案，如果方案不存在，则使用默认方案1
  if (color_scheme > length(color_schemes) || color_scheme < 1) {
    warning("Invalid color_scheme option. Using default color scheme (Option 1).")
    color_scheme <- 1
  }
  colors <- color_schemes[[color_scheme]]
  
  # 输入一个数据框，进行点图绘制，默认选择p.adjust前10的通路
  dotplot(x, showCategory = showCategory) +
    ggplot2::theme(
      #text = element_text(family = "serif"),
      plot.title = element_text(size = 20, colour = "red"),
      axis.text.x = element_text(angle = angle, hjust = 1),
      axis.title.x = element_blank()) +
    ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=70)) +
    ggtitle(title) +
    scale_color_gradientn(colours = colors, guide = guide_colorbar(reverse = F))
} # 改
extractGroups <- function(x, group){
  x@compareClusterResult <- x@compareClusterResult[x@compareClusterResult$Cluster %in%
                                                     names(x@geneClusters[group]),]
  return(x)
} # Goterm2dotplot
geneID2SYMBOL2 <- function(x, data = NULL){
  
  # 替换 geneID 列中的 ENTREZID 为 SYMBOL
  updated_rows <- x
  
  # 如果 data 为空，则读取默认文件
  if (is.null(data)) { data <- read.csv('downloads/Motif_database/rmdup_id.csv') }
  
  entrez_to_symbol <- setNames(data$SYMBOL, data$ENTREZID)
  
  updated_rows$geneID <- sapply(x$geneID, function(gene_ids) {
    ids <- unlist(strsplit(gene_ids, "/"))
    symbols <- entrez_to_symbol[ids]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]  # 保留未匹配到的值
    paste(symbols, collapse = "/")
  })
  
  x <- updated_rows
  
  return(x)
}
my_enrichKEGG <- function(gene, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH", universe, minGSSize = 10, maxGSSize = 500, 
                          qvalueCutoff = 0.2, use_internal_data = FALSE) {
  species <- clusterProfiler:::organismMapper(organism)
  if (use_internal_data) {KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db(species)}
  else {
    astk_dir <- file.path("R/")
    kegg_data_name <- paste(organism, lubridate::today(),"kegg","RData", sep = ".")
    kegg_data_path <- file.path(astk_dir, kegg_data_name)
    if (file.exists(kegg_data_path)){
      print("load local catch")
      load(kegg_data_path)
    }else {
      print("online downloading... ")
      KEGG_DATA <- clusterProfiler:::prepare_KEGG(organism,"KEGG", "kegg")
      save(KEGG_DATA, file = kegg_data_path) 
    }; print(KEGG_DATA)
    
  }
  
  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff, 
                                             pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                                             maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = KEGG_DATA)
  if (is.null(res)) 
    return(res)
  res@ontology <- "KEGG"
  res@organism <- species
  res@keytype <- keyType
  return(res)
}