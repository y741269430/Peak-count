# 使用多种方法对peak进行定量  

## 目录 ####
- 1.加载macs3输出的narrowPeak（已过滤blacklists）  
- 2.基础绘图，peak占比，TSS热图等   
- 3.featurecount
- 3.1 筛选saf文件所需要的列（GeneID, Chr, Start, End and Strand）
- 3.2 输出saf文件格式的bed
- 3.3 进行featurecount定量  
- 4.Multicov 计算count
- 4.1 构建 Multicov 计算count所用的bed文件
- 4.2 将bed文件取交集并且计算count


#### 设置工作路径 ####  
```r
readpath = '/home/jjyang/yjj/ATAC/blacklist_rm/'
path = '/home/jjyang/yjj/ATAC/'
```

## 1.加载macs3输出的narrowPeak（已过滤blacklists）  
```r
peak <- lapply(list.files(readpath, "*.narrowPeak"), function(x){
  return(readPeakFile(file.path(readpath, x)))
})

names(peak) <- str_split_fixed(list.files(readpath, "*.narrowPeak"), '_p', n = 2)[,1]
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#### 注释  
peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), 
                       annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})

peakAnno_df <- lapply(peakAnno_df, function(x){
  colnames(x)[6:12] <- c('name', 'score', 'strand', 'signalValue', 'log10pValue', 'log10qValue', 'summit_peak_start')
  return(x)
})

for (i in 1:length(peakAnno_df)) {
  peakAnno_df[[i]]$Group <- str_split_fixed(peakAnno_df[[i]]$annotation, ' ', n = 2)[,1]
  peakAnno_df[[i]]$Group[which(peakAnno_df[[i]]$Group == '3\'')] = '3UTR'
  peakAnno_df[[i]]$Group[which(peakAnno_df[[i]]$Group == '5\'')] = '5UTR'
}
```

## 2.基础绘图，peak占比，TSS热图等   
```r
plotDistToTSS(peakAnnoList)
ggplot2::ggsave(paste0(path, "results/plotDistToTSS.pdf"),
                height = 5, width = 8, dpi = 300, limitsize = FALSE)

plotAnnoBar(peakAnnoList)
ggplot2::ggsave(paste0(path, "results/plotAnnoBar.pdf"),
                height = 5, width = 8, dpi = 300, limitsize = FALSE)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim =c(-3000, 3000), conf=0.95, resample=500, facet="row")

ggplot2::ggsave(paste0(path, "results/plotAvgProf.pdf"),
                height = 20, width = 8, dpi = 300, limitsize = FALSE)
```

## 3.featurecount   
参考  
- [Nat Immunol 25, 496–511 (2024)](https://www.nature.com/articles/s41590-024-01753-9)
- [Science advances vol. 10,13 (2024)](https://www.science.org/doi/10.1126/sciadv.adk0858)
- [medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.04.02.25325129v1)

### 3.1 筛选saf文件所需要的列（GeneID, Chr, Start, End and Strand）   
储存在本地results, count文件夹内，需要提前创建该文件夹   
```r
region_bed <- lapply(peakAnno_df, function(x){
  x$SYMBOL <- paste0('Peak','_', 1:nrow(x), '_', x$SYMBOL, '_', x$Group)
  x <- x[, c("SYMBOL","seqnames","start","end", "strand")]
  #x$SYMBOL <- paste0('Peak', 1:nrow(x),'_', x$SYMBOL)
  colnames(x) <- c("GeneID","Chr","Start","End","Strand")
  return(x)
})
save(peakAnno_df, peakAnnoList, region_bed, file = paste0(path, 'results/Anno_df.RData'))
```
### 3.2 输出saf文件格式的bed    
```r
for (i in 1:length(region_bed)) {
  write.table(x = region_bed[[i]],
              file = paste0(path, 'count/', names(region_bed)[i], '_all.saf'),
              sep = "\t", row.names = FALSE, col.names = colnames(region_bed[[i]]), quote = FALSE)
}
```
### 3.3 进行featurecount定量    
制作bam文件和saf文件的路径，要求这两个数据框是一一对应的   
```r
bamPath <- paste0(path, 'bam/')
bamNames <- dir(bamPath, pattern = "last.bam$") 
bamPath <- sapply(bamNames, function(x){paste0(bamPath,x)})     
bamPath <- data.frame(bamPath); rownames(bamPath) <- str_split_fixed(bamNames, "_", n = 2)[,1]

safPath <- paste0(path, 'count/')
safNames <- dir(safPath, pattern = "_all.saf$") 
safPath <- sapply(safNames, function(x){paste0(safPath,x)})   
safPath <- data.frame(safPath)

bamPath
safPath
```
featureCounts计算  
具体参考  
- https://subread.sourceforge.net/SubreadUsersGuide.pdf

> `countMultiMappingReads` A multi-mapping read is a read that maps to more than one location in the reference genome. There are multiple options for counting such reads. Users can specify the ‘-M’ option (set `countMultiMappingReads` to TRUE in R) to fully count every alignment reported for a multimapping read (each alignment carries 1 count), or specify both ‘-M’ and ‘–fraction’ options (set both `countMultiMappingReads` and `fraction` to TRUE in R) to count each alignment fractionally (each alignment carries 1/x count where x is the total number of alignments reported for the read), or do not count such reads at all (this is the default behavior in SourceForge Subread package; In R, you need to set `countMultiMappingReads` to FALSE).

> `allowMultiOverlap` A multi-overlapping read is a read that overlaps more than one meta-feature when counting reads at meta-feature level or overlaps more than one feature when counting reads at feature level. The decision of whether or not to counting these reads is often determined by the experiment type. We recommend that reads or fragments overlapping more than one gene are not counted for RNA-seq experiments, because any single fragment must originate from only one of the target genes but the identity of the true target gene cannot be confidently determined. On the other hand, we recommend that multi-overlapping reads or fragments are counted for ChIP-seq experiments because for example epigenetic modifications inferred from these reads may regulate the biological functions of all their overlapping genes.

> `allowMultiOverlap` By default, featureCounts does not count multi-overlapping reads. Users can specify the ‘-O’ option (set `allowMultiOverlap` to TRUE in R) to fully count them for each overlapping metafeature/feature (each overlapping meta-feature/feature receives a count of 1 from a read), or specify both ‘-O’ and ‘–fraction’ options (set both `allowMultiOverlap` and `fraction` to TRUE in R) to assign a fractional count to each overlapping meta-feature/feature (each overlapping meta-feature/feature receives a count of 1/y from a read where y is the total number of meta-features/features overlapping with the read).

```r
peak_counts <- list()
for (i in 1:nrow(safPath)) {
  peak_counts[[i]] <- Rsubread::featureCounts(files = as.character(bamPath[i,]),
                                              allowMultiOverlap = T, 
                                              fraction = T,
                                              countMultiMappingReads = T,
                                              countChimericFragments = F,
                                              nthreads = 8,
                                              isPairedEnd = T,
                                              annot.ext = as.character(safPath[i,]))
  rm(i)
}
```

输出count list 并拼接样本矩阵（这里可能要先取基因交集再决定用哪些样本）

```r
count_list <- lapply(peak_counts, function(x){
  x <- x$counts
  x <- data.frame(x)
  x$SYMBOL <- rownames(x)
  return(x)
})

names(count_list) <- c(  )      # 重命名
raw_counts <- Reduce(merge, count_list)
colnames(raw_counts)[-1] <- c(  )      # 重命名

save(count_list, file = paste0(path, 'results/count_list.RData'))
write.csv(raw_counts, paste0(path, 'results/raw_counts.csv'), row.names = F)
```
取交集的code如下  
```r
data_ls <- list('con1' = count_list[[1]]$SYMBOL, 
                'con2' = count_list[[2]]$SYMBOL, 
                'con3' = count_list[[3]]$SYMBOL)

ggvenn(data_ls, 
       fill_color = c("#f0d0d9", "#fbe5c8", "#70c1b3", "#ff5733", "#00bfff"),
       stroke_color = 'white',
       text_size = 10)
```

## 4.Multicov 计算count   
参考  
- [Nat Commun 14, 7712 (2023)](https://www.nature.com/articles/s41467-023-43427-4#Sec11)
- [multicov](https://bedtools.readthedocs.io/en/latest/content/tools/multicov.html)    
- [Bedtools-Multicov---reads计数](https://www.jianshu.com/p/641c2c2cfd41)   

### 4.1 构建 Multicov 计算count所用的bed文件
```r
count_bed <- lapply(peakAnno_df, function(x){
  x$unique_name <- paste('Peak', 1:nrow(x), x$SYMBOL, x$Group, sep = '_')
  x <- x[, c("seqnames","start","end", "unique_name")]
  return(x)
})

for (i in 1:length(count_bed)) {
  write.table(x = count_bed[[i]],
              file = paste0(path, 'count/', names(count_bed)[i], '_Multicov.bed'),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
```
### 4.2 将bed文件取交集并且计算count   
```bash
bedtools intersect -a CON_1_Multicov.bed -b CON_2_Multicov.bed > CON.bed
bedtools intersect -a Tre_1_Multicov.bed -b Tre_2_Multicov.bed > Tre.bed
bedtools intersect -a CON.bed -b Tre.bed > Multicov_input.bed

bedtools multicov -bams CON_1.last.bam CON_2.last.bam Tre_1.last.bam Tre_2.last.bam -bed Multicov_input.bed > Multicov_count.csv &
```

---
## 5.构建meme-chip所需的bed文件  
采用启动子区peak中心位置左右扩展100bp作为motif预测的region，因为两个核小体之间大概50-80bp，而测序片段约200-300bp，长度越长预测越不精确。  
参考  
- [MEME-ChIP](https://meme-suite.org/meme/doc/meme-chip.html?man_type=web)  
- [MEMEsuite](https://github.com/y741269430/MEMEsuite)   
> homer软件里面的[peakMotifs](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)里面是这么写的：  
> `Region Size` ("-size <#>", "-size <#>,<#>", "-size given", default: 200)
> The size of the region used for motif finding is important.  If analyzing ChIP-Seq peaks from a transcription factor, Chuck would recommend 50 bp for establishing the primary motif bound by a given transcription factor and 200 bp for finding both primary and "co-enriched" motifs for a transcription factor.  When looking at histone marked regions, 500-1000 bp is probably a good idea (i.e. H3K4me or H3/H4 acetylated regions).  In theory, HOMER can work with very large regions (i.e. 10kb), but with the larger the regions comes more sequence and longer execution time.  These regions will be based off the center of the peaks.  If you prefer an offset, you can specify "-size -300,100" to search a region of size 400 that is centered 100 bp upstream of the peak center (useful if doing motif finding on putative TSS regions).  If you have variable length regions, use the option "-size given" and HOMER will use the exact regions that were used as input.  

```r
pm_peak500_region <- lapply(peakAnno_df, function(x){
  x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
  x$start100 <- x$start + x$summit_peak_start - 250
  x$start200 <- x$start + x$summit_peak_start + 250
  x$unique_name <- paste(x$SYMBOL, c(x$start200 - 250) ,sep = '_')
  x <- x[, c("seqnames","start100","start200", "unique_name")]
  return(x)
})
```

> 这里我为什么又要用500bp作区间呢，原因是后续我们做meme分析的时候，在peak intersect这一步，得到的peak的长度是不一致的，导致meme出现bug。
> 于是取500bp，然后将交集得到的区间，补全到所需的300bp就行了。  
> 参考 [2.peak取重叠区域做meme分析](https://github.com/y741269430/MEMEsuite?tab=readme-ov-file#2peak%E5%8F%96%E9%87%8D%E5%8F%A0%E5%8C%BA%E5%9F%9F%E5%81%9Ameme%E5%88%86%E6%9E%90)

```r
# 预先创建peak500文件夹，该输出文件将使用bedtools进行构建fasta
dir.create(paste0(path, 'peak500/'))

for (i in 1:length(pm_peak500_region)) {
  write.table(x = pm_peak500_region[[i]],
              file = paste0(path, 'peak500/', names(pm_peak500_region)[i], '_equal_p.bed'),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

```
---







