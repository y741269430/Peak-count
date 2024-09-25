# featurecounts

## 1. Read the bed files and annotated in R   
```r  
readpath = '/home/jjyang/yjj/ATAC/blacklist_rm/'
path = '/home/jjyang/yjj/ATAC/'

peak <- lapply(list.files(readpath, "*.narrowPeak"), function(x){
  return(readPeakFile(file.path(readpath, x)))
})

names(peak) <- str_split_fixed(list.files(readpath, "*.narrowPeak"), '_p', n = 2)[,1]
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


names(peak) <- c('CTRL_2', 'CTRL_4', 'CTRL_5', 'CTRL_6', 'CTRL_7', 'CTRL_8', 'CTRL_9', 'CFA_31', 'CFA_33', 'CFA_34', 'CFA_1', 'CFA_2')


peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), 
                       annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})

peakAnno_df <- lapply(peakAnno_df, function(x){
  colnames(x)[6:12] <- c('name', 'score', 'strand', 'signalValue', 'log10pValue', 'log10qValue', 'summit_peak_start')
  return(x)
})


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

#############

region_bed <- lapply(peakAnno_df, function(x){
  x <- x[, c("SYMBOL","seqnames","start","end", "strand")]
  x$SYMBOL <- paste0('Peak', 1:nrow(x),'_', x$SYMBOL)
  colnames(x) <- c("GeneID","Chr","Start","End","Strand")
  return(x)
})

pm_bed <- lapply(peakAnno_df, function(x){
  x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
  x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
  x <- x[, c("SYMBOL","seqnames","start","end", "strand")]
  x$SYMBOL <- paste0('Peak', 1:nrow(x),'_', x$SYMBOL)
  colnames(x) <- c("GeneID","Chr","Start","End","Strand")
  return(x)
})

gb_bed <- lapply(peakAnno_df, function(x){
  x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
  x <- x[grep("Promoter|Distal Intergenic", ignore.case = F, x$annotation), ]
  x <- x[, c("SYMBOL","seqnames","start","end", "strand")]
  x$SYMBOL <- paste0('Peak', 1:nrow(x),'_', x$SYMBOL)
  colnames(x) <- c("GeneID","Chr","Start","End","Strand")
  return(x)
})


save(peakAnno_df, peakAnnoList, region_bed, pm_bed, gb_bed, 
     file = paste0(path, 'count/Anno_df.RData'))


for (i in 1:length(pm_bed)) {
  write.table(x = pm_bed[[i]],
              file = paste0(path, 'count/', names(pm_bed)[i], '_pm.bed'),
              sep = "\t", row.names = FALSE, col.names = colnames(pm_bed[[i]]), quote = FALSE)
}

for (i in 1:length(gb_bed)) {
  write.table(x = gb_bed[[i]],
              file = paste0(path, 'count/', names(gb_bed)[i], '_gb.bed'),
              sep = "\t", row.names = FALSE, col.names = colnames(gb_bed[[i]]), quote = FALSE)
}
```





## 2.BED to saf in linux  

The BED files can be used to convert to saf files for featurecount, it also can be used to fimo analysis.

    vim bed2saf.sh

    #!/bin/bash
    ## make saf (bedtools) for featurecount ##

    path=./pm_saf

    cat filenames | while read i; 
    do
    nohup bedtools sort -i $path/${i}_allpeak.bed > $path/${i}.rm.bed && bedtools merge -c 4,6 -o first -i $path/${i}.rm.bed |awk 'BEGIN{print "GeneID" "\t"  "Chr" "\t" "Start" "\t" "End" "\t" "Strand"}{print $4"\t"$1"\t"strtonum($2)"\t"strtonum($3)"\t"$5}' > $path/${i}.saf && rm $path/${i}.rm.bed -rf &
    done

or  

    #!/bin/bash
    ## make saf (bedtools) for featurecount ##

    cat filenames | while read i; 
    do
    nohup bedtools sort -i ./pm_saf/${i}_allpeak.bed > ./pm_saf/${i}.rm.bed && bedtools merge -c 4,6 -o first -i ./pm_saf/${i}.rm.bed |awk 'BEGIN{print "GeneID" "\t"  "Chr" "\t" "Start" "\t" "End" "\t" "Strand"}{print $4"\t"$1"\t"strtonum($2)"\t"strtonum($3)"\t"$5}' > ./pm_saf/${i}.saf && rm ./pm_saf/${i}.rm.bed -rf &

    nohup bedtools sort -i ./gb_saf/${i}_allpeak.bed > ./gb_saf/${i}.rm.bed && bedtools merge -c 4,6 -o first -i ./gb_saf/${i}.rm.bed |awk 'BEGIN{print "GeneID" "\t"  "Chr" "\t" "Start" "\t" "End" "\t" "Strand"}{print $4"\t"$1"\t"strtonum($2)"\t"strtonum($3)"\t"$5}' > ./gb_saf/${i}.saf && rm ./gb_saf/${i}.rm.bed -rf &

    nohup bedtools sort -i ./dis_saf/${i}_allpeak.bed > ./dis_saf/${i}.rm.bed && bedtools merge -c 4,6 -o first -i ./dis_saf/${i}.rm.bed |awk 'BEGIN{print "GeneID" "\t"  "Chr" "\t" "Start" "\t" "End" "\t" "Strand"}{print $4"\t"$1"\t"strtonum($2)"\t"strtonum($3)"\t"$5}' > ./dis_saf/${i}.saf && rm ./dis_saf/${i}.rm.bed -rf &

    done

## 3.Perform featurecounts in R  

    pkc <- function(name){
      bamPath <- paste0("chip-nt-rawdata/bam/", his)
      bamNames <- dir(bamPath, pattern = ".bam$") 
      bamPath <- sapply(bamNames, function(x){paste(bamPath, x, sep='/')})   
      bamPath <- data.frame(bamPath); rownames(bamPath) <- str_split_fixed(bamNames, "_", n = 2)[,1]
      #
      safPath <- paste0("chip-nt-rawdata/", his, "/" , name, "_saf")
      safNames <- dir(safPath, pattern = ".saf$") 
      safPath <- sapply(safNames, function(x){paste(safPath, x, sep='/')})   
      safPath <- data.frame(safPath) 

      peak_counts <- list()
      for (i in 1:5) {
        peak_counts[[i]] <- Rsubread::featureCounts(files = as.character(bamPath[c(2*i-1, 2*i),]),
                                                    allowMultiOverlap = T,
                                                    countMultiMappingReads = F,
                                                    countChimericFragments = F,
                                                    nthreads = 8,
                                                    isPairedEnd = T,
                                                    annot.ext = as.character(safPath[i,]))
        rm(i)
      }

      count_list <- lapply(peak_counts, function(x){
        x <- x$counts
        x <- data.frame(x)
        x$SYMBOL <- rownames(x)
        return(x)
      })

      names(count_list) <- c('e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5')
      raw_counts <- Reduce(merge, count_list)

      colnames(raw_counts)[-1] <- c('nt_e11.5_1', 'nt_e11.5_2', 
                                    'nt_e12.5_1', 'nt_e12.5_2', 
                                    'nt_e13.5_1', 'nt_e13.5_2', 
                                    'nt_e14.5_1', 'nt_e14.5_2', 
                                    'nt_e15.5_1', 'nt_e15.5_2')

      save(count_list, file = paste0("chip-nt-rawdata/", his, "/", his, "_", name, "_count_list.RData"))
      write.csv(raw_counts, paste0("chip-nt-rawdata/", his, "/", his, "_", name, "_counts.csv"), row.names = F)
    }
    his <- 'H3K9me3'
    pkc('pm')
    pkc('gb')
    pkc('dis')
