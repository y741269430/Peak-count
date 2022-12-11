# featurecounts

## 1. Read the bed files and annotated in R  
    pka <- function(his){
      peak <- lapply(list.files(paste0('chip-nt-rawdata/', his), "*.bed"), function(x){
        return(readPeakFile(file.path(paste0('chip-nt-rawdata/', his), x)))
      })

      names(peak) <- c('e11.5', 'e12.5', 'e13.5', 'e14.5', 'e15.5')
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 0), 
                             annoDb="org.Mm.eg.db", verbose=FALSE, overlap="all")
      peakAnno_df <- lapply(peakAnnoList, function(x){x <- as.data.frame(x)})

      region_bed <- lapply(peakAnno_df, function(x){
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
      pm_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[grep("Promoter", ignore.case = F, x$annotation), ]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
      gb_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL),]
        x <- x[c(grep("5' UTR", ignore.case = F, x$annotation),
                 grep("Intron", ignore.case = F, x$annotation),
                 grep("Exon", ignore.case = F, x$annotation),
                 grep("Downstream", ignore.case = F, x$annotation),
                 grep("3' UTR", ignore.case = F, x$annotation)), ]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })
      dis_bed <- lapply(peakAnno_df, function(x){
        x <- x[-grep("Rik$", ignore.case = F, x$SYMBOL), ]
        x <- x[grep("Distal Intergenic", ignore.case = F, x$annotation),]
        colnames(x)[6:12] <- c('name', 'score', 'strand2', 'signalValue', 'pValue', 'qValue', 'peak')
        x <- x[, c(1,2,3,23,7,8)]
        return(x)
      })

      save(peakAnno_df, peakAnnoList, region_bed, pm_bed, gb_bed, dis_bed, 
           file = paste0('chip-nt-rawdata/', his, '/chip_Anno_df.RData'))

      for (i in 1:length(pm_bed)) {
        write.table(pm_bed[i],
                    paste(paste0("chip-nt-rawdata/", his, "/pm_saf/"), names(pm_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
      for (i in 1:length(gb_bed)) {
        write.table(gb_bed[i],
                    paste(paste0("chip-nt-rawdata/", his, "/gb_saf/"), names(gb_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
      for (i in 1:length(dis_bed)) {
        write.table(dis_bed[i],
                    paste(paste0("chip-nt-rawdata/", his, "/dis_saf/"), names(dis_bed[i]), "_allpeak.bed", sep = ''),
                    sep = "\t", row.names = F, col.names = F, quote = F)
      }
    }
    pka('H3K9me3')

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
