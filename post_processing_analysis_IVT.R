library('compEpiTools')
library('GenomicRanges')
library('ggplot2')

ELIGOS_results <- function(path_input, path_output) {
  
  # create a grange object with all the nucleotides that have been analysed by ELIGOS, so those with a coverage of at least 20x
  analysed_sites <- function(dir) {
    list_gr <- list()
    
    # for each file in the directory produce a GRanges object containing the chromosome, start and end position, strand
    # and name of the comparison
    for (file in dir) {
      t <- read.table(file, sep= '\t', header=TRUE)
      
      gr <- GRanges(seqnames = Rle(t$chrom), 
                    ranges = IRanges(start = t$start_loc, end=t$end_loc), 
                    strand = Rle(t$strand),
                    sample = strsplit(strsplit(file, split = '/')[[1]][7],'\\.')[[1]][1])
      
      list_gr <- c(list_gr, gr)
    }
    
    return(list_gr)
  }

  analysed_sites_comparison_1 <- analysed_sites(path_input)[[1]]
  analysed_sites_reduced_comparison_1 <- reduce(analysed_sites_comparison_1)
  print(paste(unique(analysed_sites_comparison_1$sample), 'has', as.character(sum(width(unique(analysed_sites_reduced_comparison_1)))), 'analysed nucleotides, distributed on',
              as.character(length(unique(as.vector(seqnames(analysed_sites_reduced_comparison_1))))), 'transcripts'))
 
  analysed_sites_comparison_2 <- analysed_sites(path_input)[[2]]
  analysed_sites_reduced_comparison_2 <- reduce(analysed_sites_comparison_2)
  print(paste(unique(analysed_sites_comparison_2$sample), 'has', as.character(sum(width(unique(analysed_sites_reduced_comparison_2)))), 'analysed nucleotides, distributed on', 
              as.character(length(unique(as.vector(seqnames(analysed_sites_reduced_comparison_2))))), 'transcripts')) 
 
  smaller_grange <- analysed_sites_reduced_comparison_1
  bigger_grange <- analysed_sites_reduced_comparison_2
  if (sum(width(unique(analysed_sites_reduced_comparison_1))) > sum(width(unique(analysed_sites_reduced_comparison_2)))) {
    smaller_grange <- analysed_sites_reduced_comparison_2
    bigger_grange <- analysed_sites_reduced_comparison_1
  }
  
  # identify the nucleotides analysed in both comparisons
  analysed_sites_both_comparisons <- findOverlaps(smaller_grange,bigger_grange,type = 'any')
  analysed_sites_both_comparisons <- smaller_grange[unique(queryHits(analysed_sites_both_comparisons))]
  analysed_sites_both_comparisons <- reduce(analysed_sites_both_comparisons)
  print(paste(as.character(sum(width(unique(analysed_sites_both_comparisons)))), 
              'nucleotides were analysed both in WT MinION vs IVT PromethION comparison and WT PromethION vs IVT PromethION one. These nucleotides are on',
              as.character(length(unique(as.vector(seqnames(analysed_sites_both_comparisons))))), 'transcripts'))
  
  # create a grange object with the nucleotides that fulfill the conditions to be considered modifications
  hits <- function(dir,expansion) {
    list_gr <- list()
    
    # for each file in the directory produce a GRanges object containing the chromosome, start and end position, strand, 
    # pvalue, pvalue adjusted, odds ratio and name of the comparison
    for (file in dir) {
      t <- read.table(file, sep= '\t', header=TRUE)
      
      gr <- GRanges(seqnames = Rle(t$chrom), 
                    ranges = IRanges(start = t$start_loc-expansion, end=t$end_loc+expansion), 
                    strand = Rle(t$strand),
                    pval = t$pval,
                    pvalAdj = t$adjPval, 
                    oddR = t$oddR,
                    sample = strsplit(strsplit(file, split = '/')[[1]][7],'\\.')[[1]][1])
      
      # filter for the modifications with a pvalue <= 0.05, an adjusted pvalue lower than 0.0001 and an oddR higher than 2.5
      gr <- gr[gr$pval <= 0.05 & gr$pvalAdj <= 1e-04 & gr$oddR >= 2.5]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  
  # each range is expanded of 10 nucleotides centered around the corresponding hit
  hits_10nt_condition_1 <- hits(path_input,5)[[1]]
  hits_10nt_reduced_condition_1 <- reduce(hits_10nt_condition_1)
  print(paste(unique(hits_10nt_condition_1$sample), 'has', as.character(length(unique(hits_10nt_reduced_condition_1))), 'hits, distributed on',
              as.character(length(unique(as.vector(seqnames(hits_10nt_reduced_condition_1))))), 'transcripts'))
  
  hits_10nt_condition_2 <- hits(path_input,5)[[2]]
  hits_10nt_reduced_condition_2 <- reduce(hits_10nt_condition_2)
  print(paste(unique(hits_10nt_condition_2$sample), 'has', as.character(length(unique(hits_10nt_reduced_condition_2))), 'hits, distributed on',
              as.character(length(unique(as.vector(seqnames(hits_10nt_reduced_condition_2))))), 'transcripts'))
  
  # identify the single nucleotides that fulfill the parameters without range expansion
  hits_condition_1 <- hits(path_input,0)[[1]]
  hits_reduced_condition_1 <- reduce(hits_condition_1)
  hits_condition_2 <- hits(path_input,0)[[2]]
  hits_reduced_condition_2 <- reduce(hits_condition_2)
  
  # keep only the modified nucleotides that are in sites that have been analysed in both comparisons
  hits_reduced_condition_1_common <- findOverlaps(hits_reduced_condition_1, analysed_sites_both_comparisons, type = 'any')
  hits_reduced_condition_1_common <- hits_reduced_condition_1[unique(queryHits(hits_reduced_condition_1_common))]
  hits_reduced_condition_1_common <- reduce(hits_reduced_condition_1_common)
  
  hits_reduced_condition_2_common <- findOverlaps(hits_reduced_condition_2, analysed_sites_both_comparisons, type = 'any')
  hits_reduced_condition_2_common <- hits_reduced_condition_2[unique(queryHits(hits_reduced_condition_2_common))]
  hits_reduced_condition_2_common <- reduce(hits_reduced_condition_2_common)
  
  # expand the hits in the analysed sites by both the comparisons creating a range of 10 nucleotides centered around the corresponding hit
  start(hits_reduced_condition_1_common) <- start(hits_reduced_condition_1_common) -5
  end(hits_reduced_condition_1_common) <- end(hits_reduced_condition_1_common) +5
  hits_reduced_condition_1_common <- reduce(hits_reduced_condition_1_common)
  print(paste(unique(hits_condition_1$sample), 'has', as.character(length(unique(hits_reduced_condition_1_common))), 
              'hits in sites analysed by both the comparisons, distributed on', as.character(length(unique(as.vector(seqnames(hits_reduced_condition_1_common))))),
              'transcripts'))
  
  start(hits_reduced_condition_2_common) <- start(hits_reduced_condition_2_common) -5
  end(hits_reduced_condition_2_common) <- end(hits_reduced_condition_2_common) +5
  hits_reduced_condition_2_common <- reduce(hits_reduced_condition_2_common)
  print(paste(unique(hits_condition_2$sample), 'has', as.character(length(unique(hits_reduced_condition_2_common))), 
              'hits in sites analysed by both the comparisons, distributed on', as.character(length(unique(as.vector(seqnames(hits_reduced_condition_2_common))))),
              'transcripts'))
  
  smaller_grange <- hits_reduced_condition_1_common
  bigger_grange <- hits_reduced_condition_2_common
  if (length(unique(hits_reduced_condition_1_common)) > length(unique(hits_reduced_condition_2_common))) {
    smaller_grange <- hits_reduced_condition_2_common
    bigger_grange <- hits_reduced_condition_1_common
  }
  
  intersection_hits <- findOverlaps(smaller_grange,bigger_grange, type = 'any')
  intersection_hits <- smaller_grange[unique(queryHits(intersection_hits))]
  intersection_hits <- reduce(intersection_hits)
  print(paste('WT MinION vs IVT PromethION comparison and WT PromethION vs IVT PromethION one have', as.character(length(unique(intersection_hits))), 
        'hits in common among those that are in sites analysable in both comparisons. These hits are on', as.character(length(unique(as.vector(seqnames(intersection_hits))))), 'transcripts'))
  
  intersection_hits <- list(hits_reduced_condition_1_common, hits_reduced_condition_2_common)
  names(intersection_hits) <- c(unique(hits_condition_1$sample), unique(hits_condition_2$sample))
  
  pdf(file = path_output, width = 4, height = 4)
  overlapOfGRanges(intersection_hits,plot = TRUE)
  dev.off()
  
  hits_common_granges <- list(hits_reduced_condition_1_common,hits_reduced_condition_2_common)
  names(hits_common_granges) <- c(unique(hits_condition_1$sample),unique(hits_condition_2$sample))
  
  hits_common_values <- c(length(unique(hits_reduced_condition_1_common)),length(unique(hits_reduced_condition_2_common)))
  names(hits_common_values) <- c(unique(hits_condition_1$sample),unique(hits_condition_2$sample))
  
  hits_common <- list(hits_common_granges, hits_common_values)
  return(hits_common)
}

Nanocompore_results <- function(path_input, path_output, path_bed_3utr_top_tx) {
  
  # create a grange object with all the nucleotides that have been analysed by Nanocompore, so those with a coverage of at least 30x
  analysed_sites <- function(dir) {
    list_gr <- list()
    
    for (file in dir) {
      # for each file in the directory produce a GRanges object containing the transcript name, start and end position, strand and name of the comparison
      t <- read.table(file, header=TRUE)
      
      gr <- GRanges(seqnames = t$ref_id, 
                    ranges = IRanges(start = t$pos, end=t$pos), 
                    strand = Rle(t$strand), 
                    sample = strsplit(strsplit(file, split = '/')[[1]][7],'\\.')[[1]][1])
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)
  }
  
  analysed_sites_comparison_1 <- analysed_sites(path_input)[[1]]
  analysed_sites_reduced_comparison_1 <- reduce(analysed_sites_comparison_1)
  
  analysed_sites_comparison_2 <- analysed_sites(path_input)[[2]]
  analysed_sites_reduced_comparison_2 <- reduce(analysed_sites_comparison_2)
  
  # load the coordinates of the 3' UTR of the selected transcripts 
  coordinates_3utr_tx <- read.table(path_bed_3utr_top_tx)
  
  coordinates_3utr_tx <- GRanges(seqnames = coordinates_3utr_tx$V1,
                                   ranges = IRanges(start = coordinates_3utr_tx$V2, end=coordinates_3utr_tx$V3), 
                                   strand = Rle(coordinates_3utr_tx$V6))
  
  # keep only the analysed sites that are in the 3' UTR of the selected transcripts
  analysed_sites_reduced_comparison_1_3UTR_tx <- findOverlaps(analysed_sites_reduced_comparison_1,coordinates_3utr_tx, type = 'any')
  analysed_sites_reduced_comparison_1_3UTR_tx <- analysed_sites_reduced_comparison_1[unique(queryHits(analysed_sites_reduced_comparison_1_3UTR_tx))]
  analysed_sites_reduced_comparison_1_3UTR_tx <- reduce(analysed_sites_reduced_comparison_1_3UTR_tx)
  print(paste(unique(analysed_sites_comparison_1$sample), 'has', as.character(sum(width(unique(analysed_sites_reduced_comparison_1_3UTR_tx)))), 'analysed nucleotides, distributed on',
              as.character(length(unique(as.vector(seqnames(analysed_sites_reduced_comparison_1_3UTR_tx))))), 'transcripts'))
  
  analysed_sites_reduced_comparison_2_3UTR_tx <- findOverlaps(analysed_sites_reduced_comparison_2,coordinates_3utr_tx, type = 'any')
  analysed_sites_reduced_comparison_2_3UTR_tx <- analysed_sites_reduced_comparison_2[unique(queryHits(analysed_sites_reduced_comparison_2_3UTR_tx))]
  analysed_sites_reduced_comparison_2_3UTR_tx <- reduce(analysed_sites_reduced_comparison_2_3UTR_tx)
  print(paste(unique(analysed_sites_comparison_2$sample), 'has', as.character(sum(width(unique(analysed_sites_reduced_comparison_2_3UTR_tx)))), 'analysed nucleotides, distributed on',
              as.character(length(unique(as.vector(seqnames(analysed_sites_reduced_comparison_2_3UTR_tx))))), 'transcripts'))
  
  smaller_grange <- analysed_sites_reduced_comparison_1_3UTR_tx
  bigger_grange <- analysed_sites_reduced_comparison_2_3UTR_tx
  if (sum(width(unique(analysed_sites_reduced_comparison_1_3UTR_tx))) > sum(width(unique(analysed_sites_reduced_comparison_2_3UTR_tx)))) {
    smaller_grange <- analysed_sites_reduced_comparison_2_3UTR_tx
    bigger_grange <- analysed_sites_reduced_comparison_1_3UTR_tx
  }
  
  # identify the nucleotides analysed in both comparisons
  analysed_sites_both_comparisons <- findOverlaps(smaller_grange,bigger_grange,type = 'any')
  analysed_sites_both_comparisons <- smaller_grange[unique(queryHits(analysed_sites_both_comparisons))]
  analysed_sites_both_comparisons <- reduce(analysed_sites_both_comparisons)
  print(paste(as.character(sum(width(unique(analysed_sites_both_comparisons)))), 
              'nucleotides were analysed both in WT MinION vs IVT PromethION comparison and WT PromethION vs IVT PromethION one. These nucleotides are on',
              as.character(length(unique(as.vector(seqnames(analysed_sites_both_comparisons))))), 'transcripts'))
  
  # create a grange object with the nucleotides that fulfill the conditions to be considered modifications
  hits <- function(dir,expansion) {
    
    list_gr <- list()
    
    # for each file in the directory produce a GRanges object containing the transcript name, start and end position, strand,
    # pvalue, odds ratio and name of the comparison
    for (file in dir) {
      t <- read.table(file, header=TRUE)
      t <- t[which(t$Logit_LOR !='NC'),]
      
      gr <- GRanges(seqnames = t$ref_id,
                    ranges = IRanges(start = t$pos-expansion, end=t$pos+expansion), 
                    strand = Rle(t$strand), 
                    p_value = t$GMM_logit_pvalue,
                    logit_LOR = as.numeric(t$Logit_LOR),
                    sample = strsplit(strsplit(file, split = '/')[[1]][7],'\\.')[[1]][1])
      
      # filter for the modifications with a pvalue lower than 0.01
      gr <- gr[which(gr$p_value <= 0.01)]
      
      # filter for the modifications with an absolute value of logit_LOR higher than 0.5
      gr <- gr[which(abs(gr$logit_LOR) >= 0.5)]
      
      list_gr <- c(list_gr, gr)
    }
    return(list_gr)

  }
  
  # each range is expanded of 10 nucleotides centered around the corresponding hit
  hits_10nt_condition_1 <- hits(path_input,5)[[1]]
  hits_10nt_reduced_condition_1 <- reduce(hits_10nt_condition_1)
  
  hits_10nt_condition_2 <- hits(path_input,5)[[2]]
  hits_10nt_reduced_condition_2 <- reduce(hits_10nt_condition_2)
  
  # keep only the expanded hits that are in the 3' UTR of the selected transcripts
  hits_10nt_condition_1_3UTR_tx <- findOverlaps(hits_10nt_reduced_condition_1,coordinates_3utr_tx, type = 'any')
  hits_10nt_condition_1_3UTR_tx <- hits_10nt_reduced_condition_1[unique(queryHits(hits_10nt_condition_1_3UTR_tx))]
  hits_10nt_condition_1_3UTR_tx <- reduce(hits_10nt_condition_1_3UTR_tx)
  print(paste(unique(hits_10nt_condition_1$sample), 'has', as.character(length(unique(hits_10nt_condition_1_3UTR_tx))), 'hits, distributed on',
              as.character(length(unique(as.vector(seqnames(hits_10nt_condition_1_3UTR_tx))))), 'transcripts'))
  
  hits_10nt_condition_2_3UTR_tx <- findOverlaps(hits_10nt_reduced_condition_2,coordinates_3utr_tx, type = 'any')
  hits_10nt_condition_2_3UTR_tx <- hits_10nt_reduced_condition_2[unique(queryHits(hits_10nt_condition_2_3UTR_tx))]
  hits_10nt_condition_2_3UTR_tx <- reduce(hits_10nt_condition_2_3UTR_tx)
  print(paste(unique(hits_10nt_condition_2$sample), 'has', as.character(length(unique(hits_10nt_condition_2_3UTR_tx))), 'hits, distributed on',
              as.character(length(unique(as.vector(seqnames(hits_10nt_condition_2_3UTR_tx))))), 'transcripts'))
  
  # identify the single nucleotides that fulfill the parameters without range expansion
  hits_condition_1 <- hits(path_input,0)[[1]]
  hits_reduced_condition_1 <- reduce(hits_condition_1)
  hits_condition_2 <- hits(path_input,0)[[2]]
  hits_reduced_condition_2 <- reduce(hits_condition_2)
  
  # keep only the single base hits that are in the 3' UTR of the selected transcripts
  hits_reduced_condition_1_3UTR_tx <- findOverlaps(hits_reduced_condition_1,coordinates_3utr_tx, type = 'any')
  hits_reduced_condition_1_3UTR_tx <- hits_reduced_condition_1[unique(queryHits(hits_reduced_condition_1_3UTR_tx))]
  hits_reduced_condition_1_3UTR_tx <- reduce(hits_reduced_condition_1_3UTR_tx)
  
  hits_reduced_condition_2_3UTR_tx <- findOverlaps(hits_reduced_condition_2,coordinates_3utr_tx, type = 'any')
  hits_reduced_condition_2_3UTR_tx <- hits_reduced_condition_2[unique(queryHits(hits_reduced_condition_2_3UTR_tx))]
  hits_reduced_condition_2_3UTR_tx <- reduce(hits_reduced_condition_2_3UTR_tx)
  
  # keep only the modified nucleotides that are in sites that have been analysed in both comparisons
  hits_reduced_condition_1_3UTR_tx_common <- findOverlaps(hits_reduced_condition_1_3UTR_tx, analysed_sites_both_comparisons, type = 'any')
  hits_reduced_condition_1_3UTR_tx_common <- hits_reduced_condition_1_3UTR_tx[unique(queryHits(hits_reduced_condition_1_3UTR_tx_common))]
  hits_reduced_condition_1_3UTR_tx_common <- reduce(hits_reduced_condition_1_3UTR_tx_common)
  
  hits_reduced_condition_2_3UTR_tx_common <- findOverlaps(hits_reduced_condition_2_3UTR_tx, analysed_sites_both_comparisons, type = 'any')
  hits_reduced_condition_2_3UTR_tx_common <- hits_reduced_condition_2_3UTR_tx[unique(queryHits(hits_reduced_condition_2_3UTR_tx_common))]
  hits_reduced_condition_2_3UTR_tx_common <- reduce(hits_reduced_condition_2_3UTR_tx_common)
  
  # expand the hits in the analysed sites by both the comparisons creating a range of 10 nucleotides centered around the corresponding hit
  start(hits_reduced_condition_1_3UTR_tx_common) <- start(hits_reduced_condition_1_3UTR_tx_common) -5
  end(hits_reduced_condition_1_3UTR_tx_common) <- end(hits_reduced_condition_1_3UTR_tx_common) +5
  hits_reduced_condition_1_3UTR_tx_common <- reduce(hits_reduced_condition_1_3UTR_tx_common)
  print(paste(unique(hits_condition_1$sample), 'has', as.character(length(unique(hits_reduced_condition_1_3UTR_tx_common))), 
              'hits in sites analysed by both the comparisons, distributed on', as.character(length(unique(as.vector(seqnames(hits_reduced_condition_1_3UTR_tx_common))))),
              'transcripts'))
  
  start(hits_reduced_condition_2_3UTR_tx_common) <- start(hits_reduced_condition_2_3UTR_tx_common) -5
  end(hits_reduced_condition_2_3UTR_tx_common) <- end(hits_reduced_condition_2_3UTR_tx_common) +5
  hits_reduced_condition_2_3UTR_tx_common <- reduce(hits_reduced_condition_2_3UTR_tx_common)
  print(paste(unique(hits_condition_2$sample), 'has', as.character(length(unique(hits_reduced_condition_2_3UTR_tx_common))), 
              'hits in sites analysed by both the comparisons, distributed on', as.character(length(unique(as.vector(seqnames(hits_reduced_condition_2_3UTR_tx_common))))),
              'transcripts'))
  
  smaller_grange <- hits_reduced_condition_1_3UTR_tx_common
  bigger_grange <- hits_reduced_condition_2_3UTR_tx_common
  if (length(unique(hits_reduced_condition_1_3UTR_tx_common)) > length(unique(hits_reduced_condition_2_3UTR_tx_common))) {
    smaller_grange <- hits_reduced_condition_2_3UTR_tx_common
    bigger_grange <- hits_reduced_condition_1_3UTR_tx_common
  }
  
  intersection_hits <- findOverlaps(smaller_grange,bigger_grange, type = 'any')
  intersection_hits <- smaller_grange[unique(queryHits(intersection_hits))]
  intersection_hits <- reduce(intersection_hits)
  print(paste('WT MinION vs IVT PromethION comparison and WT PromethION vs IVT PromethION one have', as.character(length(unique(intersection_hits))), 
              'hits in common among those that are in sites analysable in both comparisons. These hits are on', as.character(length(unique(as.vector(seqnames(intersection_hits))))), 'transcripts'))
  
  intersection_hits <- list(hits_reduced_condition_1_3UTR_tx_common, hits_reduced_condition_2_3UTR_tx_common)
  names(intersection_hits) <- c(unique(hits_condition_1$sample), unique(hits_condition_2$sample))
  
  pdf(file = path_output, width = 4, height = 4)
  overlapOfGRanges(intersection_hits,plot = TRUE)
  dev.off()
  
  hits_common_granges <- list(hits_reduced_condition_1_3UTR_tx_common,hits_reduced_condition_2_3UTR_tx_common)
  names(hits_common_granges) <- c(unique(hits_condition_1$sample),unique(hits_condition_2$sample))
  
  hits_common_values <- c(length(unique(hits_reduced_condition_1_3UTR_tx_common)),length(unique(hits_reduced_condition_2_3UTR_tx_common)))
  names(hits_common_values) <- c(unique(hits_condition_1$sample),unique(hits_condition_2$sample))
  
  hits_common <- list(hits_common_granges, hits_common_values)
  return(hits_common)
}

# path to the directory containing the result of ELIGOS analysis for WT MinION vs IVT PromethION comparison and  
# WT PromethION vs IVT PromethION one
path_ELIGOS_results <- list.files(path = path_ELIGOS_output, pattern = "txt", full.names = TRUE)

hits_ELIGOS <- ELIGOS_results(path_ELIGOS_results, path_output)

hits_ELIGOS_granges <- hits_ELIGOS[[1]]
hits_ELIGOS_values <- hits_ELIGOS[[2]]

# path to the directory containing the result of Nanocompore analysis for WT MinION vs IVT PromethION comparison and  
# WT PromethION vs IVT PromethION one
path_Nanocompore_results <- list.files(path = path_Nanocompore_output, pattern = "tsv", full.names = TRUE)

hits_Nanocompore <- Nanocompore_results(path_Nanocompore_results, path_output, path_bed_3utr_top_tx)

hits_Nanocompore_granges <- hits_Nanocompore[[1]]
hits_Nanocompore_values <- hits_Nanocompore[[2]]

# create a dataset with the number of hits in sites analysed by both the comparisons
tool <- c(rep("ELIGOS" , 2) , rep("Nanocompore" , 2))
Comparison <- c(rep(c(names(hits_ELIGOS_values)[1],names(hits_ELIGOS_values)[2]),2))
number_of_hits <- c(hits_ELIGOS_values[1],
                    hits_ELIGOS_values[2],
                    hits_Nanocompore_values[names(hits_ELIGOS_values)[1]],
                    hits_Nanocompore_values[names(hits_ELIGOS_values)[2]])
data <- data.frame(tool,number_of_hits,Comparison)

# plot the data in a grouped histogram
pdf(file = '', width = 5, height = 2)
ggplot(data, aes(fill=Comparison, y=number_of_hits, x=tool)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  coord_flip()
dev.off()

# overlap between the filtered hits in common analyzed sites in the same comparison between Nanocompore and ELIGOS
overlap_Nano_ELIGOS_condition_1 <- list(hits_ELIGOS_granges[[1]],hits_Nanocompore_granges[names(hits_ELIGOS_granges)[1]][[1]])
names(overlap_Nano_ELIGOS_condition_1) <- c('ELIGOS', 'Nanocompore')

pdf(file = '', width = 4, height = 4)
overlapOfGRanges(overlap_Nano_ELIGOS_condition_1,plot = TRUE)
dev.off()

overlap_Nano_ELIGOS_condition_2 <- list(hits_ELIGOS_granges[[2]],hits_Nanocompore_granges[names(hits_ELIGOS_granges)[2]][[1]])
names(overlap_Nano_ELIGOS_condition_2) <- c('ELIGOS', 'Nanocompore')

pdf(file = '', width = 4, height = 4)
overlapOfGRanges(overlap_Nano_ELIGOS_condition_2,plot = TRUE)
dev.off()


