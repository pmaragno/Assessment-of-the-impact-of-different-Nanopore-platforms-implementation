args = commandArgs(trailingOnly=TRUE)

suppressMessages(library('GenomicFeatures'))
suppressMessages(library('GenomicAlignments'))
suppressMessages(library("Rsamtools"))

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

txdb <- makeTxDbFromGFF(gtf_file)
three_UTR_tx <- threeUTRsByTranscript(txdb)
tx <- transcriptsBy(txdb, by = "gene")

# create a data frame with the tx_id and corresponding tx_name considering all the possible transcripts
tx_id_name <- data.frame(tx_id = unlist(tx)$tx_id, tx_name = unlist(tx)$tx_name)
rownames(tx_id_name) <- tx_id_name$tx_id
tx_id_name$tx_id <- NULL

# create a data frame with tx_name, length of the transcript and strand considering all the possible transcripts
tx_3UTR <- data.frame(tx_name= unlist(tx)$tx_name, end=width(unlist(tx)), strand=as.vector(strand(unlist(tx))))
rownames(tx_3UTR) <- tx_3UTR$tx_name
tx_3UTR$tx_name <- NULL

# create a data frame with the tx_id, the length of the 3' UTR and the tx_name
length_3UTR <- data.frame(tx_id=names(unlist(three_UTR_tx)), width=width(unlist(three_UTR_tx)))
length_3UTR <- cbind(length_3UTR, tx_name= tx_id_name[as.vector(length_3UTR$tx_id),])

# compute the sum of the length of all the 3' UTR regions of each transcript
sum_3UTR_by_tx <- lapply(split(length_3UTR$width, length_3UTR$tx_name), sum)

# compute the start of the 3' UTR region as the length of the whole transcript minus the length of its 3' UTR region. 
# The end of the 3' UTR region is the length of the whole transcript 
tx_3UTR <- cbind(tx_3UTR, start = NA)
tx_3UTR[names(sum_3UTR_by_tx),3] <- tx_3UTR[names(sum_3UTR_by_tx),1] - unlist(sum_3UTR_by_tx)
tx_3UTR <- tx_3UTR[-which(tx_3UTR[,1]==tx_3UTR[,3]),]

# save a bed file with the coordinates of the 3' UTR of each transcript
bed_3utr_tx <- cbind(rownames(tx_3UTR), tx_3UTR$start, tx_3UTR$end, '.', as.vector(tx_3UTR$strand))
colnames(bed_3utr_tx) <- c("chr", "start", "end", "score", "strand")
write.table(bed_3utr_tx, path_bed_3utr_tx, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# function to compute the number of sites with coverage at least 30x in a given region, required in input:
# bam_file with the mapping of the reads on the transcriptome
# regions_file is a file in bed format with the coordinates of the 3' UTR of each transcript
# output_dir is the path to the output directory
Calculate_num_sites_coverage_higher_30 <- function(bam_file, regions_file, output_dir) {
  # import the regions_file as a GRanges object
  regions_df <- read.table(regions_file, quote = "")
  regions <- GRanges(seqnames = regions_df$V1,
                     ranges = IRanges(start = regions_df$V2,
                                      end = regions_df$V3),
                     strand = regions_df$V5)
  
  pu_par <- PileupParam(max_depth=1000000, min_base_quality=0, min_mapq=0,
                        min_nucleotide_depth=0, min_minor_allele_depth=0,
                        distinguish_strands=FALSE, distinguish_nucleotides=FALSE,
                        ignore_query_Ns=TRUE, include_deletions=FALSE, include_insertions=FALSE,
                        left_bins=NULL, query_bins=NULL, cycle_bins=NULL)
  
  # import bam file
  pileup_tot <- pileup(file = bam_file, index = paste0(bam_file, ".bai"), pileupParam = pu_par)
  
  # this command must be executed only if the transcript names in the bam file contain additional information, 
  # like the gene name and the Ensembl gene id, and you want to remove these additional data taking only the 
  # Ensembl transcript id
  pileup_tot$seqnames <- gsub(x=as.vector(pileup_tot$seqnames), pattern='\\|.*', replacement = '')
  
  pileup_GR <- GRanges(seqnames = pileup_tot$seqnames,
                       ranges = IRanges(start = pileup_tot$pos, end = pileup_tot$pos))
  mcols(pileup_GR) <- pileup_tot$count
  colnames(mcols(pileup_GR)) <- "coverage"
  
  hits <- findOverlaps(query = pileup_GR, subject = regions, ignore.strand = TRUE)
  coverage_split <- split(pileup_GR[queryHits(hits)], subjectHits(hits))
  regions_names <- unname(unlist(lapply(split(as.character(seqnames(regions[subjectHits(hits)])), subjectHits(hits)), '[[', 1)))
  names(coverage_split) <- regions_names
  regions_names_all <- as.vector(seqnames(regions))
  
  num_sites_higher_30X_tmp <- lapply(coverage_split, function(x) length(x$coverage[x$coverage >=30]))
  num_sites_higher_30X <- rep(0, length(regions_names_all))
  names(num_sites_higher_30X) <- regions_names_all
  # coverage is only computed for sequences available in the BAM file
  num_sites_higher_30X[regions_names] <- num_sites_higher_30X_tmp
  names(num_sites_higher_30X) <- gsub(pattern = "_$", x = gsub(pattern = ":|-", x = names(num_sites_higher_30X), replacement = "_"), replacement = "-")
  num_sites_higher_30X <- unlist(num_sites_higher_30X)
  return(num_sites_higher_30X)
}

regions_file <- path_bed_3utr_tx

# compute the number of sites with coverage at least 30x on the 3' UTR of each transcript in IVT PromethION
num_sites_higher_30X_IVTprom <- Calculate_num_sites_coverage_higher_30(path_to_IVTprom, regions_file, output_dir)
summary(num_sites_higher_30X_IVTprom)

# compute the number of sites with coverage at least 30x on the 3' UTR of each the transcript in WT MinION
num_sites_higher_30X_WTmin <- Calculate_num_sites_coverage_higher_30(path_to_WTmin, regions_file, output_dir)
summary(num_sites_higher_30X_WTmin)

# compute the number of sites with coverage at least 30x on the 3' UTR of each the transcript in WT PromethION
num_sites_higher_30X_WTprom <- Calculate_num_sites_coverage_higher_30(path_to_WTprom, regions_file, output_dir)
summary(num_sites_higher_30X_WTprom)

# create a data frame that for each transcript reports the number of sites with coverage at least 30x in the different conditions
num_sites_higher_30X_WTmin_order_IVTprom <- num_sites_higher_30X_WTmin[names(num_sites_higher_30X_IVTprom)]
num_sites_higher_30X_WTprom_order_IVTprom <- num_sites_higher_30X_WTprom[names(num_sites_higher_30X_IVTprom)]
num_sites_higher_30X <- data.frame(num_sites_higher_30X_IVTprom, num_sites_higher_30X_WTmin_order_IVTprom, num_sites_higher_30X_WTprom_order_IVTprom)
colnames(num_sites_higher_30X) <- c('number_sites_higher_30X_IVTprom','number_sites_higher_30X_WTmin', 'number_sites_higher_30X_WTprom')

# compute for each transcript the minimum number of sites with coverage at least 30x between the different conditions and 
# sort them in decreasing order
min_number_sites_higher_30X <- apply(num_sites_higher_30X,1,min)
ind <- sort(min_number_sites_higher_30X, decreasing = TRUE, index.return=TRUE)$ix
min_number_sites_higher_30X_ordered <- (min_number_sites_higher_30X[ind])

# select the top 50 transcripts that have the higher minimum number of sites with coverage at least 30x between the different conditions
selected_transcripts <- names(min_number_sites_higher_30X_ordered)[1:50]

# save a bed file with the coordinates of the 3' UTR of the 50 selected transcripts 
# since the annotation gtf file I used contains only the Ensembl transcript id while the transcriptome fasta file contains as transcript name also 
# other information, like the gene name and the Ensembl gene id, convert the short version in the more complete one
bam_file <- readGAlignments(path_to_IVTprom, use.names = TRUE)
tx_complete <- unique(as.vector(seqnames(bam_file)))
names(tx_complete) <- gsub(x=tx_complete, pattern='\\|.*', replacement = '')

selected_transcripts_complete <- unname(tx_complete[selected_transcripts])

bed_3utr_top50_tx <- cbind(selected_transcripts_complete, tx_3UTR[selected_transcripts,]$start, tx_3UTR[selected_transcripts,]$end, selected_transcripts_complete, '.', as.vector(tx_3UTR[selected_transcripts,]$strand))
colnames(bed_3utr_top50_tx) <- c("chr", "start", "end","id","score","strand")
write.table(bed_3utr_top50_tx,path_bed_3utr_top50_tx,quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
