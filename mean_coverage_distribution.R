args = commandArgs(trailingOnly=TRUE)

suppressMessages(library('GenomicFeatures'))
suppressMessages(library('GenomicAlignments'))
suppressMessages(library("Rsamtools"))

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

# function that returns the index of the sites along the 3' UTR of each transcript with coverage at least 30x
Calculate_index_sites_coverage_higher_30 <- function(bam_file, regions_file, output_dir) {
  #import the regions_file as a GRanges object
  regions_df <- read.table(regions_file, quote = "")
  regions <- GRanges(seqnames = regions_df$V1,
                     ranges = IRanges(start = regions_df$V2,
                                      end = regions_df$V3),
                     strand = regions_df$V6)
  
  pu_par <- PileupParam(max_depth=1000000, min_base_quality=0, min_mapq=0,
                        min_nucleotide_depth=0, min_minor_allele_depth=0,
                        distinguish_strands=FALSE, distinguish_nucleotides=FALSE,
                        ignore_query_Ns=TRUE, include_deletions=FALSE, include_insertions=FALSE,
                        left_bins=NULL, query_bins=NULL, cycle_bins=NULL)
  
  pileup_tot <- pileup(file = bam_file, index = paste0(bam_file, ".bai"), pileupParam = pu_par)
    
  pileup_GR <- GRanges(seqnames = pileup_tot$seqnames,
                       ranges = IRanges(start = pileup_tot$pos, end = pileup_tot$pos))
  mcols(pileup_GR) <- pileup_tot$count
  colnames(mcols(pileup_GR)) <- "coverage"
  
  hits <- findOverlaps(query = pileup_GR, subject = regions, ignore.strand = TRUE)
  coverage_split <- split(pileup_GR[queryHits(hits)], subjectHits(hits))
  regions_names <- unname(unlist(lapply(split(as.character(seqnames(regions[subjectHits(hits)])), subjectHits(hits)), '[[', 1)))
  names(coverage_split) <- regions_names
  regions_names_all <- as.vector(seqnames(regions))
  
  index_analysed_sites <- lapply(coverage_split, function(x) {which(x$coverage>=30)}) 
  return(index_analysed_sites)
}

regions_file <- path_bed_3utr_top_tx

# for each transcript return the index of the sites with at least 30 reads mapping on them
index_sites_cover_higher30_IVTprom <- Calculate_index_sites_coverage_higher_30(path_to_IVTprom, regions_file, output_dir)

index_sites_coverage_higher30_WTmin <- Calculate_index_sites_coverage_higher_30(path_to_WTmin, regions_file, output_dir)

index_sites_coverage_higher30_WTprom <- Calculate_index_sites_coverage_higher_30(path_to_WTprom, regions_file, output_dir)

# extract only the transcripts with analyzed sites in all the conditions
tx_IVTprom <- names(index_sites_cover_higher30_IVTprom)
tx_WTmin <- names(index_sites_coverage_higher30_WTmin)
tx_WTprom <- names(index_sites_coverage_higher30_WTprom)

tx_common <- intersect(tx_IVTprom,tx_WTmin)
tx_common <- intersect(tx_common,tx_WTprom)

# produce a list that for each transcript reports only the sites analyzed in all the conditions 
common_sites <- lapply(as.list(tx_common), function(x) {
  sites_IVTprom <- unname(index_sites_cover_higher30_IVTprom[x])[[1]]
  sites_WTmin <- unname(index_sites_coverage_higher30_WTmin[x])[[1]]
  sites_WTprom <- unname(index_sites_coverage_higher30_WTprom[x])[[1]]
  common_sites <- intersect(sites_IVTprom,sites_WTmin)
  common_sites <- intersect(common_sites,sites_WTprom)
  if (length(common_sites) == 0) {
    common_sites <- NA
  }
  return(common_sites)
})

names(common_sites) <- tx_common

# remove the transcripts for which there aren't analysed sites in common between the different conditions (NA)
common_sites <- common_sites[-unname(which(is.na(common_sites)))]
print(paste('There are', as.character(length(common_sites)), 'transcripts with at least one site with coverage higher than 30x shared by all the conditions'))

# function that returns the mean coverage along the 3' UTR of each transcript only considering the sites that can be analyzed in all the conditions 
Calculate_mean_cov_common_sites <- function(bam_file, regions_file, output_dir) {
  #import the regions_file as a GRanges object
  regions_df <- read.table(regions_file, quote = "")
  regions <- GRanges(seqnames = regions_df$V1,
                     ranges = IRanges(start = regions_df$V2,
                                      end = regions_df$V3),
                     strand = regions_df$V6)
  
  pu_par <- PileupParam(max_depth=1000000, min_base_quality=0, min_mapq=0,
                        min_nucleotide_depth=0, min_minor_allele_depth=0,
                        distinguish_strands=FALSE, distinguish_nucleotides=FALSE,
                        ignore_query_Ns=TRUE, include_deletions=FALSE, include_insertions=FALSE,
                        left_bins=NULL, query_bins=NULL, cycle_bins=NULL)
  
  pileup_tot <- pileup(file = bam_file, index = paste0(bam_file, ".bai"), pileupParam = pu_par)
    
  pileup_GR <- GRanges(seqnames = pileup_tot$seqnames,
                       ranges = IRanges(start = pileup_tot$pos, end = pileup_tot$pos))
  mcols(pileup_GR) <- pileup_tot$count
  colnames(mcols(pileup_GR)) <- "coverage"
  
  hits <- findOverlaps(query = pileup_GR, subject = regions, ignore.strand = TRUE)
  coverage_split <- split(pileup_GR[queryHits(hits)], subjectHits(hits))
  regions_names <- unname(unlist(lapply(split(as.character(seqnames(regions[subjectHits(hits)])), subjectHits(hits)), '[[', 1)))
  names(coverage_split) <- regions_names
  regions_names_all <- as.vector(seqnames(regions))
  
  # extract the transcripts for which there are analyzable sites in all the conditions
  coverage_split2 <- coverage_split[names(common_sites)]

  mean_cov_common_sites <- lapply(coverage_split2, function(x){
    tx <- unique(as.vector(seqnames(x)))
    mean(x[unlist(unname(common_sites[tx]))]$coverage)
  })
  return(mean_cov_common_sites)
}

# compute the mean coverage along the 3' UTR of each transcript only considering the sites that can be analyzed in all the conditions
mean_cov_common_sites_IVTprom <- Calculate_mean_cov_common_sites(path_to_IVTprom, regions_file, output_dir)

mean_cov_common_sites_WTmin <- Calculate_mean_cov_common_sites(path_to_WTmin, regions_file, output_dir)

mean_cov_common_sites_WTprom <- Calculate_mean_cov_common_sites(path_to_WTprom, regions_file, output_dir)

pdf(output_pdf_file)
par(mfrow=c(2,2))
hist(unlist(mean_cov_common_sites_IVTprom), breaks = 100, main = paste0('IVT PromethION\nmin_mean_coverage=',as.character(round(min(unlist(mean_cov_common_sites_IVTprom)),0))),xlab = 'Mean coverage',ylim = c(0,12),col='gray')
hist(unlist(mean_cov_common_sites_WTmin), breaks = 100, main = paste0('WT MinION\nmin_mean_coverage=',as.character(round(min(unlist(mean_cov_common_sites_WTmin)),0))),xlab = 'Mean coverage',ylim = c(0,12),col='gray')
hist(unlist(mean_cov_common_sites_WTprom), breaks = 100, main = paste0('WT PromethION\nmin_mean_coverage=',as.character(round(min(unlist(mean_cov_common_sites_WTprom)),0))),xlab = 'Mean coverage',ylim = c(0,12),col='gray')
dev.off()

# saturation
mean_cov_common_sites_IVTprom_sat <- unlist(mean_cov_common_sites_IVTprom)
mean_cov_common_sites_IVTprom_sat[which(mean_cov_common_sites_IVTprom_sat > as.numeric(saturation_level))] <- as.numeric(saturation_level)
mean_cov_common_sites_WTmin_sat <- unlist(mean_cov_common_sites_WTmin)
mean_cov_common_sites_WTmin_sat[which(mean_cov_common_sites_WTmin_sat > as.numeric(saturation_level))] <- as.numeric(saturation_level)
mean_cov_common_sites_WTprom_sat <- unlist(mean_cov_common_sites_WTprom)
mean_cov_common_sites_WTprom_sat[which(mean_cov_common_sites_WTprom_sat > as.numeric(saturation_level))] <- as.numeric(saturation_level)

pdf(output_pdf_file_saturation)
par(mfrow=c(2,2))
hist(unlist(mean_cov_common_sites_IVTprom_sat), breaks = 100, main = paste0('IVT PromethION\nmin_mean_coverage=', as.character(round(min(mean_cov_common_sites_IVTprom_sat),0))),xlab = 'Mean coverage',ylim = c(0,10),col='gray')
hist(unlist(mean_cov_common_sites_WTmin_sat), breaks = 100, main = paste0('WT MinION\nmin_mean_coverage=', as.character(round(min(mean_cov_common_sites_WTmin_sat),0))),xlab = 'Mean coverage',ylim = c(0,10),col='gray')
hist(unlist(mean_cov_common_sites_WTprom_sat), breaks = 100, main = paste0('WT PromethION\nmin_mean_coverage=',as.character(round(min(mean_cov_common_sites_WTprom_sat),0))),xlab = 'Mean coverage',ylim = c(0,10),col='gray')
dev.off()





