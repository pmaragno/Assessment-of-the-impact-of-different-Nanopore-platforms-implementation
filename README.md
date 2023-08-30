# Assessment of the impact of different Nanopore platforms implementation

This analysis is inserted in the context of the description of the In Vitro Transcribed (IVT) RNA sample and its use as reference sample in Nanopore comparative tools for the detection of RNA modifications on the base of Nanopore direct RNA sequencing data. 
The advantage of the use of the IVT sample as reference is that it is completely unmodified RNA that allows an effective lack of RNA modifications with 
respect to a sample subjected to the knockdown (KD) or knockout (KO) for a modification writer enzyme.

This analysis consists of the detection of RNA modifications on the 3' UTR of 50 transcripts in a wild type sample using as reference an IVT sample, all coming from K562 cells. To study if there is an effect due to the implementation of different Nanopore platforms for the sequencing of the test and the reference samples, different comparisons were analysed:
* WT sample sequenced on a MinION flow-cell of GridION platform vs IVT sample sequenced on PromethION platform;
* WT sample vs IVT sample, both sequenced on PromethION platform.

The detection of RNA modifications was performed in parallel with two different comparative tools, ELIGOS and Nanocompore, to address the reproducibility of the results using different machine learning approaches.

Furthermore, the analysis was also implemented capping the coverage on the 3’ UTR of the 50 selected transcripts to 100 to assess if there is an impact linked to the different level of transcript coverage.  

## Getting started
link to gtf annotation file: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
link to transcriptome fasta file:

## Workflow

Raw fast5 files sequenced on PromethION platform for K562 IVT and K562 WT samples were base-called with command “guppy_basecaller -i <fast5 directory> -r -x 'auto' -s <output directory> --fast5_out -c rna_r9.4.1_70bps_hac_prom.cfg” using Guppy v6.2.1 and Guppy v6.4.6, respectively. Similarly, raw fast5 files sequenced on a MinION flow-cell of GridION platform for K562 WT sample were base-called with Guppy v6.2.1 with command “guppy_basecaller -i <fast5 directory> -r -x 'auto' -s <output directory> --fast5_out -c rna_r9.4.1_70bps_hac.cfg”. 
Reads from each sample were then aligned to the transcriptome with minimap2 v2.17 and only the reads aligned to the transcript strand were retained using samtools v1.6 with commands “minimap2 -ax map-ont -k 14 <reference_transcriptome.fasta> <reads.fastq> | samtools view -h -F2324 | samtools sort -o <filtered_reads_mapping_on_transcriptome_F2324.bam>”.

In selection_50_transcript.R script, filtered bam files for the three conditions were imported in R with pileup function of Rsamtools package v2.14.0 to compute the number of nucleotides on the 3’ UTR of each transcript with a coverage at least equal to 30x. For each transcript, the minimum number of nucleotides with coverage at least equal to 30x across the three conditions was evaluated, transcripts were ordered in decreasing order based on the minimum value and the top 50 transcripts were selected. 

### ELIGOS
At this point, we set out to call RNA modifications on K562 WT sample sequenced either on MinION or PromethION flow-cells, using IVT data sequenced on PromethION flow-cells as a baseline and using both ELIGOS v2.1.0 and Nanocompore v1.0.4 comparative tools.
In particular, for running ELIGOS, all the reads from each sample were aligned to the trascriptome (link to the trascriptome) with minimap2 v2.17 and alignments were filtered with samtools v1.6 with command “minimap2 -ax map-ont -k 14 <reference_transcriptome.fasta> <reads.fastq> | samtools view -hSb -F 2324| samtools sort -o <output.bam>”. 
ELIGOS was then run with command “eligos2 pair_diff_mod -tbam <output.bam> -cbam <IVT_filtered.bam> -reg <3UTR_50genes.bed> -ref <reference_ transcriptome.fasta> --pval 1 --oddR 0 --esb 0.2 --adjPval 1 -o <output directory>”. 

In postprocessing_analysis.R script, the output file <eligos2.combine.txt> obtained comparing the WT sample - sequenced either on MinION or on PromethION flow-cells - to the IVT baseline were imported in R using GRanges function of GenomicRanges package v1.50.2. The analysed sites were evaluated for both comparative analyses and the sites analysed in both comparisons were identified. Then, the hits were identified using as parameters cut-offs pval <= 0.05 and pvalAdj <= 1e-04 and oddR >= 2.5. This set of hits was then intersected with the nucleotides analysed in both comparisons, to obtain a set of filtered hits, which were eventually expanded to ranges of 10 nucleotides centered around each hit. The overlap between these ranges in the two comparisons was computed using overlapOfGRanges function of CompEpiTools package v1.32.0. Scripts for running ELIGOS pipeline in Nextflow8 framework are available at: https://github.com/MaestSi/nf-eligos.

### Nanocompore
While, for running Nanocompore, reads mapping to the 3’ UTR of the 50 selected transcripts were extracted using a combination of samtools v1.62 and seqtk v1.39  with commands “samtools view <filtered_reads_mapping_on_transcriptome_F2324.bam> -L <3UTR_50genes.bed> | cut -f1 | sort| uniq > <reads_3UTR_50tx.txt>” followed by “seqtk subseq <reads.fastq> <reads_3UTR_50tx.txt> > <reads_3UTR_50tx.fastq>”.
Nanocompore pipeline was then run using minimap2 v2.17, samtools v1.6, f5c v0.7 and Nanocompore v1.0.4 using the following scripts (link to the two Logan’s Nanocompore scripts).

In postprocessing_analysis.R script, the output files <outnanocompore_results.tsv> obtained comparing the WT sample - sequenced either on MinION or on PromethION flow-cells - to the IVT baseline were imported in R and the corresponding GRanges objects containing the genomic coordinates of each hit called by Nanocompore were created using GRanges function of GenomicRanges package v1.50.2. The same steps done for ELIGOS to identify analysed sites in both comparisons, all the hits and those in sites analysed in both comparisons were done also for Nanocompore - for which the analysable sites are those with at least 30x coverage - paying attention to limit the analysis only at the 3’ UTR of the 50 selected transcripts and using as parameters to identify the hits a GMM_logit_pvalue < 0.05 and |logit_LOR| > 0.5. 
Eventually, the overlap between hits called by ELIGOS and Nanocompore using either of the WT test conditions was computed as previously described.

The analysis was then repeated on a random subset of reads, obtained capping the coverage on the 3’ UTR of the 50 selected transcripts to 100x, using samNormalise.pl script with the commands: “minimap2 -ax map-ont -k 14 <reference_transcriptome.fasta> <reads_3UTR_50tx.fastq> | samtools view -h -F2324 | samtools sort -o <filtered_reads_3UTR_50tx_mapping_on_transcriptome_F2324.bam>”; “samtools view <filtered_reads_3UTR_50tx_mapping_on_transcriptome_F2324.bam> | samNormalise.pl -coverage 100 -format fastq > <reads_3UTR_50tx_cov100.fastq>”. 

Coverage was fixed to 100X on the base of the analysis done in setting_max_coverage.R script.

## Results
### ELIGOS results using all the reads on the 3' UTR of the 50 selected transcripts

| ELIGOS | WT MinION vs IVT PromethION | WT PromethION vs IVT PromethION | Common |
| ------------- | ------------- |------------- | ------------- |
| # analysed nucleotides in 3UTR | 3,952 (25 tx)  | 3,121 (25 xt) | 2,722 (25 tx) |
| # filtered hits in 3UTR | 103 (18 tx) | 78 (19 tx) | - | 
| # filtered hits in 3UTR analysed in both comparisons | 88 (18 tx) | 77 (18 tx) | 70 (17 tx) |  

![WTmin_vs_IVTprom_WTprom_vs_IVTprom_overlap_hits_in_common_analysed_sites_eligos copy-1](https://github.com/pmaragno/Assessment-of-the-impact-of-different-Nanopore-platforms-implementation/assets/103447655/bd57f212-80da-4785-88ac-99ece2306755)

### Nanocompore results using all the reads on the 3' UTR of the 50 selected transcripts

| Nanocompore | WT MinION vs IVT PromethION | WT PromethION vs IVT PromethION | Common |
| ------------- | ------------- |------------- | ------------- |
| # analysed nucleotides in 3UTR | 13,776 (28 tx) | 19,701 (36 tx) | 13,776 (28 tx) |
| # filtered hits in 3UTR | 81 (21 tx) | 169 (30 tx) | - | 
| # filtered hits in 3UTR analysed in both comparisons | 81 (21 tx) | 137 (24 tx) | 64 (18 tx) |  

<img src="https://github.com/pmaragno/Assessment-of-the-impact-of-different-Nanopore-platforms-implementation/assets/103447655/706bd1cd-a4a5-49da-8b31-e7c6258cb800" width="100" height="100" />

![WTmin_vs_IVTprom_WTprom_vs_IVTprom_overlap_hits_in_common_analysed_sites_nano copy-1](https://github.com/pmaragno/Assessment-of-the-impact-of-different-Nanopore-platforms-implementation/assets/103447655/706bd1cd-a4a5-49da-8b31-e7c6258cb800)

### Overlap between ELIGOS and Nanocompore using all the reads on the 3' UTR of the 50 selected transcripts 
WT MinION vs IVT PromethION

WT PromethION vs IVT PromethION

