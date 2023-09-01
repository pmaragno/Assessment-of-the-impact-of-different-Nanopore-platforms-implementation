# Assessment of the impact of different Nanopore platforms implementation
This analysis is inserted in the context of the description of the In Vitro Transcribed (IVT) RNA sample and its use as reference sample in Nanopore comparative tools for the detection of RNA modifications on the base of Nanopore direct RNA sequencing data. 
The advantage of the use of the IVT sample as reference is that it is completely unmodified RNA that allows an effective lack of RNA modifications with 
respect to a sample subjected to the knockdown (KD) or knockout (KO) for a modification writer enzyme.

This analysis consists of the detection of RNA modifications on the 3' UTR of 50 transcripts in a wild type sample using as reference an IVT sample, both coming from K562 cells. To study if there is an effect due to the implementation of different Nanopore platforms for the sequencing of the test and the reference samples, different comparisons were analysed:
* WT sample sequenced on a MinION flow-cell of GridION platform vs IVT sample sequenced on PromethION platform;
* WT sample vs IVT sample, both sequenced on PromethION platform.

The detection of RNA modifications was performed in parallel with two different comparative tools, ELIGOS and Nanocompore, to address the reproducibility of the results using different machine learning approaches.

Furthermore, the analysis was implemented both using all the reads mapping on the 3' UTR of the 50 selected transcripts and capping the coverage on the 3’ UTR of the these transcripts to 100x to assess if there is an impact linked to the different transcript coverage.  

## Usage
link to gtf annotation file: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz

link to transcriptome fasta file:

```
Usage of selection_transcripts.R:
Rscript selection_transcripts.R gtf_file path_bed_3utr_tx num_transcripts path_bed_3utr_top_tx output_dir path_to_IVTprom path_to_WTmin path_to_WTprom

gtf_file                                  path to the gtf annotation file
path_bed_3utr_tx                          path to the directory where saving the bed file with the 3' UTR coordinates of all the transcripts
num_transcripts                           number of top trascripts to select
path_bed_3utr_top_tx                      path to the directory where saving the bed file with the 3' UTR coordinates of the selected transcripts
output_dir                                path to the output directory where saving the Rda intermediates files
path_to_IVTprom                           path to IVT PromethION bam file
path_to_WTmin                             path to WT MinION bam file
path_to_WTprom                            path to WT PromethION bam file
```

```
Usage of post_processing_analysis_IVT.R:
hits_ELIGOS <- ELIGOS_results(path_input, path_output)

path_input                                          path to the directory containing the <eligos2.combine.txt> output files of ELIGOS analysis
path_output                                         path to the directory where saving the plots with the overlapping between different comparisons

hits_ELIGOS_granges <- hits_ELIGOS[[1]]             For each comparison, a GRange with the hits in sites analysable in both comparisons, each range is expanded of 10 nucleotides centered around the corresponding hit
hits_ELIGOS_values <- hits_ELIGOS[[2]]              For each comparison, the number of hits in sites analysable in both comparisons

hits_Nanocompore <- Nanocompore_results(path_input, path_output, path_bed_3utr_top_tx)

path_input                                          path to the directory containing the <outnanocompore_results.tsv> output files of Nanocompore analysis
path_output                                         path to the directory where saving the plots with the overlapping between different comparisons
path_bed_3utr_top_tx                                path to the directory containing the bed file with the 3' UTR coordinates of the selected transcripts

hits_Nanocompore_granges <- hits_Nanocompore[[1]]   For each comparison, a GRange with the hits in sites analysable in both comparisons, each range is expanded of 10 nucleotides centered around the corresponding hit
hits_Nanocompore_values <- hits_Nanocompore[[2]]    For each comparison, the number of hits in sites analysable in both comparisons
```

```
Usage of mean_coverage_distribution.R:
Rscript mean_coverage_distribution.R path_to_IVTprom path_to_WTmin path_to_WTprom path_bed_3utr_top_tx output_dir output_pdf_file saturation_level output_pdf_file_saturation

path_to_IVTprom                           path to IVT PromethION bam file
path_to_WTmin                             path to WT MinION bam file
path_to_WTprom                            path to WT PromethION bam file
path_bed_3utr_top_tx                      path to the directory containing the bed file with the 3' UTR coordinates of the selected transcripts
output_dir                                path to the output directory where saving the Rda intermediates files
output_pdf_file                           path to the pdf file where saving the histogram with the mean coverage distribution for each condition
saturation_level                          value of maximum mean coverage for the saturation
output_pdf_file_saturation                path to the pdf file where saving the histogram with the mean coverage distribution for each condition after saturation
```
