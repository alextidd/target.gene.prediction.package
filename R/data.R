#' Fine-mapped breast cancer SNPs from BCAC (user-input data example)
#'
#' A GRanges object from a BED file of the strong fine-mapped breast cancer SNPs within credible set at breast cancer signals.
#' This is a template for the user-inputted trait variants file.
#'
#' @format A GRanges object with 5369 ranges and 4 metadata columns:
#' \describe{
#'   \item{variant}{Unique SNP ID (rsID / position and allele information)}
#'   \item{cs}{Credible set that the SNP belongs to}
#'   \item{trait}{Trait associated with the SNP from GWAS (breast cancer)}
#'   \item{status}{Specific breast cancer (estrogen receptor) phenotype associated with the SNP}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/Traits/BC/BC.VariantList.bed}
"BCVars"

#' ReMap TFBSs (internal data example)
#'
#' A GRanges object imported from BED files of the transcription factor binding sites from different ReMap 2020 experiments
#'
#' @format A GRanges object with 118576 ranges and 3 metadata columns:
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start-end}
#'   \item{strand}{*/+/-}
#'   \item{Experiment}{Accession number of experiment from which the data was generated (e.g. GENCODE)}
#'   \item{TranscriptionFactor}{Transcription factor assayed}
#'   \item{CellType}{Cell type assayed}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/}
"remap"

#' DHSs binned by H3K27ac/H3K4me1 mark signal/specificity (internal data example)
#'
#' A GRanges object imported from BED files of the DHSs binned by H3K27ac/H3K4me1 mark signal/specificity for each cell type
#' Data only given for MCF-7 currently
#'
#' @format A GRanges object with 5369 ranges and 4 metadata columns:
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{start-end}
#'   \item{strand}{*/+/-}
#'   \item{CellType}{Cell type assayed}
#'   \item{Mark}{Histone modification mark assayed}
#'   \item{Bin}{What bin does this DHS rank in, according to the signal/specificity level of its marks?}
#'   \item{Binning}{Number of bins into which data is split (e.g. quartiles, deciles)}
#'   \item{Method}{Signal or specificity ranking of DHS sites}
#' }
#' @source {/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/bin_regions_by_specificity/sorted/H3K27ac/deciles/}
"bins"

#' Paired-end BEDPE list object c(dataframe of first end, dataframe of last end)
#'
#' A paired (first + last) list object representation of the MCF-7 HiChIP BEDPE file with an ID column matching pairs and metadata columns
#' Data only given for MCF-7 currently
#' This is not a GRanges object, as the GenomicRanges package does not appear to support paired-ends BED representation
#'
#' @format A list of paired BED-format dataframes with 359,970 pairs (719940, ends) and 2 metadata columns:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{group}{Group column}
#'   \item{value}{Reads}
#'   \item{PairID}{Numeric ID that links the paired ends in each element (first-last)}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe}
"hichip"
