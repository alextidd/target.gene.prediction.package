#' Fine-mapped breast cancer SNPs from BCAC (user-input data example)
#'
#' A valr-compatible tibble BED object from a BED file of the strong fine-mapped breast cancer SNPs within credible set at breast cancer signals.
#' This is a template for the user-inputted trait variants file.
#'
#' @format A tibble BED object with 5369 intervals and 4 metadata columns:
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{variant}{Unique SNP ID (rsID / position and allele information)}
#'   \item{cs}{Credible set that the SNP belongs to}
#'   \item{trait}{Trait associated with the SNP from GWAS (breast cancer)}
#'   \item{status}{Specific breast cancer (estrogen receptor) phenotype associated with the SNP}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/Traits/BC/BC.VariantList.bed}
"BCVars"

#' ReMap TFBSs (internal data example)
#'
#' A valr-compatible tibble BED object imported from BED files of the transcription factor binding sites from different ReMap 2020 experiments
#'
#' @format A tibble BED object with 118576 ranges and 3 metadata columns:
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{Experiment}{Accession number of experiment from which the data was generated (e.g. GENCODE)}
#'   \item{TranscriptionFactor}{Transcription factor assayed}
#'   \item{CellType}{Cell type assayed}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/}
"remap"

#' DHSs binned by H3K27ac/H3K4me1 mark signal/specificity (internal data example)
#'
#' A valr-compatible tibble BED object imported from BED files of the DHSs binned by H3K27ac/H3K4me1 mark signal/specificity for each cell type
#' Data only given for MCF-7 currently
#'
#' @format A tibble BED object with 5369 ranges and 4 metadata columns:
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
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
#'
#' @format A list of paired BED-format tibbles with 359,970 interaction pairs (719,940 intervals) and 2 metadata columns:
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{InteractionID}{ID that links the paired ends between elements (first-last)}
#'   \item{OldPairID}{Old ID}
#'   \item{value}{Reads}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe}
"hichip"

#' Data frame of list(s) of driver gene symbols, which are treated as the positive control causal gene sets
#'
#' @format Tibble
#' \describe{
#'   \item{DriversList}{Name of the drivers list}
#'   \item{symbol}{Gene symbol of causal genes}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/}
"drivers"

#' A tibble BED object of transcription start sites of all transcripts (from GENCODE)
#'
#' @format Tibble
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{ensg}{ENSEMBL gene ID}
#'   \item{symbol}{Gene symbol of causal genes}
#'   \item{enst}{ENSEMBL transcript ID}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/GENCODE/gencode.v34lift37.basic.tss.bed.gz}
"TSSs"

#' A valr-compatible genome dataframe of chromosome sizes from hg19
#'
#' @format Tibble
#' \describe{
#'   \item{chrom}{chromosome}
#'   \item{size}{chromosome length}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/hg19.genome}
"ChrSizes"
