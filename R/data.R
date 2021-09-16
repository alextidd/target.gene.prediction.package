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
