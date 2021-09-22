# Documentation for deleted internal data objects

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

#' Mark specificity-/signal-binned DHSs across cell types
#'
#' A list of valr-compatible tibble BED objects imported from BED files of the DHSs binned by H3K27ac/H3K4me1
#' mark signal/specificity for each cell type.
#' (Internal data example)
#'
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
#' @source {/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/}
"DHSs"

#' Metadata for `DHSs`
#'
#' Metadata matching each DHS celltype code to its name, tissue and developmental stage
#'
#' \describe{
#'   \item{code}{Shortened cell type code in `DHSs` CellType column}
#'   \item{name}{Full, descriptive name of the cell type}
#'   \item{tissue}{Tissue of origin of the sample}
#'   \item{stage}{Developmental stage of the sample (embryo / fetus / adult)}
#' }
"DHSs_metadata"

#' Metadata for the contact data samples that accompany the package
#'
#' Metadata matching each contact assay in the contact list to its study, celltype and tissue
#'
#' \describe{
#'   \item{code}{Shortened cell type code in `DHSs` CellType column}
#'   \item{name}{Full, descriptive name of the cell type}
#'   \item{tissue}{Tissue of origin of the sample}
#'   \item{stage}{Developmental stage of the sample (embryo / fetus / adult)}
#' }
"contact_metadata"

#' Data frame of list(s) of driver gene symbols, which are treated as the positive control causal gene sets
#'
#' @format Tibble
#' \describe{
#'   \item{DriversList}{Name of the drivers list}
#'   \item{symbol}{Gene symbol of causal genes}
#' }
#' @source {/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/}
"drivers"
