#' Import variants BED
#'
#' Import a BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#'
#' @param x A tab-delimited file in bed format with columns: chr, start, stop, snpid, locus
#'
#' @return A GRanges object representation of the variant BED file
#' @export
import_variants_GRanges <- function(x) {
  message(paste("Importing variants from", x))
  x <- target.gene.prediction.package::read_tibble(x, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
  v <- GenomicRanges::GRanges(
    # Intervals
    seqnames = x$V1,
    ranges = IRanges::IRanges(start = x$V2,
                              end = x$V3),
    # Metadata:
    snpid = x$V4,
    locus = x$V5
    )
  return(v)
}
