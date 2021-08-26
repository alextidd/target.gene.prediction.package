## code to prepare `BCVars` dataset goes here

# Import a BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants

BCVarsBedFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/Traits/BC/BC.VariantList.bed"
BCVars <- target.gene.prediction.package::import_BED(bedfile = BCVarsBedFile,
                                                     metadata_cols = c("variant", "cs", "trait", "status"))

usethis::use_data(BCVars, overwrite = TRUE)
