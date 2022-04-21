TGP Tour
================
Alexandra Tidd
October 28, 2021

## predict\_target\_genes()

TGP is a package to predict the target genes of fine-mapped variants of a trait.

`predict_target_genes()` is the master, user-facing function of this package. All other functions are helper functions called by `predict_target_genes()`. The functions are all documented, so once the package is loaded you can access the help pages for individual functions, which explain all the arguments.

``` r
?predict_target_genes()
```

### Package data

This package uses both reference genomic annotation datasets and user-provided trait-specific datasets. Some data that accompanies the package must be downloaded separately. Genomic coordinates use the hg19 build. The package's BED handling follows UCSC's `bedtools` conventions, so it expects 0-based start positions and 1-based end positions.

The reference panels and example trait data are saved locally...

``` bash
package_data_dir=/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/
reference_panels_dir=$package_data_dir/reference_panels/output/
traits_dir=$package_data_dir/traits/output/
```

#### Reference data

##### Internal reference data

Smaller generic reference datasets, including chromosome sizes, GENCODE annotations and REVEL annotations (`ChrSizes`, `TSSs`, `exons`, `introns`, `promoters`, `missense`, `nonsense`, `splicesite`) are stored internally as parsed objects in `R/sysdata.R`. They are accessible when the package is loaded, but not visible due to lazy loading. The raw data from which these objects are derived are in the reference panels data directory and the scripts to generate the output are in the reference panels code directory.

##### External reference data

Larger cell-type-specific reference panels are stored externally as local files in the reference panels output directory. These are too large to upload to GitHub and would make the package too bulky, so they will be published in a directory to be downloaded alongside the package. The reproducible scripts to generate these files are in the reference panels code directory.

``` r
list.files("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/output/")
```

    ## [1] "DHSs"         "expression"   "GENCODE"      "H3K27ac"      "HiChIP"      
    ## [6] "metadata.tsv" "REVEL"        "TADs"

#### User-provided data

There is one required user-provided file for the `predict_target_genes()` function: the trait variants (the `variants_file` argument). Known genes for the trait (the `known_genes_file` argument) are needed only if `do_performance = T`. Example inputs can be found in the traits directory for several tested traits.

``` bash
(cd /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/ ; ls -d */)
```

    ## BC_Michailidou2017_FM/
    ## BC_Michailidou2017_LD/
    ## BloodCancer_Law2017andWent2018_LD/
    ## CRC_Law2019_LD/
    ## EC_Wang2022_LD/
    ## EOC_Jones2020_FM/
    ## IBD_Huang2017_FM/
    ## IBD_Huang2017_LD/
    ## Lymphoma_Sud2018_LD/
    ## PrCa_Dadaev2018_FM/
    ## PrCa_Giambartolomei2021_FM/
    ## PrCa_Giambartolomei2021_LD/
    ## PrCa_GWASCatalog_LD/
    ## PrCa_Schumacher2018_LD/

##### Trait variants

The variants file should be a BED file with metadata columns for the variant name and the credible set to which it belongs.

``` bash
head /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/BC_Michailidou2017_FM/variants.bed
```

    ## chr1 10551762    10551763    rs657244:10551763:A:G   BCAC_FM_1.1
    ## chr1 10563363    10563364    rs202087283:10563364:G:A    BCAC_FM_1.1
    ## chr1 10564674    10564675    chr1_10564675_A_G   BCAC_FM_1.1
    ## chr1 10566521    10566522    rs617728:10566522:C:T   BCAC_FM_1.1
    ## chr1 10569000    10569000    rs60354536:10569000:C:CT    BCAC_FM_1.1
    ## chr1 10569257    10569258    rs2480785:10569258:G:A  BCAC_FM_1.1
    ## chr1 10579544    10579545    rs1411402:10579545:G:T  BCAC_FM_1.1
    ## chr1 10580890    10580891    chr1_10580891_C_T   BCAC_FM_1.1
    ## chr1 10581050    10581051    rs2506885:10581051:A:T  BCAC_FM_1.1
    ## chr1 10581657    10581658    rs2056417   BCAC_FM_1.1

##### Trait known genes

The known genes file should be a text file with a single column of known gene symbols. These symbols must be GENCODE-compatible.

``` bash
head /working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/BC_Michailidou2017_FM/known_genes.txt
```

    ## AKT1
    ## ARID1A
    ## ATM
    ## BRCA1
    ## BRCA2
    ## CBFB
    ## CDH1
    ## CDKN1B
    ## CHEK2
    ## CTCF

### Running predict\_target\_genes()

To run `predict_target_genes()`...

``` r
package_data_dir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/"
MA <- predict_target_genes(
  trait = "BC_Michailidou2017_FM",
  variants_file = paste0(package_data_dir, "traits/output/BC_Michailidou2017_FM/variants.bed"),
  known_genes_file = paste0(package_data_dir, "traits/output/BC_Michailidou2017_FM/known_genes.txt"),
  reference_panels_dir = paste0(package_data_dir, "reference_panels/output/")
  )
```

Unless an `out_dir` argument is passed, the results will be saved to "out/${trait}/${celltypes}/". If `do_timestamp = T`, then the run results will be saved to a timestamped subdirectory.

If you are calling `predict_target_genes()` repeatedly in the same session, you can load the large reference objects `H3K27ac` and `HiChIP` into the global environment once, and then pass them to the function pre-loaded. This prevents redundant re-loading with each call to `predict_target_genes()`.

``` r
# 1. load large reference panel objects
package_data_dir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/"
HiChIP <- readRDS(paste0(package_data_dir, "reference_panels/output/HiChIP/HiChIP.rds"))
H3K27ac <- readRDS(paste0(package_data_dir, "reference_panels/output/H3K27ac/H3K27ac.rds"))
# 2. pass objects to predict_target_genes for a quicker runtime
MA <- predict_target_genes(
  trait = "BC_Michailidou2017_FM",
  variants_file = paste0(package_data_dir, "traits/output/BC_Michailidou2017_FM/variants.bed"),
  known_genes_file = paste0(package_data_dir, "traits/output/BC_Michailidou2017_FM/known_genes.txt"),
  reference_panels_dir = paste0(package_data_dir, "reference_panels/output/"),
  HiChIP = HiChIP,
  H3K27ac = H3K27ac
  )
```

## Interactive development

To simulate installing and loading the package during interactive development...

``` r
setwd("/working/lab_jonathb/alexandT/tgp/")
devtools::load_all()
```

This package was written using package development conventions from <https://r-pkgs.org/>.

All of the arguments passed to `predict_target_genes()` in a given run are written to `arguments_for_predict_target_genes.R` in the output directory. If you wish to restore the internal environment of a run in your global environment to run the `R/predict_target_genes.R` script...

``` r
args <- dget("out/BC_Michailidou2017_FM/enriched_tissues/arguments_for_predict_target_genes.R") 
list2env(args, envir=.GlobalEnv)
```

All functions in the package, unless generic, use the same object names as `predict_target_genes()`, so you can also run the internal code of the helper functions directly as long as you have already run the internal code of `predict_target_genes()` up to the point at which that helper function is called. This allows you to run the complete pipeline line-by-line for package debugging and development.
