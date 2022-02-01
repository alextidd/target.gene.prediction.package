TGP Tour
================
Alexandra Tidd
October 28, 2021

## Interactive development

To simulate installing and loading the package during interactive development...

``` r
setwd("/working/lab_jonathb/alexandT/tgp/")
devtools::load_all()
```

This package was written using package development conventions from <https://r-pkgs.org/>. The functions are all documented, so once the package is loaded you can access the help pages for individual functions, which explain all the arguments. `predict_target_genes()` is the user-facing master function of this package. All other functions are helper functions called by `predict_target_genes()`.

``` r
?predict_target_genes()
```

If you also want the external reference data within the package directory, you can copy it over to mirror my structure...

``` bash
MyPackageDir=/working/lab_jonathb/alexandT/tgp/
YourPackageDir=path/to/your/copy/of/tgp/
mkdir -p $YourPackageDir/{reference_data,example_data}/data
cp -r $MyPackageDir/reference_data/data/* $YourPackageDir/reference_data/data/
cp -r $MyPackageDir/example_data/data/* $YourPackageDir/example_data/data/
```

To run `predict_target_genes()`...

``` r
# These paths work if you copied across my structure, as above. The default paths are full paths to my files, so should work the same.
MAE <- predict_target_genes(trait = "BC",
                            outDir = "out/",
                            variantsFile = "example_data/data/BC/BC.VariantList.bed",
                            driversFile = "example_data/data/BC/BC.Drivers.txt",
                            do_scoring = T,
                            do_performance = T,
                            do_XGBoost = T)
```

If you wish to run the internal code of the `R/predict_target_genes.R` script step-by-step for development and debugging, you can add its default arguments to your global environment...

``` r
setwd("/working/lab_jonathb/alexandT/tgp/") ; library(devtools) ; load_all() 
tissue_of_interest = NULL 
trait="BC" 
outDir = "out/" 
variantsFile="/working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.VariantList.bed" 
driversFile = "/working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.Drivers.txt" 
variant_to_gene_max_distance = 2e6 
min_proportion_of_variants_in_top_DHSs = 0.05 
include_all_celltypes_in_the_enriched_tissue = T 
do_all_cells = F 
do_manual_weighting = F 
n_unique_manual_weights = NULL 
do_scoring = F 
do_performance = F 
do_XGBoost = F 
contact = NULL 
DHSs = NULL
```

All functions in the package, unless generic, use the same object names as `predict_target_genes()`, so you can also run the internal code of the helper functions directly as long as you have already run the internal code of `predict_target_genes()` up to the point at which that helper function is called. This allows you to run the complete code line-by-line for debugging package functions.

If you are calling `predict_target_genes()` repeatedly in the same session, you can load the large reference objects `DHSs` and `contact` into the global environment once, and then pass them to the function pre-loaded. This prevents redundant re-loading with each call to `predict_target_genes()`.

``` r
# 1. load referenceDir objects
referenceDir = "/working/lab_jonathb/alexandT/tgp/reference_data/data/"
contact <- readRDS(paste0(referenceDir, "contact.rda"))
DHSs <- readRDS(paste0(referenceDir, "DHSs.rda"))
# 2. pass objects to predict_target_genes for quicker runtime
MAE <- predict_target_genes(trait = "BC",
                            outDir = "out/",
                            variantsFile = "example_data/data/BC/BC.VariantList.bed",
                            driversFile = "example_data/data/BC/BC.Drivers.txt",
                            do_scoring = T,
                            do_performance = T,
                            do_XGBoost = T,
                            contact = contact,
                            DHSs = DHSs)
```

## Package data

This package will use both reference genomic annotation datasets and user-provided trait-specific datasets.

### Reference data

#### Internal reference data

Smaller generic reference datasets, including chromosome sizes, GENCODE annotations and REVEL annotations (`ChrSizes`, `TSSs`, `exons`, `introns`, `promoters`, `missense`, `nonsense`, `splicesite`) are stored internally as parsed objects in `R/sysdata.R`. They are accessible when the package is loaded, but not visible due to lazy loading. The reproducible code to generate these objects is in `data-raw/sysdata.R`.

#### External reference data

Larger cell-type-specific reference datasets are stored as local files in `reference_data/data/`. These are too large to upload to GitHub and would make the package too bulky, so they will be published in a directory to be downloaded alongside the package. The reproducible code to generate these files is in `reference_data/data-raw/`

``` r
list.files("/working/lab_jonathb/alexandT/tgp/reference_data/data/")
```

    ## [1] "all_metadata.tsv"                        
    ## [2] "contact.rda"                             
    ## [3] "DHSs.rda"                                
    ## [4] "expression.tsv"                          
    ## [5] "specific_DHSs_closest_specific_genes.rda"
    ## [6] "TADs.rda"

### User-provided data

There are two required user-provided files for the `predict_target_genes()` function: the list of trait drivers (the `driversFile` argument) and the list of trait variants (the `variantsFile` argument). Example inputs can be found in my local `example_data/data/` directory. Full paths to these are set as the default arguments of the function.

#### Trait variants

The variants file should be a BED file with metadata columns for the variant name and the credible set to which it belongs.

``` bash
head /working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.VariantList.bed
```

    ## chr1 10551762    10551763    rs657244:10551763:A:G   BCAC_FM_1ichav1
    ## chr1 10563363    10563364    rs202087283:10563364:G:A    BCAC_FM_1ichav1
    ## chr1 10564674    10564675    chr1_10564675_A_G   BCAC_FM_1ichav1
    ## chr1 10566521    10566522    rs617728:10566522:C:T   BCAC_FM_1ichav1
    ## chr1 10569000    10569000    rs60354536:10569000:C:CT    BCAC_FM_1ichav1
    ## chr1 10569257    10569258    rs2480785:10569258:G:A  BCAC_FM_1ichav1
    ## chr1 10579544    10579545    rs1411402:10579545:G:T  BCAC_FM_1ichav1
    ## chr1 10580890    10580891    chr1_10580891_C_T   BCAC_FM_1ichav1
    ## chr1 10581050    10581051    rs2506885:10581051:A:T  BCAC_FM_1ichav1
    ## chr1 10581657    10581658    rs2056417   BCAC_FM_1ichav1

#### Trait drivers

The drivers file should be a file with a single column of driver gene symbols. These symbols should be GENCODE-compatible.

``` bash
head /working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.Drivers.txt
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
