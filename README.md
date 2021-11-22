TGP Tour
================
Alexandra Tidd
October 28, 2021

## Interactive development

This package was written using package development conventions from <https://r-pkgs.org/>. The functions are all documented, so once the package is loaded you can access the help pages for individual functions, which explain all the arguments. `predict_target_genes()` is the user-facing master function of this package. All other functions are helper functions called by `predict_target_genes()`.

``` r
?predict_target_genes()
```

To simulate installing and loading the package during interactive development...

``` r
setwd("/working/lab_georgiat/alexandT/tgp/")
library(devtools)
load_all()
```

If you also want the external reference data within the package directory, you can copy it over to mirror my structure...

``` bash
MyPackageDir=/working/lab_georgiat/alexandT/tgp/
YourPackageDir=path/to/your/copy/of/tgp/
mkdir -p $YourPackageDir/{reference_data,example_data}/data
cp -r $MyPackageDir/reference_data/data/* $YourPackageDir/reference_data/data/
cp -r $MyPackageDir/example_data/data/* $YourPackageDir/example_data/data/
```

To run `predict_target_genes()`...

``` r
# These paths work if you copied across my structure, as above. The default paths are full paths to my files, so should work the same.
MAE <- predict_target_genes(trait = "BC",
                            outDir = "out/BC_enriched_cells/",
                            variantsFile = "example_data/data/BC.VariantList.bed",
                            driversFile = "example_data/data/BC.VariantList.bed",
                            do_scoring = T,
                            do_performance = T,
                            do_XGBoost = T)
```

If you wish to run the internal code of the `R/predict_target_genes.R` script step-by-step for development and debugging, you can add its default arguments to your global environment...

``` r
setwd("/working/lab_georgiat/alexandT/tgp/") ; library(devtools) ; load_all() 
tissue_of_interest = NULL 
trait="BC" 
outDir = "out/BC_enriched_cells/" ; variantsFile="/working/lab_georgiat/alexandT/tgp/example_data/data/BC.VariantList.bed" 
driversFile = "/working/lab_georgiat/alexandT/tgp/example_data/data/breast_cancer_drivers_2021.txt" 
referenceDir = "/working/lab_georgiat/alexandT/tgp/reference_data/data/" 
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
referenceDir = "/working/lab_georgiat/alexandT/tgp/external_data/reference/"
contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))
DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda"))
# 2. pass objects to predict_target_genes for quicker runtime
MAE <- predict_target_genes(trait = "BC",
                            outDir = "out/BC_enriched_cells/",
                            variantsFile = "example_data/data/BC.VariantList.bed",
                            driversFile = "example_data/data/BC.VariantList.bed",
                            do_scoring = T,
                            do_performance = T,
                            do_XGBoost = T,
                            contact = contact,
                            DHSs = DHSs)
```

## Package data

This package will both use reference genomic annotation datasets and user-provided trait-specific datasets.

### Internal reference data

Smaller reference datasets (`ChrSizes`, `TSSs`, `exons`, `introns`, `promoters`) are stored internally as parsed objects in `R/sysdata.R`. They are accessible when the package is loaded, but not visible due to lazy loading. The reproducible code to generate these objects is in `data-raw/sysdata.R`.

### External reference data

Larger reference datasets are stored as local files in `/working/lab_georgiat/alexandT/tgp/reference_data/data/`. These are too large to upload to GitHub and would make the package too bulky, so they will be published in a directory to be downloaded alongside the package. The reproducible code to generate these files is in `/working/lab_georgiat/alexandT/tgp/reference_data/data-raw/`

``` r
list.files("/working/lab_georgiat/alexandT/tgp/reference_data/data/")
```

## User-provided data

`/working/lab_georgiat/alexandT/tgp/example_data/data/`
