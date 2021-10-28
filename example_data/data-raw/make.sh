WKDIR=/working/lab_georgiat/alexandT/target.gene.prediction.package/example_data/data/

# ======================================================================
# Example data:
# -> variants
# BC
BCVarsFile="/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/Traits/BC/BC.VariantList.bed"
cut -f1-5 $BCVarsFile > $WKDIR/BC.VariantList.bed

# IBD
IBDVarsFile="/working/lab_georgiat/alexandT/target_gene_prediction/data/IBD/CCVs/IBD.variant.list.txt"
cat $IBDVarsFile | sed '1d' | awk -F"\t" 'BEGIN{OFS="\t";} {print $1,$2,$2,$3,$4}' \
> $WKDIR/IBD.VariantList.bed

# CRC
CRCVarsFile="/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/data/reference_data/ccv_beds/CRC.bed"
cat $CRCVarsFile \
> $WKDIR/CRC.VariantList.bed

# -> drivers
driversFile="/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/breast_cancer_drivers_2021.txt"
cp $driversFile $WKDIR/breast_cancer_drivers_2021.txt
