# SeuratConverter

## Description
SeuratConverter is an R package to convert between Seurat objects and AnnData objects.




## Installation and Usage :

Install from CRAN-Compatible Source:

```r

install.packages(c("Seurat", "Matrix", "hdf5r"))
devtools::install_github("https://github.com/MJBioInfo/SeuratConverter")

```

## Convert Seurat → H5AD:

```r

library(SeuratConverter)
seurat_to_anndata(
  seurat_obj = pbmc_small,
  output_file = "pbmc.h5ad",
  reductions = c("pca", "umap")
)

```
## Convert H5AD → Seurat:

```r
seurat_obj <- anndata_to_seurat("pbmc.h5ad", output_file = "pbmc.rds" , assay_name = "RNA", layers = c("counts", "data"), reductions = NULL)

```
