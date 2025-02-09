#' Convert Seurat Object to AnnData H5AD File
#'
#' Converts a Seurat object to an H5AD file format compatible with Scanpy/AnnData,
#' preserving counts, metadata, and dimensionality reductions.
#'
#' @param seurat_obj A Seurat object to convert
#' @param output_file Path to output H5AD file (e.g., "data.h5ad")
#' @param assay Name of assay to convert (default: "RNA")
#' @param slots Seurat assay slots to include as H5AD layers (default: c("counts", "data", "scale.data"))
#' @param reductions Dimensionality reductions to include (default: c("pca", "umap"))
#' @param metadata_fields Metadata columns to include (NULL = all columns)
#'
#' @return Invisibly returns NULL, writes H5AD file to disk
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' seurat_to_anndata(pbmc_small, "pbmc.h5ad")
#' }
#'
#' @seealso \code{\link{anndata_to_seurat}}
seurat_to_anndata <- function(
    seurat_obj,
    output_file,
    assay = "RNA",
    slots = c("counts", "data", "scale.data"),
    reductions = c("pca", "umap"),
    metadata_fields = NULL
) {
  # Validate dependencies
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Install hdf5r: install.packages('hdf5r')")
  }
  
  # Extract data from Seurat
  layers <- list()
  for (slot in slots) {
    slot_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
    if (inherits(slot_data, "dgCMatrix")) {
      layers[[slot]] <- slot_data
    } else {
      layers[[slot]] <- as(as.matrix(slot_data), "CsparseMatrix")
    }
  }
  
  # Metadata
  obs <- seurat_obj@meta.data[, metadata_fields, drop = FALSE]
  var <- Seurat::GetAssay(seurat_obj, assay)@meta.data
  
  # Reductions
  obsm <- list()
  for (reduction in reductions) {
    if (reduction %in% names(seurat_obj@reductions)) {
      obsm[[paste0("X_", reduction)]] <- seurat_obj@reductions[[reduction]]@cell.embeddings
    }
  }
  
  # Write H5AD using hdf5r
  h5 <- hdf5r::H5File$new(output_file, mode = "w")
  tryCatch({
    # Write X (use 'data' slot as default)
    x <- layers[["data"]]
    h5$create_dataset("X", x@x, dims = dim(x), chunk_dims = c(1000, 1000))
    
    # Write layers
    if (length(layers) > 0) {
      h5_layers <- h5$create_group("layers")
      for (layer_name in names(layers)) {
        layer_data <- layers[[layer_name]]
        h5_layers$create_dataset(layer_name, layer_data@x, dims = dim(layer_data))
      }
    }
    
    # Write obs/var
    h5$create_group("obs")$write(obs)
    h5$create_group("var")$write(var)
    
    # Write reductions (obsm)
    if (length(obsm) > 0) {
      h5_obsm <- h5$create_group("obsm")
      for (obsm_name in names(obsm)) {
        h5_obsm$create_dataset(obsm_name, obsm[[obsm_name]])
      }
    }
  }, finally = h5$close_all())
  
  message("Saved to H5AD: ", output_file)
}
