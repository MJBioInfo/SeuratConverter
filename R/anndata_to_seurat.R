#' Convert AnnData H5AD File to Seurat Object
#'
#' Converts an H5AD file to a Seurat object, preserving counts, metadata,
#' and dimensionality reductions.
#'
#' @param h5ad_file Path to input H5AD file
#' @param output_file Optional path to save Seurat object as RDS file
#' @param assay_name Name for the created assay (default: "RNA")
#' @param layers H5AD layers to include as Seurat slots (default: c("counts", "data"))
#' @param reductions Dimensionality reductions to import (NULL = all available)
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- anndata_to_seurat("pbmc.h5ad")
#' }
#'
#' @seealso \code{\link{seurat_to_anndata}}
anndata_to_seurat <- function(
    h5ad_file,
    output_file = NULL,
    assay_name = "RNA",
    layers = c("counts", "data"),
    reductions = NULL
) {
  # Read H5AD
  h5 <- hdf5r::H5File$new(h5ad_file, mode = "r")
  
  # Read X, obs, var, obsm
  x <- h5[["X"]][, ]
  obs <- h5[["obs"]]$read()
  var <- h5[["var"]]$read()
  obsm <- list()
  if ("obsm" %in% names(h5)) {
    for (obsm_name in names(h5[["obsm"]])) {
      obsm[[obsm_name]] <- h5[["obsm"]][[obsm_name]][, ]
    }
  }
  
  # Read layers
  layer_data <- list()
  if ("layers" %in% names(h5)) {
    for (layer_name in layers) {
      layer_data[[layer_name]] <- h5[["layers"]][[layer_name]][, ]
    }
  }
  h5$close_all()
  
  # Build Seurat object
  counts <- layer_data$counts
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = t(counts),
    assay = assay_name,
    meta.data = obs
  )
  
  # Add layers
  for (slot in setdiff(names(layer_data), "counts")) {
    seurat_obj <- Seurat::SetAssayData(
      seurat_obj,
      slot = slot,
      new.data = t(layer_data[[slot]]),
      assay = assay_name
    )
  }
  
  # Add reductions
  for (reduction in names(obsm)) {
    key <- gsub("X_", "", reduction)
    seurat_obj@reductions[[key]] <- Seurat::CreateDimReducObject(
      embeddings = obsm[[reduction]],
      key = paste0(key, "_")
    )
  }
  
  # Save to RDS if requested
  if (!is.null(output_file)) {
    saveRDS(seurat_obj, output_file)
  }
  
  return(seurat_obj)
}