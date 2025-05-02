# data-raw/so.R
# Generate a compact Shennong demo object with multiple assays, layers, reductions, and metadata

library(ShennongObject)
set.seed(717)

# Define gene and sample names
genes <- paste0("gene", 1:200)
samples <- paste0("sample", 1:20)

# --- Construct RNA assay ---
rna_counts <- matrix(rpois(4000, lambda = 10), nrow = 200, dimnames = list(genes, samples))
rna_data <- log1p(rna_counts)
rna_scale <- scale(rna_data)

# --- Construct ATAC assay ---
atac_counts <- matrix(rpois(4000, lambda = 5), nrow = 200, dimnames = list(genes, samples))

# --- Create Shennong object with RNA assay and layers ---
so <- sn_create_shennong_object(counts = rna_counts, assay = "RNA", project = "ExampleProject")
sn_layer_data(so, layer = "data", assay = "RNA") <- rna_data
sn_layer_data(so, layer = "scale.data", assay = "RNA") <- rna_scale

sn_active_layer(so) <- "data"

# --- Add ATAC assay ---
so <- sn_add_assay(so, data = atac_counts, assay = "ATAC")

# --- Add sample-level metadata ---
meta <- data.frame(
  condition = rep(c("Normal", "Tumor"), each = 10),
  gender = rep(c("Male", "Female"), 10),
  age = round(runif(20, 40, 70)),
  row.names = samples
)
so <- sn_add_metadata(so, meta)

# --- Run PCA on scaled RNA data ---
pca_res <- prcomp(t(rna_scale), center = TRUE, scale. = FALSE)
pca_embed <- pca_res$x
pca_loadings <- pca_res$rotation
pca_stdev <- pca_res$sdev

so[["pca"]] <- sn_create_reduction(
  embedding = pca_embed,
  loadings = pca_loadings,
  stdev = pca_stdev,
  assay_used = "RNA",
  key = "PC"
)

# --- Run MDS on scaled RNA data ---
dist_mat <- dist(t(rna_scale))
mds_embed <- cmdscale(dist_mat, k = 2)

so[["mds"]] <- sn_create_reduction(
  embedding = mds_embed,
  assay_used = "RNA",
  key = "MDS"
)

# --- Save to package ---
usethis::use_data(so, overwrite = TRUE)
