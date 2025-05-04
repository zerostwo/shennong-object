# data-raw/so.R
# Generate a compact Shennong demo object with realistic assays, layers, reductions, and metadata

library(ShennongObject)
set.seed(717)

# Define gene/peak/sample names
genes <- paste0("gene", 1:200)
peaks <- paste0("chr1:", seq(1, 2000, by = 10), "-", seq(100, 2090, by = 10)) # 200 peaks
samples <- paste0("sample", 1:20)

# Construct RNA assay with sparsity
rna_counts <- matrix(
  ifelse(runif(4000) < 0.3, 0, rpois(4000, lambda = 10)),
  nrow = 200,
  dimnames = list(genes, samples)
)
rna_data <- log1p(rna_counts)
rna_scale <- scale(rna_data)

# Construct ATAC assay with peak-style rownames
atac_counts <- matrix(
  ifelse(runif(4000) < 0.4, 0, rpois(4000, lambda = 5)),
  nrow = 200,
  dimnames = list(peaks, samples)
)

# Construct sample-level metadata
metadata <- data.frame(
  condition = rep(c("Normal", "Tumor"), each = 10),
  gender = rep(c("Male", "Female"), 10),
  age = round(runif(20, 40, 70)),
  row.names = samples
)

# Create Shennong object
so <- sn_initialize_shennong_object(
  counts = rna_counts,
  metadata = metadata,
  assay = "RNA", project = "Shennong")

# Add RNA layers
sn_layer_data(so, layer = "data", assay = "RNA") <- rna_data
sn_layer_data(so, layer = "scale.data", assay = "RNA") <- rna_scale
sn_active_layer(so) <- "data"

# Add ATAC assay
so <- sn_add_assay(so, data = atac_counts, assay = "ATAC")

# Set sample identities by condition
sn_active_ident(so) <- "condition"

# PCA reduction on scaled RNA data
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

# MDS reduction on scaled RNA data
dist_mat <- dist(t(rna_scale))
mds_embed <- cmdscale(dist_mat, k = 2)

so[["mds"]] <- sn_create_reduction(
  embedding = mds_embed,
  assay_used = "RNA",
  key = "MDS"
)

# Save to package
usethis::use_data(so, overwrite = TRUE)
