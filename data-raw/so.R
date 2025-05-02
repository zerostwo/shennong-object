# data-raw/so.R

library(ShennongObject)

counts <- matrix(
  rpois(200, lambda = 10),
  nrow = 20,
  dimnames = list(paste0("gene", 1:20), paste0("sample", 1:10))
)

so <- sn_create_shennong_object(counts = counts, assay = "RNA", project = "demo")

sn_layer_data(so, layer = "data") <- log1p(counts)[rev(rownames(counts)), ]
sn_layer_data(so, layer = "scale.data") <- scale(log1p(counts))

atac_counts <- matrix(
  rbinom(200, size = 1, prob = 0.2),
  nrow = 20,
  dimnames = list(paste0("peak", 1:20), paste0("sample", 1:10))
)
so <- sn_add_assay(so, data = atac_counts, assay = "ATAC")
sn_active_assay(so) <- "RNA"

sn_active_layer(so) <- "data"

sn_active_ident(so) <- factor(rep(c("Tumor", "Normal"), each = 5))

so <- sn_add_metadata(so, metadata = rep(c("M", "F"), 5), col_name = "sex")

usethis::use_data(so, overwrite = TRUE)
