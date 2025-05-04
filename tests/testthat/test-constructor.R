test_that("Shennong object constructor works with minimal input", {
  counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))

  so <- sn_initialize_shennong_object(counts = counts)
  expect_s4_class(so, "Shennong")
  expect_equal(sn_assays(so), "RNA")
  expect_equal(dim(sn_layer_data(so, layer = "counts")), dim(counts))
})

test_that("Shennong object constructor handles metadata correctly", {
  counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))

  metadata <- data.frame(
    condition = sample(c("A", "B"), 10, replace = TRUE),
    batch = sample(c("X", "Y"), 10, replace = TRUE),
    row.names = colnames(counts)
  )

  so <- sn_initialize_shennong_object(counts = counts, metadata = metadata)
  expect_s4_class(so, "Shennong")
  expect_equal(dim(sn_metadata(so)), c(10, 5)) # 2 metadata + 2 QC + 1 sample_id
})

test_that("Filtering by min_counts and min_features works", {
  counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  counts[, 1:2] <- 0  # make first two samples bad
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))

  so <- sn_initialize_shennong_object(counts = counts, min_counts = 100, min_features = 10)
  expect_s4_class(so, "Shennong")
  expect_lt(ncol(sn_layer_data(so, "counts")), 10)
})

test_that("Shennong constructor gives informative errors", {
  bad_counts <- matrix(rpois(1000, 10), 100, 10)
  # Missing colnames
  expect_error(sn_initialize_shennong_object(counts = bad_counts[, -1]),
               "Column names are required")

  bad_counts_dup <- bad_counts
  colnames(bad_counts_dup) <- rep("sample1", 10)
  expect_error(sn_initialize_shennong_object(counts = bad_counts_dup),
               "Duplicate sample names")
})

test_that("sn_create_reduction returns valid ShennongReduction object", {
  emb <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(emb) <- c("PC1", "PC2")
  red <- sn_create_reduction(embedding = emb, assay_used = "RNA", key = "PC")
  expect_s4_class(red, "ShennongReduction")
  expect_equal(ncol(sn_embeddings(red)), 2)
  expect_true(grepl("^PC_", colnames(sn_embeddings(red))[1]))
})

test_that("sn_create_reduction handles key inference and formatting", {
  emb <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(emb) <- c("UMAP1", "UMAP2")
  red <- sn_create_reduction(embedding = emb, assay_used = "RNA")
  expect_s4_class(red, "ShennongReduction")
  expect_match(colnames(red@embedding)[1], "^UMAP_")
})

test_that("sn_create_reduction errors when key can't be inferred", {
  emb <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(emb) <- NULL
  expect_error(sn_create_reduction(embedding = emb), "must be provided")
})
