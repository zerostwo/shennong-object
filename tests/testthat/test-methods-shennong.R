# Test for methods in methods-shennong.R

test_that("sn_assays returns correct assay names", {
  # Create a mock Shennong object using sn_create_shennong_object
  data("so")
  # Test
  expect_equal(sn_assays(so), c("RNA","ATAC"))
})

test_that("sn_active_assay gets and sets active assay", {
  data("so")
  # Test getter
  expect_equal(sn_active_assay(so), "RNA")

  # Test setter
  sn_active_assay(so) <- "ATAC"
  expect_equal(sn_active_assay(so), "ATAC")

  # Test invalid assay
  expect_error(sn_active_assay(so) <- "Invalid", "Assay 'Invalid' not found")
})

test_that("sn_layers returns correct layers", {
  data("so")
  # Test
  expect_equal(sn_layers(so), c("counts","data","scale.data"))
})

test_that("sn_active_layer gets and sets active layer", {
  data("so")

  # Test getter
  expect_equal(sn_active_layer(so), "data")

  # Test setter with valid layer
  sn_active_layer(so) <- "counts"
  expect_equal(sn_active_layer(so), "counts")

  # Test invalid layer
  expect_error(sn_active_layer(so) <- "invalid_layer", "Layer 'invalid_layer' not found")
})

test_that("sn_samples and sn_features return correct names", {
  data("so")

  # Test
  expect_equal(sn_samples(so), colnames(so))
  expect_equal(sn_features(so), rownames(so))
})

test_that("sn_active_ident gets and sets active ident", {
  data("so")

  # Test getter
  expect_equal(unname(sn_active_ident(so)), factor(rep("ExampleProject", 20)))

  # Test setter
  sn_active_ident(so) <- factor(rep("NewIdent", 20))
  expect_equal(unname(sn_active_ident(so)), factor(rep("NewIdent", 20)))
})

test_that("sn_project_name gets and sets project name", {
  data("so")

  # Test getter
  expect_equal(sn_project_name(so), "ExampleProject")

  # Test setter
  sn_project_name(so) <- "NewProject"
  expect_equal(sn_project_name(so), "NewProject")
})

test_that("sn_layer_data gets and sets layer data", {
  data("so")

  new_data <- log1p(sn_layer_data(so, layer = "counts"))
  sn_layer_data(so, layer = "data") <- new_data
  expect_equal(sn_layer_data(so, layer = "data"), new_data)

  rownames(new_data) <- paste0("InvalidGene", 1:200)
  colnames(new_data) <- paste0("InvalidSample", 1:20)
  expect_error(sn_layer_data(so, layer = "data") <- new_data, "Row or column names do not match.")
})

test_that("sn_metadata and sn_add_metadata work correctly", {
  data("so")

  # Test sn_metadata
  expect_equal(sn_metadata(so)$condition, so$condition)

  # Test sn_add_metadata with matching row names
  new_metadata <- rep("NewCondition", 20)
  names(new_metadata) <- colnames(so)
  so <- sn_add_metadata(so, metadata = new_metadata, col_name = "new_column")
  expect_equal(sn_metadata(so)$new_column, unname(new_metadata))

  # Test error when names do not match
  invalid_metadata <- rep("InvalidCondition", 20)
  names(invalid_metadata) <- paste0("InvalidSample", 1:20)
  expect_error(
    sn_add_metadata(so, metadata = invalid_metadata, col_name = "invalid_column"),
    "Names of metadata vector must include all object samples."
  )

  # Test error when col_name is missing
  expect_error(
    sn_add_metadata(so, metadata = new_metadata),
    "Must provide 'col_name' when adding a single vector as metadata."
  )
})

test_that("sn_active_ident handles drop and replace options", {
  data("so")
  sn_active_ident(so, replace = TRUE) <- factor(rep("x", 20))
  expect_equal(levels(sn_active_ident(so)), "x")

  sn_active_ident(so) <- factor(rep(c("a", "b"), 10), levels = c("a", "b", "c"))
  expect_equal(levels(sn_active_ident(so)), c("a", "b", "c"))
  sn_active_ident(so, drop = TRUE) <- sn_active_ident(so)
  expect_equal(levels(sn_active_ident(so)), c("a", "b"))
})

test_that("sn_version returns correct version string", {
  data("so")
  expect_type(sn_version(so), "character")
})

test_that("sn_layer_data<- supports NULL (layer removal)", {
  data("so")
  expect_true("scale.data" %in% sn_layers(so))
  sn_layer_data(so, layer = "scale.data") <- NULL
  expect_false("scale.data" %in% sn_layers(so))
})

test_that("dim and dimnames return expected dimensions", {
  data("so")
  d <- unname(dim(so))
  expect_equal(length(d), 2)
  expect_equal(d[2], nrow(sn_metadata(so)))

  dims <- dimnames(so)
  expect_equal(dims$samples, sn_samples(so))
})

test_that("head and tail print summaries", {
  data("so")
  expect_invisible(head(so, n = 3))
  expect_invisible(tail(so, n = 2))
})

test_that("sn_fetch_data fetches metadata, expression, ident, and reduction", {
  data("so")
  vars <- c("condition", "gender", "gene1", "PC_1", "ident")
  df <- sn_fetch_data(so, vars = vars)
  expect_true(all(c("condition", "gender", "gene1", "PC_1", "ident") %in% colnames(df)))
  expect_equal(nrow(df), unname(ncol(so)))
})


