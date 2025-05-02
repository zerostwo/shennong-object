# Test for methods in methods-shennong.R

test_that("sn_assays returns correct assay names", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test
  expect_equal(sn_assays(mock_object), c("RNA"))
})

test_that("sn_active_assay gets and sets active assay", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test getter
  expect_equal(sn_active_assay(mock_object), "RNA")

  # Test setter
  sn_active_assay(mock_object) <- "RNA"
  expect_equal(sn_active_assay(mock_object), "RNA")

  # Test invalid assay
  expect_error(sn_active_assay(mock_object) <- "Invalid", "Assay 'Invalid' not found")
})

test_that("sn_layers returns correct layers", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test
  expect_equal(sn_layers(mock_object), c("counts"))
})

test_that("sn_active_layer gets and sets active layer", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test getter
  expect_equal(sn_active_layer(mock_object), "counts")

  # Test setter with valid layer
  sn_active_layer(mock_object) <- "counts"
  expect_equal(sn_active_layer(mock_object), "counts")

  # Test invalid layer
  expect_error(sn_active_layer(mock_object) <- "invalid_layer", "Layer 'invalid_layer' not found")
})

test_that("sn_samples and sn_features return correct names", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test
  expect_equal(sn_samples(mock_object), colnames(counts_matrix))
  expect_equal(sn_features(mock_object), rownames(counts_matrix))
})

test_that("sn_active_ident gets and sets active ident", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test getter
  expect_equal(unname(sn_active_ident(mock_object)), factor(rep("TestProject", 10)))

  # Test setter
  sn_active_ident(mock_object) <- factor(rep("NewIdent", 10))
  expect_equal(unname(sn_active_ident(mock_object)), factor(rep("NewIdent", 10)))
})

test_that("sn_project_name gets and sets project name", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test getter
  expect_equal(sn_project_name(mock_object), "TestProject")

  # Test setter
  sn_project_name(mock_object) <- "NewProject"
  expect_equal(sn_project_name(mock_object), "NewProject")
})

test_that("sn_layer_data gets and sets layer data", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Test setter with matching row and column names
  new_data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
  rownames(new_data) <- rownames(counts_matrix)
  colnames(new_data) <- colnames(counts_matrix)
  sn_layer_data(mock_object, layer = "counts") <- new_data
  expect_equal(sn_layer_data(mock_object, layer = "counts"), new_data)

  # Test error when row or column names do not match
  invalid_data <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
  rownames(invalid_data) <- paste0("InvalidGene", 1:100)
  colnames(invalid_data) <- paste0("InvalidSample", 1:10)
  expect_error(sn_layer_data(mock_object, layer = "counts") <- invalid_data, "Row or column names do not match.")
})

test_that("sn_metadata and sn_add_metadata work correctly", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  sample_metadata <- data.frame(
    condition = sample(c("Control", "Treated"), 10, replace = TRUE),
    batch = sample(c("A", "B"), 10, replace = TRUE),
    row.names = colnames(counts_matrix)
  )
  mock_object <- sn_create_shennong_object(counts = counts_matrix, metadata = sample_metadata, project = "TestProject")

  # Test sn_metadata
  expect_equal(sn_metadata(mock_object)$condition, sample_metadata$condition)

  # Test sn_add_metadata with matching row names
  new_metadata <- c("New1", "New2", "New3", "New4", "New5", "New6", "New7", "New8", "New9", "New10")
  names(new_metadata) <- colnames(counts_matrix)
  mock_object <- sn_add_metadata(mock_object, metadata = new_metadata, col_name = "new_column")
  expect_equal(sn_metadata(mock_object)$new_column, unname(new_metadata))

  # Test error when names do not match
  invalid_metadata <- c("Invalid1", "Invalid2", "Invalid3", "Invalid4", "Invalid5", "Invalid6", "Invalid7", "Invalid8", "Invalid9", "Invalid10")
  names(invalid_metadata) <- paste0("InvalidSample", 1:10)
  expect_error(
    sn_add_metadata(mock_object, metadata = invalid_metadata, col_name = "invalid_column"),
    "Names of metadata vector must include all object samples."
  )

  # Test error when col_name is missing
  expect_error(
    sn_add_metadata(mock_object, metadata = new_metadata),
    "Must provide 'col_name' when adding a single vector as metadata."
  )
})

test_that("sn_add_metadata handles vector input correctly", {
  # Create a mock Shennong object using sn_create_shennong_object
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)
  mock_object <- sn_create_shennong_object(counts = counts_matrix, project = "TestProject")

  # Create a metadata vector with matching names
  new_metadata <- c("New1", "New2", "New3", "New4", "New5", "New6", "New7", "New8", "New9", "New10")
  names(new_metadata) <- colnames(counts_matrix)

  # Add metadata vector
  mock_object <- sn_add_metadata(mock_object, metadata = new_metadata, col_name = "new_column")

  # Test
  expect_equal(sn_metadata(mock_object)$new_column, unname(new_metadata))

  # Test error when names do not match
  invalid_metadata <- c("Invalid1", "Invalid2", "Invalid3", "Invalid4", "Invalid5", "Invalid6", "Invalid7", "Invalid8", "Invalid9", "Invalid10")
  names(invalid_metadata) <- paste0("InvalidSample", 1:10)
  expect_error(
    sn_add_metadata(mock_object, metadata = invalid_metadata, col_name = "invalid_column"),
    "Names of metadata vector must include all object samples."
  )

  # Test error when col_name is missing
  expect_error(
    sn_add_metadata(mock_object, metadata = new_metadata),
    "Must provide 'col_name' when adding a single vector as metadata."
  )
})
