# Test for sn_create_shennong_object

test_that("sn_create_shennong_object creates a valid Shennong object", {
  # Create a dummy counts matrix
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)

  # Create dummy metadata
  sample_metadata <- data.frame(
    condition = sample(c("Control", "Treated"), 10, replace = TRUE),
    batch = sample(c("A", "B"), 10, replace = TRUE),
    row.names = colnames(counts_matrix)
  )

  # Create a Shennong object
  sn_obj <- sn_create_shennong_object(
    counts = counts_matrix,
    metadata = sample_metadata,
    project = "TestProject",
    organism = "Mus musculus"
  )

  # Assertions
  expect_s4_class(sn_obj, "Shennong")
  expect_equal(as.numeric(ncol(sn_obj)), 10) # Convert to numeric to avoid name issues
  expect_equal(as.numeric(nrow(sn_obj)), 100) # Convert to numeric to avoid name issues
  expect_equal(sn_obj@project_name, "TestProject")
  expect_equal(sn_obj@organism, "Mus musculus")
})

test_that("sn_create_shennong_object filters samples correctly", {
  # Create a dummy counts matrix
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)

  # Set the first sample to have low counts
  counts_matrix[, 1] <- 0

  # Create dummy metadata matching filtered samples
  sample_metadata <- data.frame(
    condition = sample(c("Control", "Treated"), 9, replace = TRUE),
    batch = sample(c("A", "B"), 9, replace = TRUE),
    row.names = colnames(counts_matrix)[-1] # Exclude Sample1
  )

  # Create a Shennong object with filtering
  sn_obj <- sn_create_shennong_object(
    counts = counts_matrix,
    metadata = sample_metadata,
    min_counts = 10
  )

  # Assertions
  expect_equal(as.numeric(ncol(sn_obj)), 9) # One sample should be filtered out
  expect_false("Sample1" %in% colnames(sn_obj))
})

test_that("sn_create_shennong_object handles metadata correctly", {
  # Create a dummy counts matrix
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:100)
  colnames(counts_matrix) <- paste0("Sample", 1:10)

  # Create dummy metadata matching the counts matrix
  sample_metadata <- data.frame(
    condition = sample(c("Control", "Treated"), 10, replace = TRUE),
    batch = sample(c("A", "B"), 10, replace = TRUE),
    row.names = colnames(counts_matrix) # Match exactly
  )

  # Create a Shennong object
  sn_obj <- sn_create_shennong_object(
    counts = counts_matrix,
    metadata = sample_metadata
  )

  # Assertions
  expect_equal(as.numeric(ncol(sn_obj)), 10) # No samples should be filtered
  expect_equal(rownames(sn_obj@colData), colnames(counts_matrix)) # Metadata should match samples
})

test_that("sn_create_shennong_object throws errors for invalid input", {
  # Create a dummy counts matrix without column names
  counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)

  # Assertions
  expect_error(
    sn_create_shennong_object(counts = counts_matrix),
    "Column names are required for the counts matrix"
  )

  # Create a counts matrix with duplicate column names
  colnames(counts_matrix) <- c(rep("Sample1", 5), paste0("Sample", 6:10))
  expect_error(
    sn_create_shennong_object(counts = counts_matrix),
    "Duplicate sample names"
  )
})
