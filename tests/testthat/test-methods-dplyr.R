library(testthat)
library(ShennongObject)

test_that("filter.Shennong filters samples correctly", {
  # Load demo Shennong object
  data(so)

  # Test filtering by metadata
  filtered_obj <- filter(so, condition == "Normal")
  expect_true(all(sn_metadata(filtered_obj)$condition == "Normal"))

  # Test filtering by expression data
  filtered_obj <- filter(so, gene1 > 2)
  expect_true(all(sn_layer_data(filtered_obj)["gene1", ] > 1))

  # Test filtering with no matches
  expect_error(filter(so, condition == "Nonexistent"), "No samples matched the filter criteria.")
})

test_that("filter handles reductions correctly", {
  # Load demo Shennong object
  data(so)

  # Test filtering by reduction
  filtered_obj <- filter(so, PC_1 > 0)
  expect_true(all(sn_embeddings(filtered_obj[["pca"]])[, "PC_1"] > 0))
})
