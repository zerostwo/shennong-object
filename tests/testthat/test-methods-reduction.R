test_that("sn_reductions returns reduction names", {
  data(so)
  expect_true("pca" %in% sn_reductions(so))
  expect_true("mds" %in% sn_reductions(so))
})

test_that("sn_active_reduction returns correct reduction", {
  data(so)
  active <- sn_active_reduction(so)
  expect_true(active %in% sn_reductions(so))

  # Make sure returned reduction is global
  r <- so[[active]]
  expect_true(r@global)
})

test_that("sn_active_reduction<- sets active reduction correctly", {
  data(so)
  sn_active_reduction(so) <- "pca"
  expect_equal(sn_active_reduction(so), "pca")
  expect_false(so[["mds"]]@global)
  expect_true(so[["pca"]]@global)
})

test_that("[[ and [[<- work for reductions", {
  data(so)

  # Extract
  pca <- so[["pca"]]
  expect_s4_class(pca, "ShennongReduction")

  # Replace and check
  new_pca <- pca
  new_pca@embedding <- new_pca@embedding + 1
  so[["pca"]] <- new_pca
  expect_equal(so[["pca"]]@embedding, new_pca@embedding)

  # Check that assigning resets global flags
  expect_true(so[["pca"]]@global)
  expect_false(so[["mds"]]@global)
})

test_that("sn_embeddings, sn_loadings, sn_stdev access correctly", {
  data(so)
  pca <- so[["pca"]]

  expect_equal(sn_embeddings(pca), pca@embedding)
  expect_equal(sn_loadings(pca), pca@loadings)
  expect_equal(sn_stdev(pca), pca@stdev)
})
