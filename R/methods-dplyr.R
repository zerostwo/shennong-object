#' Filter samples in a Shennong object using dplyr syntax
#'
#' Allows tidyverse-style filtering of sample metadata (`colData`) in a Shennong object.
#' This behaves like `dplyr::filter()`, and subsets the object to retain only the samples
#' matching the specified conditions.
#'
#' @param .data A `Shennong` object.
#' @param ... Filtering conditions passed to `dplyr::filter()`.
#' @param layer The name of the expression matrix layer to use for filtering. If `NULL`,
#'   the active layer will be used.
#' @param .by Optional grouping variable(s).
#' @param .preserve Whether to preserve the grouping structure.
#'
#' @return A filtered `Shennong` object.
#'
#' @examples
#' \dontrun{
#' filter(so, stage == "Stage I", gender == "Female")
#' }
#'
#' @name filter.Shennong
#' @method filter Shennong
#' @importFrom dplyr filter mutate
#' @importFrom tibble as_tibble
#' @export
filter.Shennong <- function(.data, ..., layer = NULL, .by = NULL, .preserve = FALSE) {
  if (!inherits(.data, "Shennong")) {
    abort("The object must be a Shennong object.")
  }

  # Capture filtering expressions
  subset_exprs <- rlang::enquos(...)

  # 1. Get colData
  meta <- as_tibble(colData(.data), rownames = "observation")

  # 2. Get expression matrix (active layer)
  active_layer <- layer %||% sn_active_layer(.data)
  expr_before <- tryCatch(
    {
      mat <- sn_layer_data(.data, layer = active_layer)
      df <- as_tibble(t(mat), rownames = "observation")
      df
    },
    error = function(e) NULL
  )

  # 3. Get all reductions without prefixing column names
  reductions <- lapply(.data@reductions, function(red) {
    emb <- red@embedding
    df <- as.data.frame(emb)
    df$observation <- rownames(df)
    df
  })
  reductions <- Reduce(function(x, y) merge(x, y, by = "observation", all = TRUE), reductions)

  # 4. Merge all into one tibble
  dfs <- list(meta, expr_before, reductions)
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  df <- Reduce(function(x, y) merge(x, y, by = "observation", all = TRUE), dfs)
  df <- tibble::as_tibble(df)

  # 5. Filter
  df_filtered <- dplyr::filter(df, !!!subset_exprs, .by = {{ .by }}, .preserve = .preserve)
  if (nrow(df_filtered) == 0) {
    abort("No samples matched the filter criteria.")
  }

  kept_samples <- df_filtered$observation

  # 6. Subset Shennong object
  n_samples_before <- ncol(.data)
  n_features_before <- if (!is.null(expr_before)) ncol(expr_before) - 1 else NA

  .data <- .data[, kept_samples]
  .data@active_ident <- droplevels(.data@active_ident[kept_samples])

  # Filter reductions
  .data@reductions <- lapply(.data@reductions, function(red) {
    red@embedding <- red@embedding[kept_samples, , drop = FALSE]
    red
  })

  # 7. After filter, record for logging
  expr_after <- tryCatch(
    {
      sn_layer_data(.data, layer = active_layer)
    },
    error = function(e) NULL
  )

  n_samples_after <- length(kept_samples)
  n_features_after <- if (!is.null(expr_after)) nrow(expr_after) else NA

  log_info <- list(
    layer = active_layer,
    assay = sn_active_assay(.data),
    samples_before = n_samples_before,
    samples_after = n_samples_after,
    features_before = n_features_before,
    features_after = n_features_after
  )

  .data <- sn_log_shennong_command(.data, .info = log_info)
  return(.data)
}
