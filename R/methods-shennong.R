#' @include zzz.R
#' @include generics.R
NULL

#' @rdname sn_add_assay
#' @export
setMethod(
  "sn_add_assay", "Shennong", function(object, data, assay, layer = "counts") {
    if (!is.character(assay) || length(assay) != 1) {
      abort("'assay' must be a single character string.")
    }

    if (assay %in% sn_assays(object)) {
      abort(glue("Assay '{assay}' already exists in the object."))
    }

    # Prepare SummarizedExperiment if needed
    if (inherits(data, "SummarizedExperiment")) {
      se <- data
    } else if (inherits(data, c("matrix", "Matrix"))) {
      if (is.null(colnames(data))) {
        abort("Column names are required for data (sample names).")
      }
      se <- SummarizedExperiment(assays = list2(!!layer := data))
      metadata(se)$active_layer <- layer
    } else {
      abort("`data` must be a matrix or a SummarizedExperiment.")
    }

    samples <- colnames(se)
    current_samples <- colnames(object)

    if (any(!samples %in% current_samples)) {
      missing <- setdiff(samples, current_samples)
      abort(glue("Samples not found in Shennong object: {paste(missing, collapse = ', ')}"))
    }

    # Update ExperimentList
    object <- c(object, setNames(list(se), assay))

    object <- sn_log_shennong_command(object)
    validObject(object)
    return(object)
  }
)

#' @rdname sn_assays
#' @export
setMethod("sn_assays", "Shennong", function(object, ...) {
  return(names(experiments(object)))
})

#' @rdname sn_active_assay
#' @export
setMethod("sn_active_assay", "Shennong", function(object, ...) {
  return(slot(object = object, name = "active_assay"))
})

#' @rdname sn_active_assay
#' @export
setMethod("sn_active_assay<-", "Shennong", function(object, ..., value) {
  if (!is_character(value) || length(value) != 1) {
    abort("`value` must be a single character string representing an assay name.")
  }
  if (!(value %in% sn_assays(object))) {
    abort(
      glue(
        "Assay '{value}' not found in the object. Available assays: {paste(sn_assays(object), collapse = ', ')}"
      )
    )
  }
  slot(object, "active_assay") <- value
  validObject(object)
  return(object)
})

#' @rdname sn_layers
#' @export
setMethod("sn_layers", "Shennong", function(object, assay = NULL, ...) {
  assay <- assay %||% sn_active_assay(object)
  if (!assay %in% sn_assays(object)) {
    abort(
      glue(
        "Assay '{assay}' not found. Available assays: {paste(sn_assays(object), collapse = ', ')}"
      )
    )
  }
  return(assayNames(experiments(object)[[assay]]))
})

#' @rdname sn_active_layer
#' @export
setMethod("sn_active_layer", "Shennong", function(object, assay = NULL, ...) {
  assay <- assay %||% sn_active_assay(object)
  se <- experiments(object)[[assay]]
  return(metadata(se)$active_layer %||% NULL)
})

#' @rdname sn_active_layer
#' @export
setMethod("sn_active_layer<-", "Shennong", function(object, assay = NULL, ..., value) {
  if (!is_character(value) || length(value) != 1) {
    abort("`value` must be a single character string.")
  }

  assay <- assay %||% sn_active_assay(object)
  available_layers <- sn_layers(object, assay = assay)
  if (!value %in% available_layers) {
    abort(
      glue(
        "Layer '{value}' not found in assay '{assay}'. Available layers: {paste(available_layers, collapse = ', ')}"
      )
    )
  }
  se <- experiments(object)[[assay]]
  metadata(se)$active_layer <- value
  experiments(object)[[assay]] <- se
  validObject(object)
  return(object)
})

#' @rdname sn_samples
#' @export
setMethod("sn_samples", "Shennong", function(object, assay = NULL, ...) {
  assay <- assay[1L] %||% sn_active_assay(object)
  colnames(experiments(object)[[assay]])
})

#' @rdname sn_features
#' @export
setMethod("sn_features", "Shennong", function(object, assay = NULL, ...) {
  assay <- assay[1L] %||% sn_active_assay(object)
  rownames(experiments(object)[[assay]])
})

#' @rdname sn_active_ident
#' @export
setMethod("sn_active_ident", "Shennong", function(object, ...) {
  return(slot(object = object, name = "active_ident"))
})

#' @rdname sn_active_ident
#' @export
setMethod("sn_active_ident<-", "Shennong", function(object,
                                                    samples = NULL,
                                                    drop = FALSE,
                                                    replace = FALSE,
                                                    ...,
                                                    value) {
  # Argument check
  if (!is.factor(value) &&
    !is.atomic(value) && !is.character(value)) {
    abort("'value' must be a factor, atomic vector, or column name",
      class = "shennong_error_invalid_ident"
    )
  }

  all_samples <- colnames(object)
  meta <- colData(object)

  # 1. Parse samples
  samples <- samples %||% names(value) %||% all_samples
  if (is.numeric(samples)) {
    samples <- all_samples[samples]
  }
  samples <- intersect(samples, all_samples)

  if (length(samples) == 0) {
    rlang::warn("No matching samples found in 'samples'")
    return(object)
  }

  # 2. Parse value
  if (is.character(value) &&
    length(value) == 1 && value %in% colnames(meta)) {
    value <- meta[[value]]
    names(value) <- rownames(meta)
    value <- value[samples]
  } else if (is.list(value)) {
    value <- unlist(value, use.names = FALSE)
  }

  # 3. Construct ident
  idents <- if (isTRUE(replace)) {
    rep(NA_character_, length(all_samples))
  } else {
    as.character(object@active_ident)
  }
  names(idents) <- all_samples

  idents[samples] <- rep_len(as.character(value), length.out = length(samples))
  idents[is.na(idents)] <- "NA"

  # Preserve factor levels if value is a factor
  if (is.factor(value)) {
    levels <- levels(value)
  } else {
    levels <- unique(na.omit(idents))
  }
  idents <- factor(idents, levels = levels)

  if (isTRUE(drop)) {
    idents <- droplevels(idents)
  }

  object@active_ident <- idents
  validObject(object)
  return(object)
})

#' @rdname sn_levels
#' @export
setMethod("sn_levels", "Shennong", function(object, ...) {
  return(levels(sn_active_ident(object)))
})

#' @rdname sn_levels
#' @export
setMethod("sn_levels<-", "Shennong", function(object, ..., value) {
  idents <- sn_active_ident(object)
  if (!all(levels(idents) %in% value)) {
    abort("NA's would be generated by missing levels", class = "shennong_error_invalid_levels")
  }
  sn_active_ident(object) <- factor(idents, levels = value)
  return(object)
})

#' @rdname sn_project_name
#' @export
setMethod("sn_project_name", "Shennong", function(object, ...) {
  return(slot(object = object, name = "project_name"))
})

#' @rdname sn_project_name
#' @export
setMethod("sn_project_name<-", "Shennong", function(object, ..., value) {
  if (!is.character(value) || length(value) != 1) {
    abort("`value` must be a single character string representing the project name.")
  }
  slot(object, "project_name") <- value
  validObject(object)
  return(object)
})

#' @rdname sn_layer_data
#' @export
setMethod("sn_layer_data", "Shennong", function(object,
                                                layer = NULL,
                                                assay = NULL,
                                                ...) {
  assay <- assay %||% sn_active_assay(object)
  layer <- layer %||% sn_active_layer(object)

  available_assays <- sn_assays(object)
  if (!assay %in% available_assays) {
    abort(
      message = glue(
        "Assay '{assay}' not found in object.\nAvailable assays: {paste(available_assays, collapse = ', ')}"
      ),
      class = "shennong_error_assay_not_found"
    )
  }

  available_layers <- sn_layers(object = object, assay = assay)
  if (!layer %in% available_layers) {
    abort(
      message = glue(
        "Layer '{layer}' not found in assay '{assay}'.\nAvailable layers: {paste(available_layers, collapse = ', ')}"
      ),
      class = "shennong_error_layer_not_found"
    )
  }

  return(assay(experiments(object)[[assay]], layer))
})

#' @rdname sn_layer_data
#' @export
setMethod("sn_layer_data<-", "Shennong", function(object,
                                                  layer = NULL,
                                                  assay = NULL,
                                                  ...,
                                                  value) {
  assay <- assay %||% sn_active_assay(object)
  layer <- layer %||% sn_active_layer(object)

  available_assays <- sn_assays(object)
  if (!assay %in% available_assays) {
    abort(
      message = glue(
        "Assay '{assay}' not found in object.\nAvailable assays: {paste(available_assays, collapse = ', ')}"
      ),
      class = "shennong_error_assay_not_found"
    )
  }

  se <- experiments(object)[[assay]]

  # --- Handle NULL case (remove the layer) ---
  if (is.null(value)) {
    if (layer %in% assayNames(se)) {
      assay(se, layer) <- NULL
      message(glue("Removed layer '{layer}' from assay '{assay}'."))
    } else {
      warning(glue(
        "Layer '{layer}' not found in assay '{assay}', nothing to remove."
      ))
    }
    experiments(object)[[assay]] <- se
    return(object)
  }

  idx_row <- match(rownames(se), rownames(value))
  idx_col <- match(colnames(se), colnames(value))

  if (anyNA(idx_row) || anyNA(idx_col)) {
    abort("Row or column names do not match.")
  }
  value <- value[idx_row, idx_col]
  # --- Assign layer ---
  assay(se, layer) <- value
  experiments(object)[[assay]] <- se

  return(object)
})

#' @rdname sn_metadata
setMethod("sn_metadata", "Shennong", function(object, cols = NULL, drop = FALSE, ...) {
  cd <- colData(object) |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tibble::as_tibble()
  if (is.null(cols)) {
    return(cd)
  } else {
    if (!all(cols %in% colnames(cd))) {
      missing_cols <- setdiff(cols, colnames(cd))
      stop("Metadata columns not found: ", paste(missing_cols, collapse = ", "))
    }
    selected_cd <- cd[, cols, drop = drop]
    return(selected_cd)
  }
})


#' @rdname sn_add_metadata
#' @export
setMethod("sn_add_metadata", "Shennong", function(object, metadata, col_name = NULL, ...) {
  current_cd <- colData(object)
  current_samples <- rownames(current_cd)

  if (is.atomic(metadata) || is.factor(metadata)) {
    # Handle adding a single vector
    if (is.null(col_name)) {
      abort("Must provide 'col_name' when adding a single vector as metadata.")
    }
    if (length(metadata) != length(current_samples)) {
      abort(
        glue(
          "Length of metadata vector must match the number of samples ({length(current_samples)})."
        )
      )
    }
    if (!is.null(names(metadata))) {
      # Reorder if named
      if (!all(current_samples %in% names(metadata))) {
        abort("Names of metadata vector must include all object samples.")
      }
      metadata <- metadata[current_samples]
    }
    new_df <- DataFrame(setNames(list(metadata), col_name), row.names = current_samples)
  } else if (inherits(metadata, "DataFrame") ||
    is.data.frame(metadata)) {
    # Handle adding a DataFrame/data.frame
    if (is.null(rownames(metadata))) {
      stop("Input 'metadata' must have rownames matching object samples.")
    }
    if (!all(current_samples %in% rownames(metadata))) {
      missing_s <- setdiff(current_samples, rownames(metadata))
      stop(
        "Metadata is missing for samples: ",
        paste(missing_s, collapse = ", ")
      )
    }
    # Subset and reorder metadata to match object samples
    metadata <- metadata[current_samples, , drop = FALSE]
    new_df <- DataFrame(metadata)
    if (!is.null(col_name)) {
      warning("'col_name' ignored when adding a DataFrame/data.frame.")
    }
  } else {
    stop("'metadata' must be a vector, factor, data.frame, or DataFrame.")
  }

  # Check for collisions and combine
  cols_to_add <- colnames(new_df)
  cols_collision <- intersect(cols_to_add, colnames(current_cd))
  if (length(cols_collision) > 0) {
    warning(
      "Overwriting existing metadata columns: ",
      paste(cols_collision, collapse = ", ")
    )
    current_cd[, cols_collision] <- new_df[, cols_collision] # Overwrite
    cols_new <- setdiff(cols_to_add, cols_collision)
    if (length(cols_new) > 0) {
      current_cd <- cbind(current_cd, new_df[, cols_new, drop = FALSE]) # Add new ones
    }
  } else {
    current_cd <- cbind(current_cd, new_df) # Add all new ones
  }

  colData(object) <- current_cd

  object <- sn_log_shennong_command(object)

  validObject(object)
  return(object)
})

#' @rdname sn_version
#' @export
setMethod("sn_version", "Shennong", function(object, ...) {
  return(slot(object = object, name = "version"))
})

# Methods for R-defined generics ------------------------------------------
#' Get Dimensions of Shennong Object
#'
#' Returns the number of features in the active assay and the number of primary samples.
#'
#' @param x A Shennong object.
#'
#' @return A numeric vector of length 2: c(n_features, n_samples).
#' @export
#' @method dim Shennong
#' @rdname sn_dim
setMethod("dim", "Shennong", function(x) {
  n_samples <- nrow(colData(x))
  active_assay_name <- sn_active_assay(x)
  if (length(active_assay_name) > 0 && active_assay_name %in% sn_assays(x)) {
    n_features <- nrow(experiments(x)[[active_assay_name]])
  } else {
    warning("Active assay not set or not found, returning 0 features.")
    n_features <- 0L
  }
  return(c(Features = n_features, Samples = n_samples))
})

#' Get Dimension Names of Shennong Object
#'
#' Returns the feature names from the active assay and the primary sample names.
#'
#' @param x A Shennong object.
#'
#' @return A list of length 2 containing:
#'   1. Feature names (from active assay)
#'   2. Sample names (primary samples)
#' @export
#' @method dimnames Shennong
#' @rdname sn_dimnames
setMethod("dimnames", "Shennong", function(x) {
  sample_names <- sn_samples(x) # Primary sample names
  active_assay_name <- sn_active_assay(x)
  if (length(active_assay_name) > 0 && active_assay_name %in% sn_assays(x)) {
    feature_names <- sn_features(x, assay = active_assay_name)
  } else {
    warning("Active assay not set or not found, returning NULL for features.")
    feature_names <- NULL
  }
  return(list(featuers = feature_names, samples = sample_names))
})

#' @title Head/Tail of Shennong object
#'
#' @description Show the top or bottom rows of `colData()` in a `Shennong` object.
#' Useful for quick inspection of sample metadata.
#'
#' @param x A \code{Shennong} object.
#' @param n Number of rows to show. Default is 6.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input \code{Shennong} object.
#'
#' @rdname sn_head_tail
#' @aliases head.Shennong tail.Shennong
#' @export
setMethod("head", "Shennong", function(x, n = 6L, ...) {
  n <- as.integer(n)
  if (n <= 0) {
    warning("n must be a positive integer; returning input object.")
    return(invisible(x))
  }
  cat(sprintf(
    "Head of Shennong object (samples x features: %d x %d):\n",
    nrow(colData(x)), ncol(colData(x))
  ))
  print(utils::head(colData(x), n = n))
  invisible(x)
})

#' @rdname sn_head_tail
#' @export
setMethod("tail", "Shennong", function(x, n = 6L, ...) {
  n <- as.integer(n)
  if (n <= 0) {
    warning("n must be a positive integer; returning input object.")
    return(invisible(x))
  }
  cat(sprintf(
    "Tail of Shennong object (samples x features: %d x %d):\n",
    nrow(colData(x)), ncol(colData(x))
  ))
  print(utils::tail(colData(x), n = n))
  invisible(x)
})

#' @rdname sn_fetch_data
setMethod("sn_fetch_data", "Shennong", function(object, vars, samples = NULL, layer = NULL, clean = TRUE, ...) {
  all_samples <- colnames(object)
  samples <- samples %||% all_samples
  if (is.numeric(samples)) {
    samples <- all_samples[samples]
  }
  samples <- intersect(samples, all_samples)
  if (length(samples) == 0) {
    stop("No valid samples found.")
  }

  df <- data.frame(row.names = samples)

  # Try colData first
  meta_vars <- intersect(vars, colnames(colData(object)))
  if (length(meta_vars)) {
    df[meta_vars] <- as.data.frame(colData(object)[samples, meta_vars, drop = FALSE])
  }

  # Add 'ident' if requested
  if ("ident" %in% vars) {
    df$ident <- sn_active_ident(object)[samples]
  }

  # Try from active assay
  expr_vars <- setdiff(vars, c(meta_vars, "ident"))
  if (length(expr_vars)) {
    assay <- sn_active_assay(object)
    # assay_mat <- assays(object)[[assay_name]]
    assay_mat <- sn_layer_data(object, layer = layer, assay = assay)
    found_vars <- intersect(expr_vars, rownames(assay_mat))
    if (length(found_vars)) {
      expr_data <- t(assay_mat[found_vars, samples, drop = FALSE])
      df[found_vars] <- expr_data
    }
  }

  # Clean rows with all NAs
  if (isTRUE(clean)) {
    df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
  }

  return(df)
})
