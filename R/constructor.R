# Shennong Object Constructor Function ------------------------------------
#' Create a Shennong Object
#'
#' Initializes a Shennong object from a primary counts matrix and optional
#' sample metadata.
#'
#' @param counts A matrix-like object (e.g., matrix, dgCMatrix) containing the
#'   raw counts. Rows typically represent features (genes) and columns represent
#'   samples. Column names are required and will be used as primary sample
#'   identifiers.
#' @param assay Character string. The name to assign to the primary assay
#'   created from the \code{counts} matrix (default: "RNA").
#' @param metadata Optional data frame or DataFrame containing sample-level
#'   metadata. Row names must match the column names of the \code{counts}
#'   matrix.
#' @param project Character string. Name for the project (default: "Shennong").
#' @param organism Optional character string specifying the organism.
#' @param min_counts Minimum number of counts required for a sample to be
#'   included. Samples with total counts below this threshold will be filtered
#'   out. Default is 0 (no filtering).
#' @param min_features Minimum number of detected features (genes with count >
#'   0) required for a sample to be included. Default is 0 (no filtering).
#' @param ... Additional arguments passed to the internal \code{new("Shennong",
#'   ...)} call.
#'
#' @return A \code{Shennong} object.
#'
#' @importFrom stats setNames
#' @export
#' @examples
#' # Create a dummy counts matrix
#' counts_matrix <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
#' rownames(counts_matrix) <- paste0("Gene", 1:100)
#' colnames(counts_matrix) <- paste0("Sample", 1:10)
#'
#' # Create dummy metadata
#' sample_metadata <- data.frame(
#'   condition = sample(c("Control", "Treated"), 10, replace = TRUE),
#'   batch = sample(c("A", "B"), 10, replace = TRUE),
#'   row.names = colnames(counts_matrix)
#' )
#'
#' # Create a Shennong object
#' sn_obj <- sn_create_shennong_object(
#'   counts = counts_matrix,
#'   metadata = sample_metadata,
#'   project = "TestProject",
#'   organism = "Mus musculus"
#' )
#'
#' print(sn_obj)
#'
#' # Create object without metadata
#' sn_obj_no_meta <- sn_create_shennong_object(counts = counts_matrix)
#' print(sn_obj_no_meta)
#'
#' # Apply basic filtering during creation
#' counts_matrix[, 1] <- 0 # Make first sample have 0 counts
#' sn_obj_filtered <- sn_create_shennong_object(
#'   counts = counts_matrix,
#'   metadata = sample_metadata,
#'   min_counts = 10, # Filter sample 1
#'   min_features = 5
#' )
#' print(sn_obj_filtered) # Should have 9 samples
#'
sn_create_shennong_object <- function(
    counts,
    assay = "RNA",
    metadata = NULL,
    project = "Shennong",
    organism = character(0),
    min_counts = 0,
    min_features = 0,
    ...) {
  if (missing(counts)) {
    stop("Input 'counts' matrix is required.")
  }
  if (!inherits(counts, "matrix") && !inherits(counts, "Matrix")) {
    stop("'counts' must be a matrix or Matrix object.")
  }
  if (is.null(colnames(counts))) {
    stop("Column names are required for the counts matrix (representing sample IDs).")
  }
  if (any(duplicated(colnames(counts)))) {
    stop("Duplicate sample names (colnames) are not allowed in the counts matrix.")
  }
  if (!is.character(assay) || length(assay) != 1 || nchar(assay) == 0) {
    stop("'assay' must be a non-empty character string.")
  }

  initial_samples <- colnames(counts)
  samples_to_keep <- rep(TRUE, ncol(counts))

  if (min_counts > 0) {
    sample_counts <- colSums(counts)
    samples_to_keep <- samples_to_keep & (sample_counts >= min_counts)
  }
  if (min_features > 0) {
    if (inherits(counts, "dgCMatrix")) {
      sample_features <- diff(counts@p)
    } else {
      sample_features <- colSums(counts > 0)
    }
    samples_to_keep <- samples_to_keep & (sample_features >= min_features)
  }

  # Ensure sample names are preserved after filtering
  if (!all(samples_to_keep)) {
    counts <- counts[, samples_to_keep, drop = FALSE]
    filtered_samples <- initial_samples[!samples_to_keep]
    current_samples <- colnames(counts) # Update current_samples to match filtered counts
    message(paste(
      "Filtered out", length(filtered_samples), "samples based on min_counts/min_features:",
      paste(filtered_samples, collapse = ", ")
    ))
    if (ncol(counts) == 0) {
      stop("No samples remaining after filtering.")
    }
  } else {
    current_samples <- initial_samples
  }

  # --- Compute QC Metrics: nFeature_<assay> and nCount_<assay> ---
  nCount <- colSums(counts)
  if (inherits(counts, "dgCMatrix")) {
    nFeature <- diff(counts@p)
  } else {
    nFeature <- colSums(counts > 0)
  }
  qc_metrics <- DataFrame(
    setNames(list(nFeature, nCount), c(paste0("nFeature_", assay), paste0("nCount_", assay)))
  )
  rownames(qc_metrics) <- current_samples

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata) && !inherits(metadata, "DataFrame")) {
      stop("'metadata' must be a data.frame or DataFrame.")
    }
    if (is.null(rownames(metadata))) {
      stop("Row names are required for the metadata (must match sample names).")
    }
    metadata_samples <- rownames(metadata)
    missing_meta <- setdiff(current_samples, metadata_samples)
    if (length(missing_meta) > 0) {
      warning("Metadata missing for samples: ", paste(missing_meta, collapse = ", "))
    }
    extra_meta <- setdiff(metadata_samples, current_samples)
    if (length(extra_meta) > 0) {
      warning(
        "Metadata provided for samples not present in (filtered) counts: ",
        paste(extra_meta, collapse = ", "), ". Removing extra metadata."
      )
    }
    metadata <- metadata[current_samples, , drop = FALSE]
    mae_coldata <- DataFrame(metadata)
    mae_coldata <- cbind(mae_coldata, qc_metrics)
  } else {
    mae_coldata <- qc_metrics
  }

  se <- SummarizedExperiment(
    assays = list(counts = counts)
  )
  metadata(se)$active_layer <- "counts"

  exp_list <- ExperimentList(setNames(list(se), assay))

  map_df <- DataFrame(
    primary = current_samples,
    colname = current_samples,
    assay = factor(rep(assay, length(current_samples)), levels = assay)
  )
  sample_map <- listToMap(setNames(list(map_df), assay))

  mae <- MultiAssayExperiment(
    experiments = exp_list,
    colData = mae_coldata,
    sampleMap = sample_map
  )

  # initial_ident <- factor(current_samples)
  initial_ident <- factor(rep(project, length(current_samples)),
    levels = unique(project)
  )
  names(initial_ident) <- current_samples

  sn_obj <- new(
    "Shennong",
    mae,
    organism = as.character(organism)[1],
    active_assay = assay,
    active_ident = initial_ident,
    reductions = list(),
    project_name = as.character(project)[1],
    misc = list(),
    version = as.character(utils::packageVersion("ShennongObject")),
    commands = list(),
    tools = list(),
    ...
  )

  return(sn_obj)
}


# Logging Function --------------------------------------------------------
#' Log a Shennong command
#'
#' Captures the current function call and parameters for provenance tracking.
#' Can return the command object or store it in the Shennong object's @commands
#' list.
#'
#' @param object A Shennong object
#' @param return_command Logical; return command object instead of storing it
#'
#' @return A Shennong object with updated @commands slot, or a ShennongCommand
#'   object
#' @export
sn_log_shennong_command <- function(object, return_command = FALSE) {
  time_stamp <- Sys.time()
  which.frame <- sys.nframe() - 1
  if (which.frame < 1) {
    stop("'LogShennongCommand' cannot be called at the top level", call. = FALSE)
  }

  if (as.character(sys.calls()[[1]])[1] == "do.call") {
    call_string <- deparse(sys.calls()[[1]])
    command_name <- as.character(sys.calls()[[1]])[2]
  } else {
    command_name <- deparse(sys.calls()[[which.frame]])
    command_name <- gsub("\\.Shennong", "", command_name)
    call_string <- command_name
    command_name <- sub("\\(.*$", "", command_name)
  }

  argnames <- names(formals(sys.function(which = sys.parent(1))))
  argnames <- setdiff(argnames, c("object", "..."))
  params <- list()
  p.env <- parent.frame(1)
  argnames <- intersect(argnames, ls(p.env))
  # for (arg in argnames) {
  #   param_value <- get(arg, envir = p.env)
  #   if (inherits(param_value, "Shennong")) next
  #   params[[arg]] <- param_value
  # }
  for (arg in argnames) {
    param_value <- get(arg, envir = p.env)
    if (inherits(param_value, "Shennong")) next

    # Skip large data objects for logging
    if (inherits(param_value, c("matrix", "Matrix", "SummarizedExperiment", "data.frame"))) {
      params[[arg]] <- glue("<{class(param_value)[1]}: {dim(param_value)[1]} x {dim(param_value)[2]}>")
    } else {
      params[[arg]] <- param_value
    }
  }

  assay <- params[["assay"]]
  cmd.assay <- assay
  if (!is.null(assay)) {
    command_name <- paste(command_name, assay, sep = ".")
  }
  command_name <- gsub("\\.+$", "", command_name)
  command_name <- sub("\\.\\.", ".", command_name)

  shen_cmd <- new(
    Class = "ShennongCommand",
    name = command_name,
    params = params,
    time_stamp = time_stamp,
    call_string = call_string,
    assay_used = cmd.assay %||% NA_character_
  )

  if (isTRUE(return_command)) {
    return(shen_cmd)
  }
  object@commands[[command_name]] <- shen_cmd
  return(object)
}


# Reduction ---------------------------------------------------------------
#' Create a dimensional reduction object for Shennong
#'
#' @param embedding A matrix of sample coordinates in reduced space.
#' @param loadings Optional matrix of feature loadings.
#' @param stdev Optional numeric vector of standard deviations.
#' @param assay_used Assay name from which the reduction is derived.
#' @param key Character prefix used in column names of the embedding matrix.
#' @param global Logical. If TRUE, this reduction persists even if the assay is dropped.
#' @param misc List of additional metadata.
#'
#' @return A `ShennongReduction` object.
#' @export
sn_create_reduction <- function(
    embedding,
    loadings = matrix(0, nrow = 0, ncol = 0),
    stdev = numeric(),
    assay_used = NULL,
    key = NULL,
    global = FALSE,
    misc = list()) {
  stopifnot(is.matrix(embedding))

  if (is.null(key)) {
    # Infer common prefix from colnames
    if (is.null(key)) {
      if (!is.null(colnames(embedding))) {
        key_prefixes <- gsub("[0-9]+$", "", colnames(embedding))
        key_prefixes <- unique(key_prefixes[nzchar(key_prefixes)])
        if (length(key_prefixes) == 1) {
          key <- paste0(key_prefixes, "_")
        } else {
          abort("Cannot infer `key` from column names. Please provide explicitly.")
        }
      } else {
        abort("`key` must be provided when `colnames(embedding)` are NULL.")
      }
    }
  } else {
    if (!grepl("_$", key)) {
      key <- paste0(key, "_")
    }
  }

  if (is.null(colnames(embedding)) ||
    !all(grepl(paste0("^", key, "[[:digit:]]+$"), colnames(embedding)))) {
    colnames(embedding) <- paste0(key, seq_len(ncol(embedding)))
  }

  if (ncol(loadings) > 0 && !all(colnames(loadings) == colnames(embedding))) {
    colnames(loadings) <- colnames(embedding)
  }

  new(
    Class = "ShennongReduction",
    embedding = embedding,
    loadings = loadings,
    stdev = stdev,
    assay_used = assay_used %||% "RNA",
    global = global,
    key = key,
    misc = misc
  )
}
