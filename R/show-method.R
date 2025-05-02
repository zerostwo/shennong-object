# Shennong-class ----------------------------------------------------------
#' Show method for Shennong objects
#'
#' Displays summary information about the Shennong object.
#'
#' @param object A Shennong object.
#' @importFrom methods callNextMethod slotNames slot
#' @export
#' @rdname Shennong-class
#' @aliases show,Shennong-method
setMethod("show", "Shennong", function(object) {
  cat("A Shennong object: ", object@project_name, "\n")
  cat("Organism:", ifelse(length(object@organism) > 0, object@organism, "Not set"), "\n")

  # Use callNextMethod to show the underlying MultiAssayExperiment information
  callNextMethod(object)

  # Add Shennong specific information
  cat("\n--- Shennong Specific Slots ---\n")
  cat("Active Assay:", ifelse(length(object@active_assay) > 0, object@active_assay, "Not set"), "\n")
  cat("Active Ident:", ifelse(length(object@active_ident) > 0,
    paste(
      nlevels(object@active_ident), "levels (first few:",
      paste(head(levels(object@active_ident)), collapse = ", "), ")"
    ),
    "Not set"
  ), "\n")
  cat(
    length(object@reductions), "dimensional reduction(s) present:",
    ifelse(length(object@reductions) > 0, paste(names(object@reductions), collapse = ", "), "None"), "\n"
  )

  # Optionally display misc, commands, tools summaries if needed
  cat("Misc data:", ifelse(length(object@misc) > 0, paste(names(object@misc), collapse = ", "), "None"), "\n")
  cat("Command history stored:", ifelse(length(object@commands) > 0, "Yes", "No"), "\n")
  cat("Tools data stored:", ifelse(length(object@tools) > 0, paste(names(object@tools), collapse = ", "), "None"), "\n")
  cat("Object created with ShennongObject version:", as.character(object@version), "\n")
})

# Command-class -----------------------------------------------------------
#' Show method for ShennongCommand
#'
#' Displays key information for a ShennongCommand object.
#'
#' @param object A ShennongCommand object
#' @return Invisibly returns the object.
#' @export
#' @rdname show.ShennongCommand
#' @aliases show,ShennongCommand-method
setMethod(
  f = "show", "ShennongCommand", function(object) {
    params <- slot(object, "params")
    params <- params[sapply(params, class) != "function"]
    cat("Command:", slot(object, "call_string"), "\n")
    cat("Time:", as.character(slot(object, "time_stamp")), "\n")
    for (p in seq_along(params)) {
      cat(names(params[p]), ":", params[[p]], "\n")
    }
  }
)

##' Show method for ShennongReduction
#'
#' @param object A ShennongReduction object.
#' @export
setMethod("show", "ShennongReduction", function(object) {
  n_samples <- nrow(object@embedding)
  n_dims <- ncol(object@embedding)
  key <- object@key
  assay <- object@assay_used
  has_loadings <- if (ncol(object@loadings) > 0) "present" else "absent"
  stdev_len <- length(object@stdev)
  is_global <- object@global
  misc_len <- length(object@misc)

  cat(glue::glue(
    "A ShennongReduction\n",
    "Assay used: {assay}\n",
    "Dimensions: {n_samples} samples x {n_dims} components (key = '{key}')\n",
    "Loadings: {has_loadings}\n",
    "Standard deviation: {stdev_len} dimension{ifelse(stdev_len == 1, '', 's')}\n",
    "Global: {is_global}\n",
    "Metadata: {misc_len} item{ifelse(misc_len == 1, '', 's')} in `misc`\n"
  ))
})
