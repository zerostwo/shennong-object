#' @include generics.R
NULL

#' @rdname sn_reductions
#' @export
setMethod("sn_reductions", "Shennong", function(object, ...) {
  return(names(slot(object = object, name = "reductions")))
})


#' @rdname sn_active_reduction
#' @export
setMethod("sn_active_reduction", signature("Shennong", "ANY"), function(object, assay = NULL) {
  reds <- object@reductions
  hits <- vapply(reds, function(r) {
    isTRUE(r@global) && (is.null(assay) || r@assay_used == assay)
  }, logical(1))

  if (sum(hits) == 0) {
    return(NA_character_)
  } else if (sum(hits) > 1 && is.null(assay)) {
    warn("Multiple global reductions found. Specify `assay` to disambiguate.")
  }

  names(reds)[which(hits)][[1]]
})

#' @rdname sn_active_reduction
#' @export
setReplaceMethod(
  f = "sn_active_reduction",
  signature = c("Shennong", "character"),
  definition = function(object, value) {
    if (!value %in% names(object@reductions)) {
      abort(glue("Reduction '{value}' not found in object."))
    }

    red_to_set <- object@reductions[[value]]
    assay <- red_to_set@assay_used

    # # update all reductions matching this assay
    # object@reductions <- purrr::imap(object@reductions, function(r, nm) {
    #   if (r@assay_used == assay) {
    #     r@global <- (nm == value)
    #   }
    #   r
    # })
    # update all reductions matching this assay
    for (nm in names(object@reductions)) {
      r <- object@reductions[[nm]]
      if (r@assay_used == assay) {
        r@global <- (nm == value)
      }
      object@reductions[[nm]] <- r
    }

    validObject(object)
    object
  }
)

#' @rdname sn_embeddings
#' @export
setMethod("sn_embeddings", "ShennongReduction", function(object) {
  object@embedding
})

#' @rdname sn_loadings
#' @export
setMethod("sn_loadings", "ShennongReduction", function(object) {
  object@loadings
})

#' @rdname sn_stdev
#' @export
setMethod("sn_stdev", "ShennongReduction", function(object) {
  object@stdev
})

#' @describeIn Shennong Assign a reduction to the object by name
#' @param x A Shennong object
#' @param i Name of the reduction to assign
#' @param j (Ignored) Required for method signature but not used.
#' @param value A ShennongReduction object
#' @param ... Additional arguments (ignored)
#' @export
setMethod(
  f = "[[",
  signature = c("Shennong", "character"),
  definition = function(x, i, ...) {
    if (i %in% sn_assays(x)) {
      return(methods::callGeneric(as(x, "MultiAssayExperiment"), i, ...))
    }

    if (i %in% names(x@reductions)) {
      return(x@reductions[[i]])
    }

    abort(glue(
      "Object has neither assay nor reduction named '{i}'.\n",
      "Available assays: {paste0(sn_assays(x), collapse = ', ')}\n",
      "Available reductions: {paste0(names(x@reductions), collapse = ', ')}"
    ))
  }
)


#' @describeIn Shennong Assign a reduction to the object by name
#' @param value A ShennongReduction object
#' @export
setReplaceMethod(
  f = "[[",
  signature = c(x = "Shennong", i = "character", j = "missing", value = "ShennongReduction"),
  definition = function(x, i, j, ..., value) {
    if (!is.character(i) || length(i) != 1) {
      abort("Reduction name must be a single character string.")
    }

    # 1. Set this reduction's global = TRUE
    assay <- value@assay_used
    value@global <- TRUE

    # # 2. Clear global flags for existing reductions with the same assay
    # x@reductions <- purrr::imap(x@reductions, function(r, nm) {
    #   if (r@assay_used == assay) {
    #     r@global <- FALSE
    #   }
    #   r
    # })
    # Clear global flags for existing reductions with the same assay
    for (nm in names(x@reductions)) {
      r <- x@reductions[[nm]]
      if (r@assay_used == assay) {
        r@global <- FALSE
      }
      x@reductions[[nm]] <- r
    }

    # 3. Assign the new reduction
    x@reductions[[i]] <- value

    validObject(x)
    x
  }
)
