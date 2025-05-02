#' @include zzz.R
NULL

# Generics for Shennong-class ---------------------------------------------
#' Get Sample Names
#' @param object A Shennong object or other class for which methods are defined.
#' @param assay Optional: Assay context.
#' @param ... Additional arguments.
#' @rdname sn_samples
#' @export
setGeneric("sn_samples", function(object, assay = NULL, ...) standardGeneric("sn_samples"))

#' Get Feature Names
#' @param object A Shennong object or other class.
#' @param assay Optional: Assay context.
#' @param ... Additional arguments.
#' @rdname sn_features
#' @export
setGeneric("sn_features", function(object, assay = NULL, ...) standardGeneric("sn_features"))

#' Get Assay Names
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_assays
#' @export
setGeneric("sn_assays", function(object, ...) standardGeneric("sn_assays"))

#' Get Default Assay
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_active_assay
#' @export
setGeneric("sn_active_assay", function(object, ...) standardGeneric("sn_active_assay"))

#' Set Default Assay<-
#' @param object A Shennong object.
#' @param ... Additional arguments.
#' @param value The assay name to set as active.
#' @rdname sn_active_assay
#' @export
setGeneric("sn_active_assay<-", function(object, ..., value) standardGeneric("sn_active_assay<-"))

#' Get Sample Identities
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_active_ident
#' @export
setGeneric("sn_active_ident", function(object, ...) standardGeneric("sn_active_ident"))

#' Set Sample Identities<-
#' @param object A Shennong object.
#' @param samples Optional: Samples to set identities for.
#' @param drop Logical; drop unused identities.
#' @param replace Logical; replace existing identities.
#' @param ... Additional arguments.
#' @param value Identities to set (factor, vector, or metadata column name).
#' @rdname sn_active_ident
#' @export
setGeneric("sn_active_ident<-", function(object, ..., value) standardGeneric("sn_active_ident<-"))

#' Get Identity Levels
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_levels
#' @export
setGeneric("sn_levels", function(object, ...) standardGeneric("sn_levels"))

#' Set Identity Levels<-
#' @param object A Shennong object.
#' @param ... Additional arguments.
#' @param value New levels (character vector).
#' @rdname sn_levels
#' @export
setGeneric("sn_levels<-", function(object, ..., value) standardGeneric("sn_levels<-"))

#' Get Project Name
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_project_name
#' @export
setGeneric("sn_project_name", function(object, ...) standardGeneric("sn_project_name"))

#' Set Project Name<-
#' @param object A Shennong object.
#' @param ... Additional arguments.
#' @param value New project name (character string).
#' @rdname sn_project_name
#' @export
setGeneric("sn_project_name<-", function(object, ..., value) standardGeneric("sn_project_name<-"))

#' Get Layer Data
#' @param object A Shennong object or other class.
#' @param assay Optional: Assay name.
#' @param layer Layer name or index.
#' @param ... Additional arguments.
#' @rdname sn_layer_data
#' @export
setGeneric("sn_layer_data", function(object, layer = NULL, assay = NULL, ...) standardGeneric("sn_layer_data"))

#' Set Layer Data<-
#' @param object A Shennong object.
#' @param assay Optional: Assay name.
#' @param layer Layer name or index.
#' @param ... Additional arguments.
#' @param value The data matrix to assign.
#' @rdname sn_layer_data
#' @export
setGeneric("sn_layer_data<-", function(object, layer, assay = NULL, ..., value) standardGeneric("sn_layer_data<-"))

#' Get Layer Names
#' @param object A Shennong object or other class.
#' @param assay Optional: Assay name.
#' @param ... Additional arguments.
#' @rdname sn_layers
#' @export
setGeneric("sn_layers", function(object, assay = NULL, ...) standardGeneric("sn_layers"))

#' Get the active layer
#'
#' @param object A Shennong object.
#' @param ... Ignored.
#'
#' @return A character scalar representing the current active layer.
#' @param assay Optional: Assay name.
#' @rdname sn_active_layer
#' @export
setGeneric("sn_active_layer", function(object, assay = NULL, ...) standardGeneric("sn_active_layer"))

#' Set the active layer
#'
#' @param object A Shennong object.
#' @param assay Optional: Assay name.
#' @param value A character string of the layer name.
#' @param ... Ignored.
#'
#' @return The modified Shennong object.
#' @rdname sn_active_layer
#' @export
setGeneric("sn_active_layer<-", function(object, assay = NULL, ..., value) standardGeneric("sn_active_layer<-"))

#' Get Sample Metadata
#' @param object A Shennong object or other class.
#' @param cols Optional: Columns to retrieve.
#' @param drop Logical: Drop dimensions if possible.
#' @param ... Additional arguments.
#' @rdname sn_metadata
#' @export
setGeneric("sn_metadata", function(object, cols = NULL, drop = FALSE, ...) standardGeneric("sn_metadata"))

#' Add Sample Metadata
#'
#' Adds new columns to the main sample metadata table (`colData`).
#'
#' @param object A Shennong object.
#' @param metadata A data frame or DataFrame containing the new metadata columns.
#'   Row names must match primary sample names in the object.
#' @param col_name Optional: If `metadata` is a single vector, specify the name
#'   for the new metadata column.
#' @param ... Ignored.
#'
#' @return The modified Shennong object.
#' @export
#' @rdname sn_add_metadata
#' @aliases AddMetaData
setGeneric("sn_add_metadata", function(object, metadata, col_name = NULL, ...) standardGeneric("sn_add_metadata"))

#' Add an assay to a Shennong object
#'
#' @param object A Shennong object.
#' @param data A matrix-like object or SummarizedExperiment.
#' @param assay Character, name of the new assay.
#' @param layer Optional. Name of the data layer (default = "counts").
#'
#' @return A Shennong object with the new assay added.
#' @export
setGeneric("sn_add_assay", function(object, data, assay, layer = "counts") {
  standardGeneric("sn_add_assay")
})

#' Fetch data from a Shennong object
#'
#' Retrieves metadata, identities, or assay-level expression data for specified variables and cells.
#'
#' @param object A Shennong object.
#' @param vars Character vector of variable names to fetch (columns in `colData`, or rownames in assays).
#' @param samples Optional character/numeric vector of sample names or indices.
#' @param layer Optional layer name if expression matrices have multiple versions (e.g., normalized).
#' @param clean Logical or character. If TRUE, remove samples with all missing vars; if "ident", keep only rows with valid identity.
#' @param ... Not used.
#'
#' @return A data.frame with samples as rows and requested variables as columns.
#' @export
setGeneric("sn_fetch_data", function(object, vars, samples = NULL, layer = NULL, clean = TRUE, ...) standardGeneric("sn_fetch_data"))

#' Get samples grouped by identity class
#' @param object A Shennong object.
#' @param idents Identity levels to include.
#' @param samples Samples to consider.
#' @param return.null Return NULL if no samples found.
#' @rdname sn_samples_by_identities
#' @export
setGeneric("sn_samples_by_identities", function(object, idents = NULL, samples = NULL, return.null = FALSE) standardGeneric("sn_samples_by_identities"))

#' Reorder Sample Identities
#' @param object A Shennong object.
#' @param var Variable to use for ordering.
#' @param reverse Reverse order.
#' @param fun Aggregation function.
#' @param reorder.numeric Rename levels numerically.
#' @param ... Additional arguments.
#' @rdname sn_reorder_ident
#' @export
setGeneric("sn_reorder_ident", function(object, var, reverse = FALSE, fun = mean, reorder.numeric = FALSE, ...) standardGeneric("sn_reorder_ident"))

#' Stash Current Identities
#' @param object A Shennong object.
#' @param save.name Column name to save identities under.
#' @param ... Additional arguments.
#' @rdname sn_stash_ident
#' @export
setGeneric("sn_stash_ident", function(object, save.name = "orig.ident", ...) standardGeneric("sn_stash_ident"))

# Generics for Command-class ----------------------------------------------
#' Get Command History
#' @param object A Shennong object or other class.
#' @param command Optional: Specific command name.
#' @param ... Additional arguments.
#' @rdname sn_commands
#' @export
setGeneric("sn_commands", function(object, command = NULL, ...) standardGeneric("sn_commands"))

#' Set Command History<-
#' @param object A Shennong object.
#' @param command Name for the command log entry.
#' @param ... Additional arguments.
#' @param value The command log data.
#' @rdname sn_commands
#' @export
setGeneric("sn_commands<-", function(object, command, ..., value) standardGeneric("sn_commands<-"))

# Generics for Reduction-class --------------------------------------------
#' Get Dimensional Reductions
#' @param object A Shennong object or other class.
#' @param reduction Optional: Specific reduction name.
#' @param ... Additional arguments.
#' @rdname sn_reductions
#' @export
setGeneric("sn_reductions", function(object, reduction = NULL, ...) standardGeneric("sn_reductions"))

#' Set Dimensional Reduction<-
#' @param object A Shennong object.
#' @param reduction Name for the reduction.
#' @param ... Additional arguments.
#' @param value The reduction data (list).
#' @rdname sn_reductions
#' @export
setGeneric("sn_reductions<-", function(object, reduction, ..., value) standardGeneric("sn_reductions<-"))

#' Get Embeddings
#' @param object A Shennong object or other class.
#' @param reduction The name of the reduction.
#' @param ... Additional arguments.
#' @rdname sn_embeddings
#' @export
setGeneric("sn_embeddings", function(object, reduction, ...) standardGeneric("sn_embeddings"))

#' Get Feature Loadings
#' @param object A Shennong object or other class.
#' @param reduction The name of the reduction.
#' @param ... Additional arguments.
#' @rdname sn_loadings
#' @export
setGeneric("sn_loadings", function(object, reduction, ...) standardGeneric("sn_loadings"))

#' Get Standard Deviations
#' @param object A Shennong object or other class.
#' @param reduction The name of the reduction.
#' @param ... Additional arguments.
#' @rdname sn_stdev
#' @export
setGeneric("sn_stdev", function(object, reduction, ...) standardGeneric("sn_stdev"))

# Generics for Tools-class ------------------------------------------------
#' Get Tool Data
#' @param object A Shennong object or other class.
#' @param tool Optional: Specific tool name.
#' @param ... Additional arguments.
#' @rdname sn_tools
#' @export
setGeneric("sn_tools", function(object, tool = NULL, ...) standardGeneric("sn_tools"))

#' Set Tool Data<-
#' @param object A Shennong object.
#' @param tool Name of the tool.
#' @param ... Additional arguments.
#' @param value The data generated by the tool.
#' @rdname sn_tools
#' @export
setGeneric("sn_tools<-", function(object, tool, ..., value) standardGeneric("sn_tools<-"))

# Others ------------------------------------------------------------------
#' Get Miscellaneous Data
#' @param object A Shennong object or other class.
#' @param slot Optional: Specific element name.
#' @param ... Additional arguments.
#' @rdname sn_misc
#' @export
setGeneric("sn_misc", function(object, slot = NULL, ...) standardGeneric("sn_misc"))

#' Set Miscellaneous Data<-
#' @param object A Shennong object.
#' @param slot Name for the misc data element.
#' @param ... Additional arguments.
#' @param value The data to store.
#' @rdname sn_misc
#' @export
setGeneric("sn_misc<-", function(object, slot, ..., value) standardGeneric("sn_misc<-"))

#' Get Version Information
#' @param object A Shennong object or other class.
#' @param ... Additional arguments.
#' @rdname sn_version
#' @export
setGeneric("sn_version", function(object, ...) standardGeneric("sn_version"))
