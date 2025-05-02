#' Example Shennong object for demonstration
#'
#' This dataset provides a minimal example of a `Shennong` object created using
#' the `sn_create_shennong_object()` function. It contains two assays ("RNA" and "ATAC")
#' across 10 samples and 20 features each. Several data layers are included
#' (e.g., raw counts, log-transformed, and scaled data), along with sample-level
#' metadata and identity annotations.
#'
#' @format A `Shennong` object with the following components:
#' \describe{
#'   \item{assays}{Two assays: "RNA" and "ATAC"}
#'   \item{layers for RNA}{\code{counts}, \code{data} (log1p-transformed), and \code{scale.data} (scaled)}
#'   \item{active assay}{Set to "RNA"}
#'   \item{active layer}{Set to "data"}
#'   \item{sample identities}{Factor with levels "Tumor" and "Normal"}
#'   \item{sample metadata}{Column \code{sex} with values "M" and "F"}
#' }
#'
#' @usage data(so)
#' @source Simulated data generated in `data-raw/so.R`
#' @seealso \code{\link{sn_create_shennong_object}}, \code{\link{sn_layer_data}}, \code{\link{sn_add_assay}}, \code{\link{sn_add_metadata}}
#' @examples
#' data(so)
#' so
"so"
