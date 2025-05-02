#' Demo Shennong Object
#'
#' A compact example `Shennong` object for demonstration and testing. Includes two assays,
#' multiple data layers, dimensionality reductions, and sample metadata.
#'
#' @format A `Shennong` object with:
#' \describe{
#'   \item{2 assays}{\code{RNA}, \code{ATAC}}
#'   \item{3 layers in RNA}{\code{counts}, \code{data} (log1p), \code{scale.data} (z-scored)}
#'   \item{2 reductions}{\code{pca} (via \code{prcomp}), \code{mds} (via \code{cmdscale})}
#'   \item{20 samples with metadata}{
#'     \code{condition} (\code{Normal} or \code{Tumor}),\cr
#'     \code{gender} (\code{Male} or \code{Female}),\cr
#'     \code{age} (integer from 40–70)
#'   }
#'   \item{200 genes}{Shared across both assays}
#' }
#'
#' @usage data(so)
#'
#' @examples
#' data(so)
#' so
#' sn_assays(so)
#' sn_layers(so, assay = "RNA")
#' sn_reductions(so)
#' sn_metadata(so)
#'
"so"
