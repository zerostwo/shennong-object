#' ggplot method for Shennong object
#'
#' Allows tidy ggplot2 visualization using sample-level metadata or expression values.
#'
#' @param data A Shennong object.
#' @param mapping Aesthetic mapping (e.g., \code{aes(x = ..., y = ..., color = ...)}).
#' @param ... Additional arguments passed to \code{ggplot()}.
#' @param environment The environment in which to evaluate the \code{mapping}. Usually left as default.
#' @param vars Optional character vector of variables to fetch using \code{sn_fetch_data()}.
#'
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes
#' @method ggplot Shennong
#' @export
ggplot.Shennong <- function(data, mapping = aes(), ..., environment = parent.frame(), vars = NULL) {
  # Try to automatically detect variables used in aes mapping
  all_vars <- unique(unlist(lapply(mapping, rlang::as_label)))
  vars <- vars %||% all_vars

  df <- sn_fetch_data(data, vars = vars)

  ggplot2::ggplot(data = df, mapping = mapping, ..., environment = environment)
}
