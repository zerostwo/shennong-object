#' Filter samples in a Shennong object using dplyr syntax
#'
#' Allows tidyverse-style filtering of sample metadata (`colData`) in a Shennong object.
#' This behaves like `dplyr::filter()`, and subsets the object to retain only the samples
#' matching the specified conditions.
#'
#' @param .data A `Shennong` object.
#' @param ... Filtering conditions passed to `dplyr::filter()`.
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
filter.Shennong <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  if (!inherits(.data, "Shennong")) {
    stop("The object must be a Shennong object.")
  }

  # 1. Extract colData and convert to tibble for filtering
  meta <- colData(.data)
  meta_tbl <- as_tibble(meta, .name_repair = "minimal") %>%
    mutate(.sample_id = rownames(meta))

  # 2. Apply dplyr filtering
  filtered_tbl <- filter(meta_tbl, ..., .by = .by, .preserve = .preserve)

  if (nrow(filtered_tbl) == 0) {
    abort("No samples matched the filter criteria!")
  }

  # 3. Get filtered sample IDs
  kept_samples <- filtered_tbl$.sample_id

  # 4. Subset Shennong object (samples = columns)
  .data <- .data[, kept_samples]

  # 5. Return the filtered object
  return(.data)
}
