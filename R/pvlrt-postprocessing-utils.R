#' @export
extract_pvalue_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  class(object) <- "matrix"
  as.matrix(object)
}

#' @export
extract_lrstat_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "lrstat")
  out
}

#' @export
extract_zi_probability <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "omega")
  out
}

#' @export
extract_significant_pairs <- function(object, significance_level = 0.05, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  p.value <- NULL

  out <- summary(object) %>%
    subset(p.value < significance_level)

  out
}

#' Extracting AE and Drug names from a pvlrt object
#'
#' @export
extract_AE_names <- function(object) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  rownames(object)
}

#' @rdname extract_AE_names
#' @export
extract_Drug_names <- function(object) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  colnames(object)
}

