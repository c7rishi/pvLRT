#' Extract various summary measures from a pvlrt object
#'
#' @inheritParams summary.pvlrt
#'
#' @details
#'
#' - \code{extract_lrstat_matrix} extracts the matrix of
#' the computed log-likelihood ratio test statistics. This produces
#' a result identical to applying \code{as.matrix}.
#'
#' - \code{extract_p_value_matrix} extracts the matrix of
#' computed p-values
#'
#' - \code{extract_zi_probability} extracts a vector of (estimated)
#' zero-inflation probabilities
#'
#' - \code{extract_n_matrix} extracts the original contingency table
#' used.
#'
#' - \code{extract_significant_pairs} extract the AE/drug
#' pairs determined to be significant at the provided significance level
#'
#' @export
extract_lrstat_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  as.matrix(object)
}

#' @rdname extract_lrstat_matrix
#' @export
extract_p_value_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "p_value") %>%
    set_dimnames(dimnames(object))
  out
}

#' @rdname extract_lrstat_matrix
#' @export
extract_zi_probability <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "omega") %>%
    setNames(colnames(object))
  out
}


#' @rdname extract_lrstat_matrix
#' @export
extract_n_matrix <- function(object, ...) {
  attr(object, "contin_table") %>%
    set_dimnames(dimnames(object))
}


#' @rdname extract_lrstat_matrix
#' @param significance_level numeric. Level of significance.
#' @export
extract_significant_pairs <- function(object, significance_level = 0.05, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  p_value <- NULL

  out <- summary(object) %>%
    subset(p_value < significance_level)

  out
}

#' Extracting and setting AE and Drug names from a pvlrt object
#' @inheritParams summary.pvlrt
#'
#' @note
#' Because a `pvlrt` object is simply a reclassified matrix, the AE (rows)
#' and Drug (columns) names can also be modified through \link{rownames} and
#' \link{colnames} respectively.
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

.setnames_generic <- function(object, old, new, dim) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  stopifnot(
    !missing(old),
    !missing(new),
    all(is.character(old)),
    all(is.character(new)),
    length(old) == length(new)
  )
  newnames <- oldnames <- dimnames(object)[[dim]]
  N <- length(old)
  for (jj in 1:N) {
    newnames[newnames == old[jj]] <- new[jj]
  }
  dimnames(object)[[dim]] <- newnames

  object
}

#' @rdname extract_AE_names
#' @param old character vector containing the old names
#' @param new character vector containing the new names
#' @export
set_AE_names <- function(object, old, new) {
  tmp <- tryCatch(
    .setnames_generic(object = object, old = old, new = new, dim = 1)
  )
  if (is(tmp, "error")) {
    stop(tmp$message)
  }
  tmp
}

#' @rdname extract_AE_names
#' @export
set_Drug_names <- function(object, old, new) {
  tmp <- tryCatch(
    .setnames_generic(object = object, old = old, new = new, dim = 2)
  )
  if (is(tmp, "error")) {
    stop(tmp$message)
  }
  tmp
}
