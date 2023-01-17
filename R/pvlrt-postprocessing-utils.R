#' Extract various summary measures from a pvlrt object
#'
#' @inheritParams summary.pvlrt
#'
#' @returns
#'
#' - \code{extract_lrstat_matrix} returns the matrix of
#' the computed log-likelihood ratio test statistics for signals. This produces
#' a result identical to applying \code{as.matrix}.
#'
#' - \code{extract_p_value_matrix} returns the matrix of
#' computed p-values associated with the likelihood ratio tests.
#'
#' - \code{extract_zi_probability} returns a vector of (estimated)
#' zero-inflation probabilities.
#'
#' - \code{extract_n_matrix} returns the original contingency table (matrix)
#' used.
#'
#' - \code{extract_significant_pairs} returns a data.table listing the AE/drug
#' pairs determined to be significant at the provided significance level. This
#' is essentially a subset of the data.table obtained through summary.pvlrt()
#' that satisfies the provided significance threshold.
#'
#' - \code{extract_run_time} returns a \link{difftime} object measuring the
#'  total CPU time needed to run the original \link{pvlrt} call.
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, test_zi = TRUE, nsim = 500)
#' extract_lrstat_matrix(test1)
#' extract_p_value_matrix(test1)
#' extract_zi_probability(test1)
#' extract_n_matrix(test1)
#' extract_significant_pairs(test1)
#'
#'
#' @seealso
#' \link{pvlrt}
#'
#' @md
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



#' @rdname extract_lrstat_matrix
#' @export
extract_run_time <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  p_value <- NULL

  out <- attr(object, "run_time")

  out
}


#' Extracting and setting AE and Drug names from a pvlrt object
#' @inheritParams summary.pvlrt
#'
#' @note
#' Because a `pvlrt` object is simply a reclassified matrix, the AE (rows)
#' and Drug (columns) names can also be extracted/modified through \link{rownames} and
#' \link{colnames} respectively.
#'
#' @returns
#'
#' - `extract_AE_names` returns a character vector of the names of the
#' AEs in the input `pvlrt` object
#'
#' - `extract_Drug_names` returns a character vector of the names of the Drugs
#' in the input `pvlrt` object
#'
#' - `set_AE_names` returns a new `pvlrt` object with updated AE names as
#' specified through the arguments `old` and `new`.
#'
#' - `set_Drug_names` returns a new `pvlrt` object with updated Drug names as
#' specified through the arguments `old` and `new`.
#'
#' @examples
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, test_zi = TRUE, nsim = 500)
#' extract_AE_names(test1)
#' extract_Drug_names(test1)
#'
#' set_AE_names(test1, old = "Rhabdomyolysis", new = "Rhabdo")
#' set_Drug_names(test1, old = "Other", new = "Other-Drugs")
#'
#' ## can be chained with pipes `%>%`:
#' test2 <- test1 %>%
#'   set_AE_names(old = "Rhabdomyolysis", new = "Rhabdo") %>%
#'   set_Drug_names(old = "Other", new = "Other-Drugs")
#'
#' # see the summary for changed labels
#' summary(test2)
#'
#'
#' @seealso
#' \link{pvlrt}
#'
#' @md
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
