#' @export
extract_pvalue_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }

  unnecessary_attrs <- attributes(object) %>%
    names() %>%
    setdiff(c("dim", "dimnames", "class"))

  for (nm in unnecessary_attrs) {
    attr(object, nm) <- NULL
  }

  as.matrix(object)
}

#' @export
extract_lrstat_matrix <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "lrstat") %>%
    set_dimnames(dimnames(object))
  out
}

#' @export
extract_zi_probability <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  out <- attr(object, "omega") %>%
    setNames(colnames(object))
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

