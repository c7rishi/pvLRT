.rdirichlet <- function(n, alpha) {
  # using representation X_i/sum(X_i)
  # where X_i are iid gamma
  n_alpha <- length(alpha)
  x <- rgamma(n = n_alpha * n, shape = alpha, rate = 1) %>%
    matrix(ncol = n_alpha, byrow = TRUE)
  rowsum_x <- .rowSums(x, m = n, n = n_alpha)
  x / as.vector(rowsum_x)
}


.r_contin_table_zip_1_samp <- function(row_marginals,
                                       col_marginals,
                                       obs_total = sum(row_marginals),
                                       Eij_mat = tcrossprod(row_marginals, col_marginals) / obs_total,
                                       signal_mat,
                                       omega_vec = rep(0, length(col_marginals)),
                                       no_zero_infl_idx = NULL,
                                       n_row,
                                       n_col,
                                       ...) {
  . <- NULL
  p_i0 <- c(.rdirichlet(1, row_marginals))
  p_0j <- c(.rdirichlet(1, col_marginals))
  p_ij <- (signal_mat * tcrossprod(p_i0, p_0j)) %>%
    {
      . / sum(.)
    }

  z_ij <- lapply(
    1:n_col,
    function(j) {
      rbinom(n_row, size = 1, prob = omega_vec[j])
    }
  ) %>%
    do.call(cbind, .) %>%
    .process_zero_inflation(no_zero_infl_idx)


  rmultinom(
    n = 1,
    size = obs_total,
    prob = c(p_ij)
  ) %>%
    c() %>%
    matrix(n_row, n_col) %>%
    {
      . * (1 - z_ij)
    } %>%
    magrittr::set_rownames(names(row_marginals)) %>%
    magrittr::set_colnames(names(col_marginals))
}


#' Generate random contingency tables for adverse effect
#' (across rows) and drug (across columns) combinations
#' given row and column marginal totals, embedded signals,
#' and possibly with zero inflation
#'
#' @param n number of random matrices to generate.
#' @param row_marginals,col_marginals (possibly named) vector of row and column
#' marginal totals. Must add up to the same total. If named, the names are passed
#' on to the randomly generated matrix/matrices.
#' @param signal_mat numeric matrix of dimension
#' length(row_marginals) x length(col_marginals). The (i, j)-th entry of
#' signal_mat determines the signal stregth of the i-th adverse effect and
#' j-th drug pair. The default is 1 for each pair, which means no signal for the pair.
#' @param omega_vec vector of zero inflation probabilities. Must be of the same length
#' as col_marginals.
#' @param no_zi_idx List of pairs {(i, j)} where zero inflation is not allowed. To
#' specify the entirety i-th row (or j-th column) use c(i, 0) (or c(0, j)). See examples.
#' @param ... Additional arguments. Currently unused.
#'
#' @return
#'
#' A list of length \code{n}, with each entry being a
#' \code{length(row_marginal)} by \code{length(col_marginal)} matrix.
#'
#' @examples
#'
#' set.seed(42)
#'
#' # first load the 46 statin data
#' data(statin46)
#' # zero inflation probabilities
#' omega_tru <- c(rep(0.15, ncol(statin46) - 1), 0)
#'
#' # signal matrix
#' signal_mat <- matrix(1, nrow(statin46), ncol(statin46))
#'
#' # "positive" signal at the (1, 1) entry
#' # the first column
#' signal_mat[1, 1] <- 10
#'
#' # Now simulate data with the same row/column marginals
#' # as in statin46, with embedded signals as specified in
#' # the above signal_mat
#'
#' # no zero inflation at (1, 1)
#' # (where signal is elicited) and the last row ("Other PT")
#' # and at the last column ("Other drugs") of the simulated matrix
#'
#' sim_statin <- r_contin_table_zip(
#'   n = 1,
#'   row_marginals = rowSums(statin46),
#'   col_marginals = colSums(statin46),
#'   signal_mat = signal_mat,
#'   omega_vec = omega_tru,
#'   no_zi_idx = list(
#'     c(1, 1),
#'     c(nrow(statin46), 0), # the entire last row
#'     c(0, ncol(statin46)) # the entire last column
#'   )
#' )[[1]]
#'
#' # now analyze the above simulated data
#'
#' # using a pseudo LRT with a ZIP model
#' test1 <- pvlrt(
#'   contin_table = sim_statin,
#'   nsim = 500
#'   # set to 500 for demonstration purposes only,
#'   # we recommend the default 10000 or bigger
#' )
#' extract_zi_probability(test1)
#' extract_significant_pairs(test1)
#'
#' # LRT with a poisson model
#' test2 <- lrt_poisson(
#'   contin_table = sim_statin,
#'   nsim = 500
#'   # set to 500 for demonstration purposes only,
#'   # we recommend the default 10000 or bigger
#' )
#' extract_zi_probability(test2)
#' extract_significant_pairs(test2)
#'
#' # LRT with true omega supplied
#' test3 <- pvlrt(
#'   contin_table = sim_statin,
#'   zi_prob = omega_tru,
#'   nsim = 500
#'   # set to 500 for demonstration purposes only,
#'   # we recommend the default 10000 or bigger
#' )
#' extract_zi_probability(test3)
#' extract_significant_pairs(test3)
#'
#' @export
r_contin_table_zip <- function(n = 1,
                               row_marginals,
                               col_marginals,
                               signal_mat = matrix(1, length(row_marginals), length(col_marginals)),
                               omega_vec = rep(0, length(col_marginals)),
                               no_zi_idx = NULL,
                               ...) {
  stopifnot(
    is.numeric(n),
    n > 0,
    sum(row_marginals) == sum(col_marginals),
    length(omega_vec) == length(col_marginals),
    nrow(signal_mat) == length(row_marginals),
    ncol(signal_mat) == length(col_marginals)
  )

  dots <- list(...)
  if (!is.null(dots$no_zero_inflation_idx)) {
    msg <- glue::glue(
      "argument 'no_zero_inflation_idx' is deprecated. Use
      'no_zi_idx' instead"
    )
    warning(msg)

    if (is.null(no_zi_idx)) {
      no_zi_idx <- dots$no_zero_inflation_idx
    } else {
      msg <- glue::glue(
        "Supply one of 'no_zi_idx' and 'no_zero_inflation_idx'"
      )
      stop(msg)
    }
  }

  obs_total <- sum(row_marginals)
  Eij_mat <- tcrossprod(row_marginals, col_marginals) / obs_total
  n_row <- nrow(Eij_mat)
  n_col <- ncol(Eij_mat)

  out <- lapply(
    1:n,
    function(ii) {
      .r_contin_table_zip_1_samp(
        row_marginals = row_marginals,
        col_marginals = col_marginals,
        obs_total = obs_total,
        Eij_mat = Eij_mat,
        signal_mat = signal_mat,
        omega_vec = omega_vec,
        no_zero_infl_idx = no_zi_idx,
        n_row = n_row,
        n_col = n_col,
        ...
      )
    }
  )

  out
}
