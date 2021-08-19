#' Pseudo Likelihood Ratio Test for determining significant AE-Drug pairs under
#' zero-inflated Poisson model
#' @inheritParams lrt_vanilla_poisson
#' @param omega_vec vector (for all drugs) of estimates of the zero-inflation.
#' If NULL, then is estimated from the data. See also the
#' description of the argument \code{grouped_omega_est}. If \code{omega_vec = rep(0, ncol(contin_table))},
#' then test reduces to an ordinary (non-zero inflated) Poisson test.
#' @param grouped_omega_est Logical. When performing a test with grouped drug classes (extended LRT),
#' should the estimated zero-inflation parameter "omega" reflect
#' the corresponding grouping? If TRUE, then the estimated omegas are obtained by combining
#' columns from the same group, and if FALSE (default), then omegas are estimated separately for each drug (column)
#' irrespective of the groups specified through  \code{drug_class_idx}. Ignored if \code{omega_vec} is
#' supplied/non-\code{NULL} (i.e., not estimated).
#' @param test_omega logical. perform a bootstrap pseudo likelihood ratio test for omega? Ignored if
#' \code{omega_vec} is supplied (is non-NULL). Defaults to FALSE.
#' @param pval_ineq_strict logical. Use a strict inequality in the definition of the p-values?  Defaults to FALSE.
#' @param use_gamma_smooth_omega logical. Use a gamma prior (smoothing) on the signals (lambdas) while estimating
#' omega from data. Defaults to FALSE. NOTE: if TRUE, then the gamma prior produces a marginal negative binomial
#' distribution for the counts, which is then optimized for estimating omega.
#' @param test_drug_idx integer vector representing the columns (drug indices) of contin_table to be tested for signal.
#' Defaults to all columns.
#'
#'
#' @note \code{lrt_poisson} is a wrapper function on \code{lrt_zi_poisson} with \code{omega_vec = rep(0, ncol(contin_table))}
#'
#' @examples
#'
#' ## uncomment the following to show progress
#' # progressr::handlers("progress")
#'
#'
#' data("lovastatin")
#' # no grouping -- each drug forms its own class
#' test1 <- lrt_zi_poisson(lovastatin)
#' ## extract the observed LRT statistic
#' attr(test1, "lrstat")
#' ## extract the estimated omegas
#' attr(test1, "omega")
#'
#' # grouped drugs --
#' # group1 : drug 1, drug 2
#' # group 2: drug 3
#' drug_groups <- list(c(1, 2), 3)
#' test2 <- lrt_zi_poisson(lovastatin, drug_class_idx = drug_groups)
#'
#' \dontrun{
#' # test omegas
#' test3 <- lrt_zi_poisson(lovastatin, test_omega = TRUE)
#' ## extract the estimate, statistic, and p-values for omegas
#' attr(test3, "omega")
#' attr(test3, "omega_lrstat")
#' attr(test3, "omega_pvalue")
#'
#'
#' # use gamma prior assumption on signals while
#' # estimating omegas
#' test3 <- lrt_zi_poisson(
#'   lovastatin,
#'   use_gamma_smooth_omega = TRUE
#' )
#' attr(test3, "omega")
#' }
#'
#' @export
lrt_zi_poisson <- function(contin_table,
                           nsim = 1e4,
                           omega_vec = NULL,
                           drug_class_idx = as.list(1:ncol(contin_table)),
                           test_drug_idx = 1:ncol(contin_table),
                           grouped_omega_est = FALSE,
                           test_omega = FALSE,
                           pval_ineq_strict = FALSE,
                           use_gamma_smooth_omega = FALSE,
                           ...)
{


  stopifnot(
    is.list(drug_class_idx),
    is.matrix(contin_table),
    # length(dim(contin_table)) == 2,
    # is.null(omega_vec), (!is.null(omega_vec) && is.numeric(omega_vec)),
    # is.list(drug_class_idx),
    is.logical(grouped_omega_est),
    length(grouped_omega_est) == 1,
    is.logical(test_omega),
    length(test_omega) == 1,
    is.logical(pval_ineq_strict),
    length(pval_ineq_strict) == 1,
    is.logical(use_gamma_smooth_omega),
    length(use_gamma_smooth_omega) == 1
  )



  I <- nrow(contin_table)
  J <- ncol(contin_table)

  n_0_0 <- sum(contin_table)
  n_i_0_all <- rowSums(contin_table)
  n_0_j_all <- colSums(contin_table)

  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(contin_table))



  dots <- list(...)

  # parallel_bootstrap <- dot$parallel_bootstrap %>%
  #   ifelse(is.null(.), TRUE, .)
  parallel_bootstrap <- FALSE

  `%>or>=%` <- if (pval_ineq_strict) `>` else `>=`


  if (!is.null(dots$omega_est_vec)) {
    warning("The argument 'omega_est_vec' is deprecated. Use 'omega_vec' instead")
    if (is.null(omega_vec)) {
      omega_vec <- dots$omega_est_vec
    }
  }

  if (!is.null(dots$skip_null_omega_estimation)) {
    warning("The argument 'skip_null_omega_estimation' is deprecated and unused.")
  }

  omega_constrained_lambda <- dots$omega_constrained_lambda

  if (is.null(omega_constrained_lambda)) {
    omega_constrained_lambda <- TRUE
  }

  if (!is.logical(omega_constrained_lambda)) {
    omega_constrained_lambda <- TRUE
  }


  len_check_1 <- sort(unlist(drug_class_idx)) %>%
    setdiff(1:ncol(contin_table)) %>%
    length() %>%
    {. == 0}

  if (!len_check_1) {
    stop("'drug_class_idx' contains more columns than 'contin_table'")
  }


  len_check_2 <- c(1:ncol(contin_table)) %>%
    setdiff(sort(unlist(drug_class_idx))) %>%
    length() %>%
    {. == 0}
  if (!len_check_2) {
    stop("'drug_class_idx' does not contain all columns of 'contin_table'")
  }


  do_omega_estimation <- is.null(omega_vec)

  # skip_omega_step <- !is.null(omega_vec) | skip_null_omega_estimation

  # skip_omega_step <- skip_null_omega_estimation

  if (!is.null(omega_vec)) {
    skip_omega_step <- TRUE
  }


  omega_lrstat_vec <- omega_pval_vec <-  NULL


  if (do_omega_estimation) {

    omega_vec_obj <- .est_zi_1tab_rrr(
      n_ij_mat = contin_table,
      grouped_omega_est = grouped_omega_est,
      use_gamma_smoothing = use_gamma_smooth_omega,
      omega_constrained_lambda = omega_constrained_lambda
    )

    omega_vec <- sapply(omega_vec_obj, "[[", "omega")
    omega_lrstat_vec <- sapply(omega_vec_obj, "[[", "lrstat")
    this_grouped_omega_est <- grouped_omega_est


    if (test_omega) {
      msg <- paste(
        "\nBootstrapping null distribution of the",
        "pseudo LRT statistic for omega",
        "and calculating p-values..\n\n"
      )

      cat(msg)


      # E_ij_adj <- pmax(Eij_mat, 1e-20)
      # lambda_ij_hat <- pmax(contin_table/E_ij_adj, 1)
      poisson_mean_hat <- if (omega_constrained_lambda) {
        # corresponds to \hat lambda_ij = max(1, n_ij/E_ij)
        pmax(contin_table, Eij_mat)
      } else {
        contin_table
      }


      gen_rand_table_for_omega <- function() {
        lapply(
          1:J,
          function(jstar) {
            nn <- length(poisson_mean_hat[, jstar])
            rpois(
              n = nn,
              lambda = poisson_mean_hat[, jstar]
            )
          }
        ) %>%
          do.call(cbind, .) %>%
          `dimnames<-`(dimnames(contin_table))
      }

      omega_lr_stat_function <- function(n_ij_table) {
        tmp <- .est_zi_1tab_rrr(
          n_ij_mat = n_ij_table,
          grouped_omega_est = this_grouped_omega_est,
          use_gamma_smoothing = use_gamma_smooth_omega,
          omega_constrained_lambda = omega_constrained_lambda,
          test_j_idx = test_drug_idx
        )
        sapply(tmp, "[[", "lrstat")
      }



      pval <- rep(0, J) %>%
        setNames(colnames(contin_table))


      if (!parallel_bootstrap) {
        pb <- progress::progress_bar$new(
          format = " simulating [:bar] :percent eta: :eta",
          total = nsim + 1,
          clear = FALSE,
          width = 60
        )

        # pb <- progressr::progressor(steps = nsim + 1)

        for (ii in 1:(nsim+1)) {
          if (ii <= nsim) {
            lr_stat_omega_null <- gen_rand_table_for_omega() %>%
              omega_lr_stat_function()
          } else {
            lr_stat_omega_null <- omega_lrstat_vec
          }

          pval <- pval + 1 * (lr_stat_omega_null %>or>=% omega_lrstat_vec)
          pb$tick()
          # pb()
        }
        pval <- pval/(nsim + 1)

      } else {

        pval <- furrr::future_map(
          1:nsim,
          function(ii) {
            out <- gen_rand_table_for_omega() %>%
              omega_lr_stat_function() %>%
              {. %>or>=% omega_lrstat_vec}
            pb()
          },
          # future.seed = TRUE
          .options = furrr::furrr_options(seed = TRUE),
          .progress = TRUE
        ) %>%
          c(list(omega_lrstat_vec %>or>=% omega_lrstat_vec)) %>%
          Reduce("+", .) %>%
          {./(nsim + 1)}
      }
      omega_pval_vec <- pval
    }



  }

  msg <- paste(
    "Calculating observed",
    ifelse(do_omega_estimation, "pseudo LRT", "LRT"),
    "statistic for global null against signals...\n"
  )

  cat(msg)

  lr_stat_func <- .lr_stat_pseudo_lik_1tab_zip_rrr

  lr_stat_obs_obj <- lr_stat_func(
    contin_table,
    n_0_0 = n_0_0,
    use_stan = use_stan,
    test_j_idx = test_drug_idx,
    omega_vec = omega_vec,
    drug_class_idx = drug_class_idx,
    grouped_omega_est = grouped_omega_est
  )
  lr_stat_obs <- lr_stat_obs_obj$log_lrt



  # generate random contingency tables
  # with fixed row and column sums from multinomial
  msg <- paste(
    "\nBootstrapping null distribution of the",
    ifelse(do_omega_estimation, "pseudo LRT", "LRT"),
    "statistics for signals & calculating p-values..\n\n"
  )
  cat(msg)

  drug_class_idx_adj <- drug_class_idx %>%
    lapply(. %>% intersect(test_drug_idx)) %>%
    .[sapply(., length) > 0]

  gen_rand_table <- function() {
    lapply(
      1:J,
      function(jstar) {
        lambda <- c(Eij_mat[, jstar])
        nn <- length(lambda)
        rzipois(
          n = nn,
          lambda = lambda,
          pi = omega_vec[jstar]
        )
      }
    ) %>%
      do.call(cbind, .)
  }

  calc_mlr <- function(xmat) {
    # n_row <- nrow(xmat)
    col_maxs <- rep(NA, J)
    col_maxs[test_drug_idx] <- xmat[, test_drug_idx, drop = FALSE] %>%
      apply(2, max)

    for (ii in 1:length(drug_class_idx_adj)) {
      drug_class <- drug_class_idx_adj[[ii]]
      if (length(drug_class) > 0) {
        col_maxs[drug_class] <- max(col_maxs[drug_class])
      }
    }

    out <- matrix(NA, I, J)
    for (jj in test_drug_idx) {
      out[, jj] <- col_maxs[jj]
    }
    out
  }

  pval <- matrix(NA, I, J) %>%
    `dimnames<-`(dimnames(contin_table))


  if (!parallel_bootstrap) {
    pb <- progress::progress_bar$new(
      format = " simulating [:bar] :percent eta: :eta",
      total = nsim + 1,
      clear = FALSE,
      width = 60
    )

    # pb <- progressr::progressor(steps = nsim + 1)

    pval[, test_drug_idx] <- 0

    for (ii in 1:(nsim+1)) {

      if (ii <= nsim) {
        lr_stat_null <- `class<-`(NA, "error")
        ntry <- 0

        while (is(lr_stat_null, "error") & ntry < 100) {
          lr_stat_null <- tryCatch(
            gen_rand_table() %>%
              lr_stat_func(
                n_0_0 = n_0_0,
                omega_vec = omega_vec,
                test_j_idx = test_drug_idx,
                drug_class_idx = drug_class_idx,
                grouped_omega_est = grouped_omega_est
              ) %>%
              .[["log_lrt"]],
            error = function(e) e
          )
          ntry <- ntry + 1
        }

        if (ntry > 1) {

          if (is(lr_stat_null, "error")) stop(lr_stat_null)

          msg <- paste0(
            "errors encountered while computing", ii,
            "th null bootstrapped test statistic. Resampled another null table.",
            " Several of such errors may indicate problems with the original data"
          )
          warning(msg)
        }

      } else {
        lr_stat_null <- lr_stat_obs
      }

      mlr_stat <- calc_mlr(lr_stat_null)

      pval[, test_drug_idx] <- pval[, test_drug_idx] +
        1L * (mlr_stat[, test_drug_idx] %>or>=% lr_stat_obs[, test_drug_idx])

      # pb()
      pb$tick()
    }

    pval[, test_drug_idx] <- pval[, test_drug_idx]/(nsim + 1)

  } else {

    pval[, test_drug_idx] <- furrr::future_map(
      1:nsim,
      function(ii) {
        out <- gen_rand_table() %>%
          lr_stat_func(
            n_0_0 = n_0_0,
            omega_vec = omega_vec,
            test_j_idx = test_drug_idx,
            drug_class_idx = drug_class_idx,
            grouped_omega_est = grouped_omega_est
          ) %>%
          .[["log_lrt"]] %>%
          calc_mlr() %>%
          {1L * (.[, test_drug_idx] %>or>=% lr_stat_obs[, test_drug_idx])}

        pb()

        out
      },
      # future.seed = TRUE
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    ) %>%
      c(
        list(
          lr_stat_obs %>%
            calc_mlr() %>%
            {1L * (.[, test_drug_idx] %>or>=% lr_stat_obs[, test_drug_idx])}
        )
      ) %>%
      Reduce("+", .) %>%
      {./(nsim + 1)}
  }

  lr_stat_pvalue <- pval


  attr(lr_stat_pvalue, "lrstat") <- lr_stat_obs
  attr(lr_stat_pvalue, "omega") <- omega_vec
  attr(lr_stat_pvalue, "omega_lrstat") <- omega_lrstat_vec
  attr(lr_stat_pvalue, "omega_pvalue") <- omega_pval_vec
  attr(lr_stat_pvalue, "do_omega_estimation") <- do_omega_estimation
  attr(lr_stat_pvalue, "parametrization") <- "rrr"
  attr(lr_stat_pvalue, "test_drug_idx") <- test_drug_idx

  class(lr_stat_pvalue) <- c("pvlrt", "matrix")

  lr_stat_pvalue
}



#' @rdname lrt_zi_poisson
#'
#' @export
lrt_poisson <- function(contin_table,
                        nsim = 1e4, ...)
{

  dots <- list(...)
  if (!is.null(dots$omega_vec) | !is.null(dots$omega_est_vec)) {
    msg <- paste(
      "omega_vec' is automatically set to 0 in lrt_poisson().",
      "For non-zero omegas use lrt_zi_poisson(), and supply as 'omega_vec'"
    )
    dots$omega_vec <- dots$omega_est_vec <- NULL
  }

  inargs <- list(
    contin_table = contin_table,
    nsim = nsim,
    omega_vec = rep(0, ncol(contin_table))
  ) %>%
    c(dots)

  do.call(lrt_zi_poisson, inargs)
}

