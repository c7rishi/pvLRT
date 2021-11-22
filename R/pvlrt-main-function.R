#' Pseudo Likelihood Ratio Test for determining significant AE-Drug pairs under
#' Poisson and zero-inflated Poisson models for pharmacovigilance
#' @inheritParams r_contin_table_zip
#' @param contin_table IxJ contingency table showing pairwise counts of adverse
#' effects for I AE and J Drugs
#' @param nsim Number of simulated null contingency table to use for computing
#' the p-value of the test
#' @param drug_class_idx a list, with the h-th entry providing the h-th group/class
#' of drugs. By default, each drug forms its own class. If more than one drug is
#' present in a class, an extended LRT is performed. See examples.
#' @param zi_prob,omega_vec Alias, determining zero inflation probabilities
#' in the model. Can be a vector, providing different zero inflation
#' probabilities for different drugs, or a scalar, providing the common zero
#' inflation probability for all drugs. If NULL (default), then is estimated from the data. See also the
#' description of the argument \code{grouped_omega_est}. If \code{omega_vec = rep(0, ncol(contin_table))},
#' then test reduces to an ordinary (non-zero inflated) Poisson test. NOTE: \code{zi_prob} and \code{omega_vec}
#' are alias.
#' @param grouped_omega_est Logical. When performing a test with grouped drug classes (extended LRT),
#' should the estimated zero-inflation parameter "omega" reflect
#' the corresponding grouping? If TRUE, then the estimated omegas are obtained by combining
#' columns from the same group, and if FALSE (default), then omegas are estimated separately for each drug (column)
#' irrespective of the groups specified through  \code{drug_class_idx}. Ignored if \code{omega_vec} is
#' supplied/non-\code{NULL} (i.e., not estimated).
#' @param test_zi,test_omega logical indicators specifying whether to perform a
#' pseudo likelihood ratio test for zero inflation. Defaults to FALSE. Ignored
#' if \code{omega_vec} is supplied (is non-NULL). Defaults to FALSE.
#' NOTE: \code{test_omega} and \code{test_zi} are aliases.
#' @param pval_ineq_strict logical. Use a strict inequality in the definition of the p-values?  Defaults to FALSE.
#' @param test_drug_idx integer vector representing the columns (drug indices) of contin_table to be tested for signal.
#' Defaults to all columns.
#' @param is_zi_structural logical. Do the inflated zeros correspond to structural
#' zeros (indicating impossible AE-Drug combination)? This determines how the
#' bootstrap null zero-inflation indicators are generated. If TRUE (default),
#' then then the null zero-inflation random indicators are randomly generated using the
#' (estimated) *conditional* probabilities of zero inflation given observed
#' counts. If FALSE, then they are generated using the *marginal* (drug-specific)
#' estimated probabilities of zero-inflation.
#' @param return_overall_loglik logical. Return overall log-likelihood for the table? This is needed
#' if \code{logLik} method is to be used.
#' @param parametrization Type of parametrization to use in the LR test. Available choices are "rrr", "lambda", "rr",
#' and "p-q". The relative reporting ratio (default) parametrization of the test is used when
#' when `parametrization %in% c("rrr", "lambda")`, and the reporting rate parametrization is used
#' otherwise. NOTE: zero inflation can be handled only for the relative reporting ratio parametrization.
#' @param null_boot_type Type of bootstrap sampling to perform for generating null resamples.
#' Available choices are "parametric" (default) and "non-parametric". NOTE: zero inflation is not
#' handled properly in a non-parametric bootstrap resampling.
#' @param ... additional arguments. Currently unused.
#'
#'
#' @returns
#'
#' The function returns an S3 object of class `pvlrt` containing  test results. A `pvlrt`
#' object is simply a re-classified matrix containing log likelihood ratio test statistics
#' for cells in the input contingency table, with various other test and input data information (including
#' p-values, estimated zero inflation parameters, overall log-likelihood etc.) embedded
#' as attributes. The following S3 methods and functions are available for an `pvlrt` object:
#'
#'
#'
#' Various postprocessing methods for `pvlrt` objects are available. This includes:
#'
#' * \link{bubbleplot_pvlrt}
#'
#' * \link{extract_AE_names}
#'
#' * \link{extract_Drug_names}
#'
#' * \link{extract_lrstat_matrix}
#'
#' * \link{extract_n_matrix}
#'
#' * \link{extract_p_value_matrix}
#'
#' * \link{extract_significant_pairs}
#'
#' * \link{extract_zi_probability}
#'
#' * \link{heatmap_pvlrt}
#'
#' * \link{lrt_poisson}
#'
#' * \link{lrt_vanilla_poisson}
#'
#' * \link{lrt_zi_poisson}
#'
#' * \link{r_contin_table_zip}
#'
#' * \link{set_AE_names}
#'
#' * \link{set_Drug_names}
#'
#' * \link{print.pvlrt}
#'
#' * \link{plot.pvlrt}
#'
#' * \link{summary.pvlrt}
#'
#' * \link{logLik.pvlrt}
#'
#' * \link{as.matrix.pvlrt}
#'
#'
#'
#'
#' @examples
#'
#' data("statin46")
#'
#' # 500 bootstrap iterations (nsim) in each example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' # no grouping -- each drug forms its own class,
#' # default "rrr" (lambda) parametrization, possible zi,
#' # null distribution evaluated using parametric bootstrap
#' test1 <- pvlrt(statin46, nsim = 500)
#' test1
#' ## extract the observed LRT statistic
#' extract_lrstat_matrix(test1)
#' ## extract the estimated omegas
#' extract_zi_probability(test1)
#'
#' # grouped drugs --
#' # group 1: drug 1, drug 2
#' # group 2: drug 3, drug 4
#' # drug 5, 6, 7 in their own groups
#' drug_groups <- list(c(1, 2), c(3, 4), 5, 6, 7)
#' test2 <- pvlrt(statin46, drug_class_idx = drug_groups, nsim = 500)
#' test2
#'
#'
#' # specify no zero inflation at the entirety of the last row and the last column
#' no_zi_idx <- list(c(nrow(statin46), 0), c(0, ncol(statin46)))
#' test3 <- pvlrt(statin46, no_zi_idx = no_zi_idx, nsim = 500)
#' test3
#'
#' \donttest{
#' # use non-parametric bootstrap to evaluate the null distribution
#' # takes longer, due to computational costs with non-parametric
#' # contigency table generation
#' test4 <- pvlrt(statin46, null_boot_type = "non-parametric", nsim = 500)
#' test4
#' }
#'
#' # test zi probabilities (omegas)
#' test5 <- pvlrt(statin46, test_omega = TRUE, nsim = 500)
#' test5
#'
#'
#' @md
#' @export
pvlrt <- function(contin_table,
                  nsim = 1e4,
                  omega_vec = NULL,
                  zi_prob = NULL,
                  no_zi_idx = NULL,
                  drug_class_idx = as.list(1:ncol(contin_table)),
                  test_drug_idx = 1:ncol(contin_table),
                  grouped_omega_est = FALSE,
                  test_zi = NULL,
                  test_omega = NULL,
                  pval_ineq_strict = FALSE,
                  parametrization = "rrr",
                  null_boot_type = "parametric",
                  is_zi_structural = TRUE,
                  return_overall_loglik = TRUE,
                  ...) {
  . <- NULL

  contin_table <- tryCatch(
    data.matrix(contin_table),
    error = function(e) e
  )

  if (is(contin_table, "error")) {
    msg <- glue::glue(
      "'contin_table' cannot be converted into a matrix:
      {contin_table$message}"
    )
  }


  nm_type <- c("row", "col")
  charc_type <- c("AE", "Drug")
  contin_dim <- dim(contin_table)

  for (kk in 1:2) {
    nm_type_fn <- match.fun(paste0(nm_type[kk], "names"))
    `nm_type_fn<-` <- match.fun(paste0(nm_type[kk], "names<-"))
    this_nms <- nm_type_fn(contin_table)
    if (any(is.null(this_nms))) {
      msg <- glue::glue(
        "'contin_table' has no {nm_type[kk]} names. \\
        Created {nm_type[kk]} names: {charc_type[kk]}_1, {charc_type[kk]}_2 etc. \\
        You can replace them posthoc using `set_{charc_type[kk]}_names()` \\
        or `{nm_type[kk]}names()` on the generated pvlrt object."
      )
      message(msg)
      nm_type_fn(contin_table) <- paste(
        charc_type[kk],
        1:contin_dim[kk],
        sep = "_"
      )
    }
  }
  stopifnot(
    is.list(drug_class_idx),
    is.logical(grouped_omega_est),
    length(grouped_omega_est) == 1,
    is.logical(pval_ineq_strict),
    length(pval_ineq_strict) == 1,
    parametrization %in% c("rrr", "lambda", "rr", "p-q"),
    length(parametrization) == 1,
    null_boot_type %in% c("parametric", "non-parametric"),
    length(null_boot_type) == 1,
    is.logical(is_zi_structural),
    length(is_zi_structural) == 1
  )

  if (!is.null(no_zi_idx)) {
    stopifnot(
      is.list(no_zi_idx)
    )
  }

  dots <- list(...)
  all_inputs <- c(as.list(environment()), dots)

  I <- nrow(contin_table)
  J <- ncol(contin_table)

  n_0_0 <- sum(contin_table)
  n_i_0_all <- rowSums(contin_table)
  n_0_j_all <- colSums(contin_table)

  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all) / n_0_0) %>%
    set_dimnames(dimnames(contin_table))


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

  # will be multiplied with zi indicator
  # during random generation
  # zi: 1 means zi 0 means no zi
  zi_idx_adj <- matrix(1, I, J) %>%
    set_dimnames(dimnames(contin_table)) %>%
    .process_zero_inflation(
      no_zero_infl_idx = no_zi_idx
    )


  # parallel_bootstrap <- dot$parallel_bootstrap %>%
  #   ifelse(is.null(.), TRUE, .)
  parallel_bootstrap <- FALSE

  `%>or>=%` <- if (pval_ineq_strict) `>` else `>=`


  for (nm in c("omega_est_vec", "omega_vec")) {
    if (!is.null(dots[[nm]])) {
      msg <- glue::glue(
        "The argument '{nm}' is deprecated. \\
        Use 'omega_vec' instead."
      )
      warning(msg)
      if (is.null(omega_vec)) {
        omega_vec <- dots[[nm]]
      }
    }
  }

  for (nm in c(
    "skip_null_omega_estimation",
    "use_gamma_smooth_omega"
  )) {
    if (!is.null(dots[[nm]])) {
      msg <- glue::glue(
        "The argument '{nm}' is deprecated and unused."
      )
    }
  }


  if (is.null(zi_prob) && !is.null(omega_vec)) {
    zi_prob <- omega_vec
  } else if (!is.null(zi_prob) && is.null(omega_vec)) {
    omega_vec <- zi_prob
  } else if (!is.null(zi_prob) && !is.null(omega_vec)) {
    if (!identical(zi_prob, omega_vec)) {
      stop("'zi_prob' and 'omega_vec' both supplied and they differ")
    }
  }

  if (!is.null(omega_vec)) {
    omega_vec <- zi_prob <- as.numeric(omega_vec)

    if (length(omega_vec) == 1) {
      omega_vec <- zi_prob <- rep(omega_vec, J)
    }
    if (is.null(names(omega_vec))) {
      names(omega_vec) <- colnames(contin_table)
    }
    if (any((omega_vec < 0 | omega_vec >= 1))) {
      msg <- glue::glue(
        "'zi_prob' (or 'omega_vec') must be betwen [0, 1)."
      )
      stop(msg)
    }
  }


  if (is.null(test_zi) && is.null(test_omega)) {
    test_zi <- test_omega <- FALSE
  } else if (is.null(test_zi) && !is.null(test_omega)) {
    test_zi <- test_omega
  } else if (!is.null(test_zi) && is.null(test_omega)) {
    test_omega <- test_zi
  } else if (!is.null(test_zi) && !is.null(test_omega)) {
    if (!identical(test_zi, test_omega)) {
      stop("'test_zi' and 'test_omega' both supplied and they differ")
    }
  }

  stopifnot(
    is.logical(test_zi),
    length(test_zi) == 1,
    is.logical(test_omega),
    length(test_omega) == 1
  )


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
    {
      . == 0
    }

  if (!len_check_1) {
    stop("'drug_class_idx' contains more columns than 'contin_table'")
  }


  len_check_2 <- c(1:ncol(contin_table)) %>%
    setdiff(sort(unlist(drug_class_idx))) %>%
    length() %>%
    {
      . == 0
    }
  if (!len_check_2) {
    stop("'drug_class_idx' does not contain all columns of 'contin_table'")
  }


  do_omega_estimation <- is.null(omega_vec) &
    parametrization %in% c("rrr", "lambda")

  if (!is.null(omega_vec)) {
    skip_omega_step <- TRUE
  }


  if (parametrization %in% c("rrr", "lambda")) {
    lr_stat_func <- .lr_stat_pseudo_lik_1tab_zip_rrr
  } else {
    unused_args <- c(
      "omega_vec", "omega_est_vec", "zi_prob",
      "test_omega", "test_zi"
    )

    for (nm in unused_args) {
      chk <- is.null(all_inputs[[nm]])
      if (!chk) {
        msg <- glue::glue(
          "{nm} is ignored if \\
          parameterization %in% c('rr', 'p-q')"
        )
        warning(msg)
      }
    }
    lr_stat_func <- .lr_stat_pq_1tab
    omega_vec <- rep(0, J) %>% setNames(colnames(contin_table))
  }

  omega_lrstat_vec <- omega_pval_vec <- NULL

  if (do_omega_estimation) {
    omega_vec_obj <- .est_zi_1tab_rrr(
      n_ij_mat = contin_table,
      grouped_omega_est = grouped_omega_est,
      # use_gamma_smoothing = use_gamma_smooth_omega,
      omega_constrained_lambda = omega_constrained_lambda
    )

    omega_vec <- sapply(omega_vec_obj, "[[", "omega")
    omega_lrstat_vec <- sapply(omega_vec_obj, "[[", "lrstat")
    this_grouped_omega_est <- grouped_omega_est



    # likelihood ratio test for zero inflation
    if (test_zi) {
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
        # lapply(
        #   1:J,
        #   function(jstar) {
        #     nn <- length(poisson_mean_hat[, jstar])
        #     rpois(
        #       n = nn,
        #       lambda = poisson_mean_hat[, jstar]
        #     )
        #   }
        # ) %>%
        #   do.call(cbind, .) %>%
        #   set_dimnames(dimnames(contin_table))

        rpois(I * J, lambda = c(poisson_mean_hat)) %>%
          matrix(
            nrow = I, ncol = J,
            dimnames = dimnames(contin_table)
          )
      }

      omega_lr_stat_function <- function(n_ij_table) {
        tmp <- .est_zi_1tab_rrr(
          n_ij_mat = n_ij_table,
          grouped_omega_est = this_grouped_omega_est,
          # use_gamma_smoothing = use_gamma_smooth_omega,
          omega_constrained_lambda = omega_constrained_lambda,
          test_j_idx = test_drug_idx
        )
        sapply(tmp, "[[", "lrstat")
      }

      pval <- rep(0, J) %>%
        setNames(colnames(contin_table))


      # if (!parallel_bootstrap) {
      pb <- progress::progress_bar$new(
        format = " simulating [:bar] :percent eta: :eta",
        total = nsim + 1,
        clear = FALSE,
        width = 60
      )

      # pb <- progressr::progressor(steps = nsim + 1)

      for (ii in 1:(nsim + 1)) {
        if (ii <= nsim) {
          lr_stat_omega_null <- gen_rand_table_for_omega() %>%
            omega_lr_stat_function()
        } else {
          lr_stat_omega_null <- omega_lrstat_vec
        }

        pval_increment <- 1 * (lr_stat_omega_null %>or>=% omega_lrstat_vec)

        pval <- pval + pval_increment

        pb$tick()
        # pb()
      }

      pval <- pval / (nsim + 1)

      #   } else {
      #
      #     pval <- furrr::future_map(
      #       1:nsim,
      #       function(ii) {
      #         out <- gen_rand_table_for_omega() %>%
      #           omega_lr_stat_function() %>%
      #           {. %>or>=% omega_lrstat_vec}
      #         pb()
      #       },
      #       # future.seed = TRUE
      #       .options = furrr::furrr_options(seed = TRUE),
      #       .progress = TRUE
      #     ) %>%
      #       c(list(omega_lrstat_vec %>or>=% omega_lrstat_vec)) %>%
      #       Reduce("+", .) %>%
      #       {./(nsim + 1)}
      #   }
      omega_pval_vec <- pval
    }
  }

  msg <- paste(
    "Calculating observed",
    ifelse(do_omega_estimation, "pseudo LRT", "LRT"),
    "statistic for global null against signals...\n"
  )

  cat(msg)


  lr_stat_obs_obj <- lr_stat_func(
    contin_table,
    n_i_0_all = n_i_0_all,
    n_0_j_all = n_0_j_all,
    n_0_0 = n_0_0,
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

  omega_mat_orig <- rep(1, I) %>%
    tcrossprod(omega_vec) %>%
    set_dimnames(dimnames(contin_table))

  omega_mat <- omega_mat_orig %>%
    {
      if (is_zi_structural) {
        # posterior probabilities of zi
        ifelse(
          contin_table == 0,
          .safe_divide(., (. + (1 - .) * exp(-Eij_mat))),
          0
        )
      } else {
        .
      }
    }


  if (null_boot_type == "parametric") {
    gen_rand_table <- function(ii) {
      nij <- rpois(I * J, c(Eij_mat))
      zij <- (runif(I * J) <= c(omega_mat)) * c(zi_idx_adj)

      (nij * (1 - zij)) %>%
        matrix(I, J) %>%
        set_dimnames(dimnames(contin_table))
    }
  } else {
    all_2d_tabs <- r2dtable(nsim, r = n_i_0_all, c = n_0_j_all) %>%
      lapply(set_dimnames, dimnames(contin_table))

    gen_rand_table <- function(ii) {
      all_2d_tabs[[ii]]
    }
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
    set_dimnames(dimnames(contin_table))


  # if (!parallel_bootstrap) {
  pb <- progress::progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = nsim + 1,
    clear = FALSE,
    width = 60
  )

  # pb <- progressr::progressor(steps = nsim + 1)

  pval[, test_drug_idx] <- 0

  for (ii in 1:(nsim + 1)) {
    if (ii <= nsim) {
      lr_stat_null <- `class<-`(NA, "error")
      ntry <- 0

      while (is(lr_stat_null, "error") & ntry < 100) {
        lr_stat_null <- tryCatch(
          gen_rand_table(ii) %>%
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

    pval_increment <- 1L * (
      mlr_stat[, test_drug_idx] %>or>=%
        lr_stat_obs[, test_drug_idx]
    )

    pval[, test_drug_idx] <- pval[, test_drug_idx] + pval_increment


    # pb()
    pb$tick()
  }

  pval[, test_drug_idx] <- pval[, test_drug_idx] / (nsim + 1)

  # } else {
  #
  #   pval[, test_drug_idx] <- furrr::future_map(
  #     1:nsim,
  #     function(ii) {
  #       out <- gen_rand_table() %>%
  #         lr_stat_func(
  #           n_0_0 = n_0_0,
  #           omega_vec = omega_vec,
  #           test_j_idx = test_drug_idx,
  #           drug_class_idx = drug_class_idx,
  #           grouped_omega_est = grouped_omega_est
  #         ) %>%
  #         .[["log_lrt"]] %>%
  #         calc_mlr() %>%
  #         {1L * (.[, test_drug_idx] %>or>=% lr_stat_obs[, test_drug_idx])}
  #
  #       pb()
  #
  #       out
  #     },
  #     # future.seed = TRUE
  #     .options = furrr::furrr_options(seed = TRUE),
  #     .progress = TRUE
  #   ) %>%
  #     c(
  #       list(
  #         lr_stat_obs %>%
  #           calc_mlr() %>%
  #           {1L * (.[, test_drug_idx] %>or>=% lr_stat_obs[, test_drug_idx])}
  #       )
  #     ) %>%
  #     Reduce("+", .) %>%
  #     {./(nsim + 1)}
  # }

  lr_stat_p_value <- pval


  # log likelihoods under 4 possible setups
  # loglik_full_zip <- loglik_full_poisson <-
  #   loglik_null_zip <- loglik_null_poisson <- NULL
  loglik_df_comb <- NULL


  if (return_overall_loglik & !(parametrization %in% c("lambda", "rrr"))) {
    msg <- glue::glue(
      "overall log likelihood is only computed for parameterization = 'rrr' or
      'lambda'. Setting return_overall_loglik = FALSE"
    )
    warning(msg)
  }

  if (return_overall_loglik & parametrization %in% c("lambda", "rrr")) {
    loglik_df_comb <- expand.grid(
      mod = c("full", "null"),
      dist = c("poisson", "zip")
    ) %>%
      data.table::setDT()

    lambda_times_Eij <- list(
      full = pmax(contin_table, Eij_mat),
      null = Eij_mat
    )

    # this is the estimated omega mat
    # NOT the posterior probabilities
    omega_mat_est <- omega_mat_orig

    log_lik_fun <- list(
      poisson = function(theta) {
        sum(
          dpois(
            x = c(contin_table),
            lambda = c(theta),
            log = TRUE
          )
        )
      },
      zip = function(theta) {
        sum(
          dzipois(
            x = c(contin_table),
            lambda = c(theta),
            pi = c(omega_mat_est),
            log = TRUE
          )
        )
      }
    )

    mod <- dist <- NULL

    loglik_df_comb[
      ,
      `:=`(
        type = paste(mod, dist, sep = "-"),
        logLik = mapply(
          function(mod, dist) {
            this_theta <- lambda_times_Eij[[mod]]
            this_log_lik_fun <- log_lik_fun[[dist]]
            this_log_lik_fun(this_theta)
          },
          mod,
          dist,
          SIMPLIFY = TRUE
        ),
        df = data.table::fcase(
          mod == "full" & dist == "poisson", as.double(I * J),
          mod == "null" & dist == "poisson", as.double(0.0),
          mod == "full" & dist == "zip", ifelse(
            do_omega_estimation,
            I * J + J,
            I * J
          ) %>% as.double(.),
          mod == "null" & dist == "zip", ifelse(
            do_omega_estimation, J, 0.0
          ) %>% as.double(.)
        ),
        N = n_0_0
      ),
    ]

    type <- logLik <- df <- N <- NULL
    loglik_df_comb <- loglik_df_comb[
      ,
      list(type, logLik, df, N)
    ]
  }

  attrs <- list(
    p_value = lr_stat_p_value,
    # lrstat = lr_stat_obs,
    omega = omega_vec,
    test_omega = test_omega,
    omega_lrstat = omega_lrstat_vec,
    omega_p_value = omega_pval_vec,
    omega_qvalue = p.adjust(omega_pval_vec),
    do_omega_estimation = do_omega_estimation,
    parametrization = parametrization,
    test_drug_idx = test_drug_idx,
    contin_table = contin_table,
    no_zi_idx = no_zi_idx,
    null_boot_type = null_boot_type,
    is_zi_structural = is_zi_structural,
    return_overall_loglik = return_overall_loglik
  ) %>%
    lapply(unname) %>%
    # need the names for loglik_df_comb
    c(list(loglik_df = loglik_df_comb))

  # output <- lr_stat_p_value
  output <- lr_stat_obs %>%
    set_dimnames(dimnames(contin_table))

  attributes(output) <- attributes(output) %>% c(attrs)

  class(output) <- c("pvlrt", class(lr_stat_p_value))

  output
}
