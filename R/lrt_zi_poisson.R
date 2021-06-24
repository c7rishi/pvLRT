.lr_stat_pseudo_lik_1tab_zip <- function(n_ij_mat,
                                         n_i_0_all = rowSums(n_ij_mat),
                                         n_0_j_all = colSums(n_ij_mat),
                                         n_0_0 = sum(n_i_0_all),
                                         test_j_idx = 1:ncol(n_ij_mat),
                                         omega_est_vec,
                                         drug_class_idx = as.list(1:ncol(n_ij_mat)),
                                         grouped_omega_est = grouped_omega_est,
                                         ...) {
  I <- nrow(n_ij_mat)
  J <- ncol(n_ij_mat)
  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(n_ij_mat))

  # will be updated separately for each pair
  lambda_ij_test_indic <- matrix(0, I, J)

  # browser()




  all_log_lrt <- mapply(
    function(jstar, do_this_test) {
      if (do_this_test) {
        theta_est_vec_list <- list(
          # corresponds to \hat lambda_ij = 1
          null =  Eij_mat[, jstar],

          # corresponds to \hat lambda_ij = max(1, n_ij/E_ij)
          alt = pmax(n_ij_mat[, jstar], Eij_mat[, jstar])
        )



        log_den_vec_list <- lapply(
          theta_est_vec_list,
          function(theta_vec) {
            dzipois(
              x = n_ij_mat[, jstar],
              lambda = theta_vec,
              pi = omega_est_vec[jstar],
              log = TRUE
            )
          }
        )
        res <- pmax(log_den_vec_list$alt - log_den_vec_list$null, 0)
      } else {
        res <- rep(0, I)
      }
      setNames(res, rownames(n_ij_mat))
    },
    1:J,
    1:J %in% test_j_idx,
    SIMPLIFY = FALSE
  )  %>%
    setNames(colnames(n_ij_mat)) %>%
    do.call(cbind, .)


  all_res <- list(
    log_lrt = all_log_lrt,
    omega = omega_est_vec
  )

  all_res
}

.est_omega_1tab <- function(n_ij_mat,
                            # E_ij_mat,
                            n_i_0_all = rowSums(n_ij_mat),
                            n_0_j_all = colSums(n_ij_mat),
                            n_0_0 = sum(n_i_0_all),
                            grouped_omega_est = grouped_omega_est,
                            use_gamma_smoothing = FALSE,
                            omega_constrained_lambda = TRUE,
                            ...){
  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(n_ij_mat))

  omega_est_fn <- if (use_gamma_smoothing) {
    .estimate_zigammapois_mle_omega
  } else {
    .estimate_zipois_mle_omega
  }

  # fit unrestricted models separately on each column
  # to get consistent estimator of omega

  drug_class_idx_final <- if (grouped_omega_est) {
    drug_class_idx
  }  else {
    as.list(1:ncol(n_ij_mat))
  }

  omega_est_vec_obj <- lapply(
    drug_class_idx_final,
    function(jstar_list) {
      est <- omega_est_fn(
        n_ij_mat[, jstar_list, drop = FALSE],
        Eij_mat[, jstar_list, drop = FALSE],
        omega_constrained_lambda = omega_constrained_lambda
      )
      nn <- length(jstar_list)

      list(
        omega = rep(est$omega, nn),
        lrstat = rep(est$lrstat, nn)
      )
    }
  ) %>%
    # unlist() %>%
    setNames(colnames(n_ij_mat)[unlist(drug_class_idx_final)]) %>%
    .[colnames(n_ij_mat)]

  omega_est_vec_obj

}


#' Pseudo Likelihood Ratio Test for determining significant AE-Drug pairs under
#' zero-inflated Poisson model
#' @inheritParams lrt_vanilla_poisson
#' @param omega_est_vec vector (for all drugs) of estimates of the zero-inflation.
#' If NULL, then is estimated from the data. See also the
#' description of the argument \code{grouped_omega_est}. If \code{omega_est_vec = rep(0, ncol(contin_table))},
#' then test reduces to an ordinary (non-zero inflated) Poisson test.
#' @param grouped_omega_est Logical. When performing a test with grouped drug classes (extended LRT),
#' should the estimated zero-inflation parameter "omega" reflect
#' the corresponding grouping? If TRUE, then the estimated omegas are obtained by combining
#' columns from the same group, and if FALSE (default), then omegas are estimated separately for each drug (column)
#' irrespective of the groups specified through  \code{drug_class_idx}. Ignored if \code{omega_est_vec} is
#' supplied/non-\code{NULL} (i.e., not estimated).
#' @param test_omega logical. perform a bootstrap pseudo likelihood ratio test for omega? Ignored if
#' \code{omega_est_vec} is supplied (is non-NULL). Defaults to FALSE.
#' @param pval_ineq_strict logical. Use a strict inequality in the definition of the p-values?  Defaults to FALSE.
#' @param skip_null_omega_estimation logical. Should estimation of the (nuisance parameter) omega be skipped while simulating
#' bootstrapped null distribution of the likelihood ratio statistic? Defaults to FALSE. If TRUE,
#' then the same global (possibly estimated) omega is used in each null booststrap draw.
#' @param use_gamma_smooth_omega logical. Use a gamma prior (smoothing) on the signals (lambdas) while estimating
#' omega from data. Defaults to FALSE. NOTE: if TRUE, then the gamma prior produces a marginal negative binomial
#' distribution for the counts, which is then optimized for estimating omega.
#'
#' @examples
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
#' }
#'
#' @export
lrt_zi_poisson <- function(contin_table,
                           nsim = 1e4,
                           omega_est_vec = NULL,
                           drug_class_idx = as.list(1:ncol(contin_table)),
                           test_drug_idx = 1:ncol(contin_table),
                           grouped_omega_est = FALSE,
                           test_omega = FALSE,
                           pval_ineq_strict = FALSE,
                           skip_null_omega_estimation = FALSE,
                           use_gamma_smooth_omega = FALSE,
                           ...)
{
  I <- nrow(contin_table)
  J <- ncol(contin_table)

  n_0_0 <- sum(contin_table)
  n_i_0_all <- rowSums(contin_table)
  n_0_j_all <- colSums(contin_table)

  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(contin_table))


  stopifnot(
    is.list(drug_class_idx)
  )


  dots <- list(...)
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


  # browser()

  # skip_omega_step <- !is.null(omega_est_vec) | skip_null_omega_estimation

  skip_omega_step <- skip_null_omega_estimation

  if (!is.null(omega_est_vec)) {
    skip_omega_step <- skip_null_omega_estimation <- TRUE
  }


  omega_lrstat_vec <- omega_pval_vec <-  NULL


  if (is.null(omega_est_vec)) {
    omega_est_vec_obj <- .est_omega_1tab(
      n_ij_mat = contin_table,
      grouped_omega_est = grouped_omega_est,
      use_gamma_smoothing = use_gamma_smooth_omega,
      omega_constrained_lambda = omega_constrained_lambda
    )

    omega_est_vec <- sapply(omega_est_vec_obj, "[[", "omega")
    omega_lrstat_vec <- sapply(omega_est_vec_obj, "[[", "lrstat")
    this_grouped_omega_est <- grouped_omega_est


    if (test_omega) {

      cat("\nBootstrapping null distribution of the pseudo LRT statistic for omega, and calculating p-values..\n")


      # E_ij_adj <- pmax(Eij_mat, 1e-20)
      # lambda_ij_hat <- pmax(contin_table/E_ij_adj, 1)
      poisson_mean_hat <- if (omega_constrained_lambda) {
        # corresponds to \hat lambda_ij = max(1, n_ij/E_ij)
        pmax(contin_table, E_ij)
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
        tmp <- .est_omega_1tab(
          n_ij_mat = n_ij_table,
          grouped_omega_est = this_grouped_omega_est,
          use_gamma_smoothing = use_gamma_smooth_omega,
          omega_constrained_lambda = omega_constrained_lambda
        )
        sapply(tmp, "[[", "lrstat")
      }



      pval <- rep(0, J) %>%
        setNames(colnames(contin_table))


      pb <- progress::progress_bar$new(
        format = " simulating [:bar] :percent eta: :eta",
        total = nsim + 1,
        clear = FALSE,
        width = 60
      )


      for (ii in 1:(nsim+1)) {
        rand_mat <- gen_rand_table_for_omega()

        if (ii <= nsim) {
          lr_stat_omega_null <- omega_lr_stat_function(rand_mat)
        } else {
          lr_stat_omega_null <- omega_lrstat_vec
        }

        pval <- if (pval_ineq_strict) {
          pval + 1 * (lr_stat_omega_null > omega_lrstat_vec)
        } else {
          pval + 1 * (lr_stat_omega_null >= omega_lrstat_vec)
        }

        pb$tick()
      }

      pval <- pval/(nsim + 1)

      omega_pval_vec <- pval
    }



  }


  cat("Calculating observed LR stat for global null for signal...\n")


  lr_stat_func <- .lr_stat_pseudo_lik_1tab_zip

  lr_stat_obs_obj <- lr_stat_func(
    contin_table,
    n_0_0 = n_0_0,
    use_stan = use_stan,
    test_j_idx = test_drug_idx,
    omega_est_vec = omega_est_vec,
    drug_class_idx = drug_class_idx,
    grouped_omega_est = grouped_omega_est
  )
  lr_stat_obs <- lr_stat_obs_obj$log_lrt
  # omega_est_vec <- lr_stat_obs_obj$omega %>%
  #   setNames(colnames(contin_table))




  # generate random contingency tables
  # with fixed row and column sums from multinomial
  cat("\nBootstrapping null distribution of the pseudo LRT statistics for signals & calculating p-values..")
  txt <- ifelse(
    skip_null_omega_estimation,
    "\nglobal omega supplied in the null distribution of signals...\n",
    "\nomega is estimated in the null distribution of signals...\n"
  )
  cat(txt)

  gen_rand_table <- function() {
    lapply(
      1:J,
      function(jstar) {
        lambda <- c(Eij_mat[, jstar])
        nn <- length(lambda)
        rzipois(
          n = nn,
          lambda = lambda,
          pi = omega_est_vec[jstar]
        )
      }
    ) %>%
      do.call(cbind, .) #%>%
    # `dimnames<-`(dimnames(contin_table))
  }

  calc_mlr <- function(xmat) {
    # n_row <- nrow(xmat)

    col_maxs <- apply(xmat, 2, max)
    for (ii in 1:length(drug_class_idx)) {
      drug_class <- drug_class_idx[[ii]]
      col_maxs[drug_class] <- max(col_maxs[drug_class])
    }
    out <- tcrossprod(rep(1, I), col_maxs)
    out
  }

  pb <- progress::progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = nsim + 1,
    clear = FALSE,
    width = 60
  )

  pval <- matrix(0, I, J) %>%
    `dimnames<-`(dimnames(contin_table))

  # browser()

  for (ii in 1:(nsim+1)) {
    rand_mat <- gen_rand_table()

    if (ii <= nsim) {

      if (skip_null_omega_estimation) {
        this_omega_est <- omega_est_vec
      } else {
        this_omega_est <- .est_omega_1tab(
          n_ij_mat = rand_mat %>%
            `dimnames<-`(dimnames(contin_table)),
          grouped_omega_est = grouped_omega_est,
          use_gamma_smoothing = use_gamma_smooth_omega
        ) %>%
          sapply("[[", "omega")
      }


      lr_stat_null <- lr_stat_func(
        rand_mat,
        n_0_0 = n_0_0,
        omega_est_vec = this_omega_est,
        test_j_idx = test_drug_idx,
        drug_class_idx = drug_class_idx,
        grouped_omega_est = grouped_omega_est
      )[["log_lrt"]]
    } else {
      lr_stat_null <- lr_stat_obs
    }

    mlr_stat <- calc_mlr(lr_stat_null)
    pval <- if (pval_ineq_strict) {
      pval + 1 * (mlr_stat > lr_stat_obs)
    } else {
      pval + 1 * (mlr_stat >= lr_stat_obs)
    }
    pb$tick()
  }

  pval <- pval/(nsim + 1)

  lr_stat_pvalue <- pval


  attr(lr_stat_pvalue, "lrstat") <- lr_stat_obs
  attr(lr_stat_pvalue, "omega") <- omega_est_vec
  attr(lr_stat_pvalue, "omega_lrstat") <- omega_lrstat_vec
  attr(lr_stat_pvalue, "omega_pvalue") <- omega_pval_vec


  lr_stat_pvalue
}

