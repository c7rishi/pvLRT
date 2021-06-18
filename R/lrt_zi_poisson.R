.lr_stat_pseudo_lik_1tab_zip <- function(n_ij_mat,
                                         n_i_0_all = rowSums(n_ij_mat),
                                         n_0_j_all = colSums(n_ij_mat),
                                         n_0_0 = sum(n_i_0_all),
                                         test_j_idx = 1:ncol(n_ij_mat),
                                         omega_est_vec = NULL,
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




  if (is.null(omega_est_vec)) {
    # fit unrestricted models separately on each column
    # to get consistent estimator of omega

    drug_class_idx_final <- if (grouped_omega_est) {
      drug_class_idx
    }  else {
      as.list(1:ncol(n_ij_mat))
    }

    omega_est_vec <- lapply(
      drug_class_idx_final,
      function(jstar_list) {
        est <- .estimate_zigammapois_mle(
          n_ij_mat[, jstar_list, drop = FALSE],
          Eij_mat[, jstar_list, drop = FALSE]
        )
        rep(est$omega, length(jstar_list))
      }
    ) %>%
      unlist() %>%
      setNames(colnames(n_ij_mat)[unlist(drug_class_idx_final)]) %>%
      .[colnames(n_ij_mat)]

  } else if (any(is.null(names(omega_est_vec)))) {
    names(omega_est_vec) <- colnames(n_ij_mat)
  }


  all_log_lrt <- mapply(
    function(jstar, do_this_test) {
      if (do_this_test) {
        theta_est_vec_list <- list(
          # corresponds to \hat lambda_ij = 1
          null =  Eij_mat[, jstar],

          # corresponds to \hat lambda_ij = max(1, n_ij/E_ij)
          alt = pmax(n_ij_mat[, jstar], Eij_mat[, jstar])
        )


        # if (any(is.na(unlist(theta_est_vec_list)) | is.null(unlist(theta_est_vec_list)))) browser()
        # if (any(is.na(unlist(omega_est_vec)) | is.null(unlist(omega_est_vec)))) browser()

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


#' Pseudo Likelihood Ratio Test for determining significant AE-Drug pairs under
#' zero-inflated Poisson model
#' @inheritParams lrt_vanilla_poisson
#' @param omega_est_vec vector (for all drugs) of estimates of the zero-inflation.
#' If NULL, then is estimated from the data under a Gamma process assumption. See also the
#' description of the argument \code{grouped_omega_est}. If \code{omega_est_vec = rep(0, ncol(contin_table))},
#' then test reduces to an ordinary (non-zero inflated) Poisson test.
#' @param grouped_omega_est Logical. When performing a test with grouped drug classes (extended LRT),
#' should the estimated zero-inflation parameter "omega" reflect
#' the corresponding grouping? If TRUE, then the estimated omegas are obtained by combining
#' columns from the same group, and if FALSE (default), then omegas are estimated separately for each drug (column)
#' irrespective of the groups specified through  \code{drug_class_idx}. Ignored if \code{omega_est_vec} is
#' supplied/non-\code{NULL} (i.e., not estimated).
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
#'
#' @export
lrt_zi_poisson <- function(contin_table,
                           nsim = 1e4,
                           omega_est_vec = NULL,
                           drug_class_idx = as.list(1:ncol(contin_table)),
                           test_drug_idx = 1:ncol(contin_table),
                           grouped_omega_est = FALSE,
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


  cat("Calculating observed LR stat...\n")


  lr_stat_func <- .lr_stat_pseudo_lik_1tab_zip

  lr_stat_obs_obj <- lr_stat_func(
    contin_table,
    # n_i_0_all = n_i_0,
    # n_0_j_all = n_0_j.all,
    n_0_0 = n_0_0,
    use_stan = use_stan,
    test_j_idx = test_drug_idx,
    omega_est_vec = omega_est_vec,
    drug_class_idx = drug_class_idx,
    grouped_omega_est = grouped_omega_est
  )
  lr_stat_obs <- lr_stat_obs_obj$log_lrt
  omega_est_vec <- lr_stat_obs_obj$omega %>%
    setNames(colnames(contin_table))


  # mlr_obs <- .compute_mlr(
  #   lr_mat = lr_stat_obs,
  #   dt_drug_by_class = dt_drug_class_idx,
  #   I = I,
  #   dim_names = dimnames(contin_table)
  # )


  # generate random contingency tables
  # with fixed row and column sums from multinomial
  cat("\nSimulating null distribution of the pseudo LRT statistics & calculating p-values..\n")

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
    # browser()
    out <- tcrossprod(rep(1, I), col_maxs) #%>%
    # `dimnames<-`(dimnames(contin_table))

    # out <- apply(xmat, 2, function(x) rep(max(x), n_row))
    # .compute_mlr(
    #   xmat,
    #   dt_drug_by_class = dt_drug_class_idx,
    #   I = I,
    #   dim_names = dimnames(contin_table)
    # )
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

  for (ii in 1:(nsim+1)) {
    rand_mat <- gen_rand_table()

    if (ii <= nsim) {
      lr_stat_null <- lr_stat_func(
        rand_mat,
        n_0_0 = n_0_0,
        omega_est_vec = omega_est_vec,
        test_j_idx = test_drug_idx,
        drug_class_idx = drug_class_idx,
        grouped_omega_est = grouped_omega_est
      )[["log_lrt"]]
    } else {
      lr_stat_null <- lr_stat_obs
    }

    mlr_stat <- calc_mlr(lr_stat_null)
    pval <- pval + 1 * (mlr_stat >= lr_stat_obs)

    pb$tick()
  }

  pval <- pval/(nsim + 1)

  lr_stat_pvalue <- pval

  attr(lr_stat_pvalue, "lrstat") <- lr_stat_obs
  attr(lr_stat_pvalue, "omega") <- omega_est_vec


  lr_stat_pvalue
}

