.lr_stat_pseudo_lik_1tab_zip_rrr <- function(n_ij_mat,
                                             n_i_0_all = rowSums(n_ij_mat),
                                             n_0_j_all = colSums(n_ij_mat),
                                             n_0_0 = sum(n_i_0_all),
                                             test_j_idx = 1:ncol(n_ij_mat),
                                             omega_vec,
                                             drug_class_idx = as.list(1:ncol(n_ij_mat)),
                                             grouped_omega_est = grouped_omega_est,
                                             ...) {
  I <- nrow(n_ij_mat)
  J <- ncol(n_ij_mat)
  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    set_dimnames(dimnames(n_ij_mat))

  Eij_mat_safe <- pmax(Eij_mat, 1e-20)

  hat_lambda_ij_mat <- all_log_lrt <- matrix(NA, I, J) %>%
    set_dimnames(dimnames(n_ij_mat))

  hat_lambda_ij_mat[, test_j_idx] <- pmax(
    n_ij_mat[, test_j_idx]/Eij_mat_safe[, test_j_idx],
    1
  )



  all_log_lrt[, test_j_idx] <-
    -(hat_lambda_ij_mat[, test_j_idx] - 1) * Eij_mat[, test_j_idx] +
    n_ij_mat[, test_j_idx] * log(hat_lambda_ij_mat[, test_j_idx])



  all_res <- list(
    log_lrt = all_log_lrt,
    omega = omega_vec
  )

  all_res
}

.est_zi_1tab_rrr <- function(n_ij_mat,
                             # E_ij_mat,
                             n_i_0_all = rowSums(n_ij_mat),
                             n_0_j_all = colSums(n_ij_mat),
                             n_0_0 = sum(n_i_0_all),
                             grouped_omega_est = grouped_omega_est,
                             use_gamma_smoothing = FALSE,
                             omega_constrained_lambda = TRUE,
                             test_j_idx = 1:ncol(n_ij_mat),
                             ...){
  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    set_dimnames(dimnames(n_ij_mat))

  omega_est_fn <- .estimate_zipois_mle_omega

  # fit unrestricted models separately on each column
  # to get consistent estimator of omega

  drug_class_idx_adj <- if (grouped_omega_est) {
    drug_class_idx
  }  else {
    as.list(1:ncol(n_ij_mat))
  }

  drug_class_idx_final <- drug_class_idx_adj %>%
    lapply(. %>% intersect(test_j_idx)) %>%
    .[sapply(., length) > 0]


  colnames_final <- colnames(n_ij_mat)[unlist(drug_class_idx_final)]

  omega_vec_obj <- lapply(
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
    setNames(colnames_final)

  if (ncol(n_ij_mat) > length(unlist(drug_class_idx_final))) {
    remain_cols <- unlist(drug_class_idx_final) %>%
      setdiff(1:ncol(n_ij_mat), .) %>%
      colnames(n_ij_mat)[.]
    omega_vec_obj[remain_cols] <- lapply(remain_cols, function(x) list(omega = NA, lrstat = NA))
  }

  omega_vec_obj
}


.lr_stat_pq_1tab <- function(n_ij_mat,
                             n_i_0_all = rowSums(n_ij_mat),
                             n_0_j_all = colSums(n_ij_mat),
                             n_0_0 = sum(n_i_0_all),
                             test_j_idx = 1:ncol(n_ij_mat),
                             drug_class_idx = as.list(1:ncol(n_ij_mat)),
                             omega_vec = rep(0, ncol(n_ij_mat)),
                             ...) {

  safe_zero <- function(x) {
    ifelse(x == 0, 1e-10, x)
  }


  safe_log <- function(x) {
    ifelse(x == 0, log(1e-10), log(x))
  }

  I <- nrow(n_ij_mat)
  J <- ncol(n_ij_mat)

  # N = column of ((nij)) for a specific j
  logLR_vec <- function(N, n_i, n_j, n, I){
    p_0 <- n_j/n

    p_i_vec <- N/n_i
    q_i_vec <- (n_j-N)/(n-n_i)

    logLR_vec <- ifelse(
      N == 0,
      0,
      (N*log(p_i_vec)
       +(n_j-N)*log(q_i_vec)
       -n_j*log(p_0))*(p_i_vec>q_i_vec)
    )

    logLR_vec
  }

  # perform the LRT separately for each j
  lrt_res <- lapply(
    1:J,
    function(jstar) {
      loglrt_vec <- logLR_vec(
        N = n_ij_mat[, jstar],
        n_i = n_i_0_all,
        n_j = n_0_j_all[jstar],
        n = n_0_0,
        I = I
      )

      loglrt_vec
    }
  ) %>%
    setNames(colnames(n_ij_mat)) %>%
    do.call(cbind, .)


  all_res <- list(
    log_lrt = lrt_res,
    omega = omega_vec
  )
}


