# .lr_stat_brute_force_1tab_zip <- function(n_ij_mat,
#                                           n_i_0_all = rowSums(n_ij_mat),
#                                           n_0_j_all = colSums(n_ij_mat),
#                                           n_0_0 = sum(n_i_0_all),
#                                           test_j_idx = 1:ncol(n_ij_mat),
#                                           ...) {
#
#   I <- nrow(n_ij_mat)
#   J <- ncol(n_ij_mat)
#   Eij_mat <- tcrossprod(n_i_0_all, n_0_j_all)/n_0_0
#
#   # will be updated separately for each pair
#   lambda_ij_test_indic <- matrix(0, I, J)
#
#
#   # jstar <- 1
#
#
#   # browser()
#   all_log_lrt <- pbapply::pbmapply(
#     function(jstar, do_this_test) {
#       if (do_this_test) {
#         res <- lapply(
#           1:I,
#           function(this_i) {
#             tmp_lambda_ij_test_indic <- lambda_ij_test_indic
#             tmp_lambda_ij_test_indic[this_i, jstar] <- 1
#
#             stan_dat <- list(
#               I = I,
#               J = 1,
#               n_ij = n_ij_mat[, jstar, drop=F],
#               E_ij = Eij_mat[, jstar, drop=F],
#               lambda_1_indic_ij = tmp_lambda_ij_test_indic[, jstar, drop=F]
#             )
#             stan_opt_null <- rstan::optimizing(
#               STAN_MODEL_NULL__,
#               data = stan_dat
#             )
#
#             stan_opt_alt <- rstan::optimizing(
#               STAN_MODEL_ALT__,
#               data = stan_dat
#             )
#
#             # log_lrt <- max(stan_opt_alt$value - stan_opt_null$value, 0)
#             # attr(log_lrt, "omega_null") <- stan_opt_null$par["omega"]
#             # attr(log_lrt, "omega_null") <- stan_opt_null$par["omega"]
#             list(
#               log_lrt = max(stan_opt_alt$value - stan_opt_null$value, 0),
#               omega_null = stan_opt_null$par["omega"] %>% unname(),
#               omega_alt = stan_opt_alt$par["omega"] %>% unname()
#             )
#           }
#         )
#       } else {
#         res <- lapply(
#           1:I,
#           function(this_i) {
#             list(
#               log_lrt = 0,
#               omega_null = NA,
#               omega_alt = NA
#             )
#           }
#         )
#       }
#       setNames(res, rownames(n_ij_mat))
#     },
#     1:J,
#     1:J %in% test_j_idx,
#     SIMPLIFY = FALSE
#   )  %>%
#     setNames(colnames(n_ij_mat))
#
#
#   all_res <-
#     c("log_lrt", "omega_null", "omega_alt") %>%
#     setNames(., .) %>%
#     lapply(
#       function(this_elem) {
#         all_log_lrt %>%
#           lapply(
#             function(x) {
#               x %>%
#                 lapply("[[", this_elem) %>%
#                 unlist()
#             }
#           ) %>%
#           do.call(cbind, .)
#       }
#     )
#
#
#   # out <- all_res$log_lrt
#   # attr(out, "omega_null") <- all_res$omega_null
#   # attr(out, "omega_alt") <- all_res$omega_alt
#   all_res
# }


.lr_stat_pseudo_lik_1tab_zip <- function(n_ij_mat,
                                         n_i_0_all = rowSums(n_ij_mat),
                                         n_0_j_all = colSums(n_ij_mat),
                                         n_0_0 = sum(n_i_0_all),
                                         test_j_idx = 1:ncol(n_ij_mat),
                                         omega_est_vec = NULL,
                                         ...) {
  I <- nrow(n_ij_mat)
  J <- ncol(n_ij_mat)
  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(n_ij_mat))

  # will be updated separately for each pair
  lambda_ij_test_indic <- matrix(0, I, J)

  if (is.null(omega_est_vec)) {
    # fit unrestricted models separately on each column
    # to get consistent estimator of omega

    fit_unrestricted <- lapply(
      1:J,
      function(jstar) {
        .estimate_zigammapois_mle(
          n_ij_mat[, jstar, drop = FALSE],
          Eij_mat[, jstar, drop = FALSE]
        )
      }
    ) %>%
      setNames(colnames(n_ij_mat))

    # map_dbl(kk, "llik")
    #

    # map_dbl(fit_unrestricted, "llik")

    omega_est_vec <- sapply(fit_unrestricted, "[[", "omega")
  }

  names(omega_est_vec) <- colnames(omega_est_vec)

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


#' Pseudo Likelihood Ratio Test for determining significant AE-Drug pairs under
#' zero-inflated Poisson model
#' @inheritParams lrt_vanilla_poisson
#' @param omega_est_vec vector (for all drugs) of estimates of the zero-inflation.
#' If NULL, then estimated from the data under a Gamma process assumption.
#' @export
lrt_zi_poisson <- function(contin_table,
                           nsim = 1e4,
                           omega_est_vec = NULL,
                           ...)
{
  I <- nrow(contin_table)
  J <- ncol(contin_table)

  n_0_0 <- sum(contin_table)
  n_i_0_all <- rowSums(contin_table)
  n_0_j_all <- colSums(contin_table)

  Eij_mat <- (tcrossprod(n_i_0_all, n_0_j_all)/n_0_0) %>%
    `dimnames<-`(dimnames(contin_table))


  # browser()


  cat("Calculating observed LR stat...\n")

  lr_stat_func <- .lr_stat_pseudo_lik_1tab_zip

  lr_stat_obs_obj <- lr_stat_func(
    contin_table,
    # n_i_0_all = n_i_0,
    # n_0_j_all = n_0_j.all,
    n_0_0 = n_0_0,
    use_stan = use_stan,
    omega_est_vec = omega_est_vec
  )
  lr_stat_obs <- lr_stat_obs_obj$log_lrt
  omega_est_vec <- lr_stat_obs_obj$omega %>%
    setNames(colnames(contin_table))


  mlr_obs <- apply(lr_stat_obs, 2, max)


  # generate random contingency tables
  # with fixed row and column sums from multinomial
  cat("simulating random null contingency tables..\n")

  # rand_contin_tab_list <- r2dtable(
  #   n = nsim,
  #   r = n_i_0_all,
  #   c = n_0_j_all
  # )

  # rand_contin_tab_list <- pbapply::pblapply(
  #   1:nsim,
  #   function(this_idx) {
  #     lapply(
  #       n_0_j_all,
  #       function(this_n_0_j) {
  #         c(rmultinom(n = 1, size = this_n_0_j, prob = n_i_0_all/n_0_0)) %>%
  #           setNames(names(n_i_0_all))
  #       }
  #     ) %>%
  #       setNames(names(n_0_j_all)) %>%
  #       do.call(cbind, .)
  #   }
  # )

  rand_contin_tab_list <- pbapply::pblapply(
    1:nsim,
    function(this_idx) {
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
        do.call(cbind, .) %>%
        `dimnames<-`(dimnames(contin_table))
    }
  )



  cat("simulating null distribution of the lr stat..\n")
  lr_stat_null <- pbapply::pblapply(
    rand_contin_tab_list,
    lr_stat_func,
    n_0_0 = n_0_0,
    omega_est_vec = omega_est_vec
  ) %>%
    lapply("[[", "log_lrt") %>%
    c(list(lr_stat_obs))

  mlr_stat_null <- lapply(
    lr_stat_null,
    function(xmat) {
      n_row <- nrow(xmat)
      apply(xmat, 2, function(x) rep(max(x), n_row))
    }
  )


  # browser()

  cat("calculating simulated p for lr test...\n")
  lr_stat_pvalue <- mlr_stat_null %>%
    lapply(function(x) c(x > lr_stat_obs)) %>%
    do.call(rbind, .) %>%
    apply(2, mean, na.rm = TRUE) %>%
    matrix(
      nrow = I,
      ncol = J,
      dimnames = dimnames(contin_table)
    )


  attr(lr_stat_pvalue, "lrstat") <- lr_stat_obs
  attr(lr_stat_pvalue, "omega") <- omega_est_vec


  lr_stat_pvalue
}


# lr_test_zip_poisson(obs_mat)

# lovastatin_res1 <- .lr_stat_brute_force_1tab_zip(obs_mat)
# lovastatin_res2 <- .lr_stat_pseudo_lik_1tab_zip(obs_mat)
# lovastatin_res %>%
#   lapply(round, 3)

# cbind(
#   brute_force = lovastatin_res1$log_lrt[, 3],
#   pseudo_llik = lovastatin_res2$log_lrt[, 3]
# ) %>%
#   as_tibble(
#     rownames = "AE"
#   ) %>%
#   mutate_if(is.numeric, round, 2) %>%
#   data.table::as.data.table()
