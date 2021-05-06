.lr_stat_1tab <- function(n_ij_mat,
                          n_i_0_all = rowSums(n_ij_mat),
                          n_0_j_all = colSums(n_ij_mat),
                          n_0_0 = sum(n_i_0_all),
                          ...) {

  `%>%` <- magrittr::`%>%`

  safe_zero <- function(x) {
    ifelse(x == 0, 1e-10, x)
  }


  safe_log <- function(x) {
    ifelse(x == 0, log(1e-10), log(x))
  }

  I <- nrow(n_ij_mat)
  J <- ncol(n_ij_mat)

  logLR <- function(N, n_i, n_j, n, I){
    p_0 <- n_j/n
    logLR <- c()
    for (i in 1:I){
      if (N[i]==0){
        logLR[i] <- 0
      } else{
        p_i <- N[i]/n_i[i]
        q_i <- (n_j-N[i])/(n-n_i[i])
        logLR[i] <- (N[i]*log(p_i)+(n_j-N[i])*log(q_i)-n_j*log(p_0))*(p_i>q_i)
      }
    }
    logLR
  }


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

      # loglrt <- logLR(
      #   N = n_ij_mat[, jstar],
      #   n_i = n_i_0_all,
      #   n_j = n_0_j_all[jstar],
      #   n = n_0_0,
      #   I = I
      # )

      loglrt_vec <- logLR_vec(
        N = n_ij_mat[, jstar],
        n_i = n_i_0_all,
        n_j = n_0_j_all[jstar],
        n = n_0_0,
        I = I
      )

      # #
      #       n_0_jstar <- n_0_j_all[jstar]
      #       # n_0_minus.jstar <- n_0_0 - n_0_jstar
      #
      #       n_i_jstar_all <- n_ij_mat[, jstar]
      # n_i_minus.jstar_all <- n_i_0_all - n_i_jstar_all


      # n_minus.i_0_all <- n_0_0 - n_i_0_all
      # n_minus.i_jstar_all <- n_0_jstar - n_i_jstar_all


      # lr_all <- n_i_jstar_all * safe_log(n_i_jstar_all/safe_zero(n_i_0_all)) +
      #   n_minus.i_jstar_all * safe_log(n_minus.i_jstar_all/safe_zero(n_minus.i_0_all)) -
      #   n_0_jstar * safe_log(n_0_jstar/n_0_0)

      # lr_all <- n_i_jstar_all * safe_log(n_i_jstar_all/safe_zero(n_i_0_all)) +
      #     (n_0_jstar - n_i_jstar_all) * safe_log(
      #       (n_0_jstar - n_i_jstar_all) /
      #         safe_zero(n_0_0 - n_i_0_all)
      #     ) -
      #     n_0_jstar * safe_log(n_0_jstar/n_0_0)

      loglrt_vec

    }
  ) %>%
    setNames(colnames(n_ij_mat)) %>%
    do.call(cbind, .)


  lrt_res
}

#' Likelihood Ratio Test for determining significant AE-Drug pairs
#' @param contin_table IxJ contingency table showing pairwise counts of adverse effects
#' for I AE and J Drugs
#' @param nsim Number of simulated contin_table to use for computing the p-value of the test
#'
#' @export
lrt_vanilla_poisson <- function(contin_table,
                                nsim = 1e4, ...)
{
  I <- nrow(contin_table)
  J <- ncol(contin_table)

  n_0_0 <- sum(contin_table)
  n_i_0_all <- rowSums(contin_table)
  n_0_j_all <- colSums(contin_table)



  cat("Calculating observed LR stat...\n")
  lr_stat_obs <- .lr_stat_1tab(
    contin_table,
    # n_i_0_all = n_i_0,
    # n_0_j_all = n_0_j.all,
    n_0_0 = n_0_0
  )


  mlr_obs <- apply(lr_stat_obs, 2, max)


  # generate random contingency tables
  # with fixed row and column sums from multinomial
  cat("simulating random null contingency tables..\n")

  # rand_contin_tab_list <- r2dtable(
  #   n = nsim,
  #   r = n_i_0_all,
  #   c = n_0_j_all
  # )

  rand_contin_tab_list <- pbapply::pblapply(
    1:nsim,
    function(this_idx) {
      lapply(
        n_0_j_all,
        function(this_n_0_j) {
          c(rmultinom(n = 1, size = this_n_0_j, prob = n_i_0_all/n_0_0)) %>%
            setNames(names(n_i_0_all))
        }
      ) %>%
        setNames(names(n_0_j_all)) %>%
        do.call(cbind, .)
    }
  )

  cat("simulating null distribution of the lr stat..\n")
  lr_stat_null <- pbapply::pblapply(
    rand_contin_tab_list,
    .lr_stat_1tab,
    n_0_0 = n_0_0,
  ) %>%
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

  lr_stat_pvalue
}
