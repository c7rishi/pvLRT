# --------------------------------------------------------------
# Zero inflated Poisson-Gamma -- density & random generation
# --------------------------------------------------------------
# taken from R package "countreg"
dzipois <- function(x, lambda, pi, log = FALSE)
{
  if (any(pi < 0) | any(pi > 1))
    warning("'pi' must be in [0, 1]")
  rval <- log(1 - pi) + dpois(x, lambda = lambda, log = TRUE)
  if (any(x0 <- (x == 0L)))
    rval[x0] <- log(exp(rval) + pi)[x0]
  if (log)
    rval
  else exp(rval)
}

rzipois <- function (n, lambda, pi)
{
  if (any(pi < 0) | any(pi > 1))
    warning("'pi' must be in [0, 1]")
  rval <- rpois(n, lambda = lambda)
  rval[runif(n) < pi] <- 0
  rval
}
# --------------------------------------------------------------


# --------------------------------------------------------------
# Poisson-Gamma mixture -- density & random generation
# --------------------------------------------------------------
dpoisgamma <- function(x, nu = 1, sigma = 1, log = FALSE) {
  gamma <- sigma/(1 + sigma)
  out <-  (nu + x) * log(gamma) + lgamma(nu + x) -
    nu * log(sigma) - lgamma(nu) - lfactorial(x)

  if (!log) out <- exp(out)

  out
}

# # use R's negative bionomial formulation
# dpoisgamma <- function(x, nu = 1, sigma = 1, log = FALSE) {
#   prob <- 1/(1 + sigma)
#   dnbinom(x, size = nu, prob = prob, log = log)
# }

rpoisgamma <- function(n = 1, nu = 1, sigma = 1) {
  theta <- rgamma(n = n, shape = nu, scale = sigma)
  out <- rpois(n, lambda = theta)
  out
}

# use R's neg binomial
# rpoisgamma <- function(n = 1, nu = 1, sigma = 1) {
#   prob <- 1/(1 + sigma)
#   rnbinom(n, size = nu, prob = prob)
# }

# --------------------------------------------------------------
# Zero inflated Poisson-Gamma -- density & random generation
# --------------------------------------------------------------
dzipoisgamma <- function(x, nu = 1, sigma = 1, pi = 0, log = FALSE) {
  idx_0 <- which(x == 0)

  den <- log(1 - pi) + dpoisgamma(
    x,
    nu = nu,
    sigma = sigma,
    log = TRUE
  )

  den[idx_0] <- log(pi + exp(den[idx_0]))

  if (!log) den <- exp(den)

  den
}

rzipoisgamma <- function(n = 1, nu = 1, sigma = 1, pi = 0, ...) {
  out <- rep(0, n)
  non0_idx <- (runif(n) <= 1 - pi) # from poisgamma w.p. 1-pi
  n_non0_idx <- sum(non0_idx)
  if (n_non0_idx > 0) {
    out[non0_idx] <- rpoisgamma(
      n = n_non0_idx,
      nu = nu,
      sigma = sigma
    )
  }
  out
}
# --------------------------------------------------------------


# # ---------------------------------------------------------------
# # MM Estimation of ZIPoissonGamma parameters
# # ---------------------------------------------------------------
# mme_dzipoisgamma <- function(r = NULL,
#                              Nr = NULL,
#                              x = NULL,
#                              return_r_Nr = FALSE,
#                              ...) {
#
#   if (is.null(r) | is.null(Nr)) {
#     tmp <- c(table(x))
#     r <- as.numeric(names(tmp))
#     Nr <- as.numeric(unname(tmp))
#   }
#
#   N <- sum(Nr)
#   m1 <- sum(Nr * r) / N
#   m2 <- sum(Nr * r * (r-1)) / N
#   m3 <- sum(Nr * r * (r-1) * (r-2)) / N
#
#   nu <- max(0, (2 * m2^2 - m1*m3) / (m1*m3 - m2^2))
#   sigma <- (m2 / m1) / (nu + 1)
#   pi <- min(1, max(1 - (m1 / (nu * sigma)), 0))
#
#   out <- c(nu = nu, sigma = sigma, pi = pi)
#
#   if (return_r_Nr) {
#     attr(out, "r") <- r
#     attr(out, "Nr") <- Nr
#   }
#
#   out
#
# }
#
#
# # ---------------------------------------------------------------
# # ML Estimation of ZIPoissonGamma parameters
# # ---------------------------------------------------------------
# mle_dzipoisgamma <- function(r = NULL,
#                              Nr = NULL,
#                              x = NULL,
#                              method = "BFGS",
#                              ...) {
#
#   # initialize with method of moments estimates
#   init_mm <- mme_dzipoisgamma(
#     r = r,
#     Nr = Nr,
#     x = x,
#     return_r_Nr = TRUE
#   )
#
#   r <- attr(init_mm, "r")
#   Nr <- attr(init_mm, "Nr")
#
#
#   # optimization of pi in logit scale
#   # nu, sigma in log scale
#   expit <- function(x) {
#     n <- length(x)
#     idx_pos <- x >= 0
#     n_idx_pos <- sum(idx_pos)
#     out <- rep(NA, n)
#     if (n_idx_pos > 0) {
#       out[idx_pos] <- 1/(1 + exp(-x[idx_pos]))
#     }
#     if (n_idx_pos < n) {
#       tmp <- exp(x[!idx_pos])
#       out[!idx_pos] <- tmp/(1+tmp)
#     }
#     setNames(out, names(x))
#   }
#
#   logit <- function(x) {
#     log(x) - log(1 - x)
#   }
#
#   neg_llik <- function(par) {
#     nu <- exp(par["nu"])
#     sigma <- exp(par["sigma"])
#     pi <- expit(par["pi"])
#
#     tmp <- Nr * dzipoisgamma(
#       x = r,
#       nu = nu,
#       sigma = sigma,
#       pi = pi,
#       log = TRUE
#     )
#     -sum(tmp)
#   }
#
#   scaled_init_mm <- c(
#     log(init_mm["nu"]),
#     log(init_mm["sigma"]),
#     logit(init_mm["pi"])
#   )
#   scaled_init_mm[is.infinite(scaled_init_mm)] <- 0
#
#   opt_list <- lapply(
#     list(
#       scaled_init_mm,
#       c(nu = 0, sigma = 0, pi = 0)
#     ),
#     function(this_init) {
#       optim(
#         par = this_init,
#         fn = neg_llik,
#         method = method,
#         ...
#       )
#     }
#   )
#
#   obj_val <- sapply(opt_list, "[[", "value")
#
#   opt <- opt_list[[which.min(obj_val)[1]]]
#
#   out <- c(
#     exp(opt$par["nu"]),
#     exp(opt$par["sigma"]),
#     expit(opt$par["pi"])
#   )
#   names(opt$par) <- c("log_nu", "log_sigma", "logit_pi")
#   attr(out, "optim") <- opt
#
#   out
# }
#
# # ---------------------------------------------------------------
#

# MLE of parameters in n_ij ~ ZIP(lambda_ij * E_ij, omega)
.estimate_zip_mle <- function(n_ij, E_ij, ...) {
  expit <- function(x) {
    n <- length(x)
    idx_pos <- x >= 0
    n_idx_pos <- sum(idx_pos)
    out <- rep(NA, n)
    if (n_idx_pos > 0) {
      out[idx_pos] <- 1/(1 + exp(-x[idx_pos]))
    }
    if (n_idx_pos < n) {
      tmp <- exp(x[!idx_pos])
      out[!idx_pos] <- tmp/(1+tmp)
    }
    setNames(out, names(x))
  }

  neg_llik <- function(par) {
    omega <- expit(par["logit_omega"])
    lambda_ij <- exp(par[names(par) != "logit_omega"])
    lambda_ij_times_E_ij <- c(lambda_ij) * c(E_ij)

    llik_contri <- dzipois(
      x = c(n_ij),
      lambda = lambda_ij_times_E_ij,
      pi = omega,
      log = TRUE
    )

    -sum(llik_contri)
  }

  par_init <- c(
    logit_omega = 0,
    log_lambda = log(pmax(c(n_ij/pmax(E_ij, 1e-10)), 1e-10))
  )

  # start with BFGS
  opt <- tryCatch(
    optim(
      par = par_init,
      fn = neg_llik,
      method = "BFGS",
      ...
    ),
    error = function(e) e
  )

  if (is(opt, "error")) {
    # try Nelder-Mead
    opt <- tryCatch(
      optim(
        par = par_init,
        fn = neg_llik,
        method = "Nelder-Mead",
        ...
      ),
      error = function(e) e
    )

    if (is(opt, "error")) stop(opt)
  }

  # opt_list <- sapply(
  #   c("BFGS", "Nelder-Mead"),
  #   function(this_method) {
  #     out <- tryCatch(
  #       optim(
  #         par = par_init,
  #         fn = neg_llik,
  #         method = this_method,
  #         # lower = c(nu = 1e-5, sigma = 1e-7, pi = 0),
  #         # upper = c(nu = 1e7, sigma = 1e7, pi = 1-1e-7),
  #         ...
  #       ),
  #       error = function(e) e
  #     )
  #
  #     if (is(out, "error")) {
  #       out$val <- Inf
  #       out$par <- NA
  #     }
  #
  #     out
  #   },
  #   simplify = FALSE,
  #   USE.NAMES = TRUE
  # )


  # return the optimization instance
  # that reached minimum across all 27 instances
  # obj_val <- sapply(opt_list, "[[", "value")
  # opt <- opt_list[[which.min(obj_val)[1]]]
  est_omega <- opt$par %>%
    .["logit_omega"] %>%
    expit() %>%
    unname()
  est_lambda <- opt$par %>%
    .[names(.) != "logit_omega"] %>%
    exp() %>%
    unname() %>%
    matrix(nrow(n_ij), ncol(n_ij)) %>%
    `dimnames<-`(dimnames(n_ij))

  out <- c(
    omega = est_omega,
    lambda = est_lambda,
    llik = -opt$value,
    optim = opt
  )

  out
}



# ------------------------------------------------------------------
# ML Estimation of parameters in the ZIGammaPoisson model
# n_ij ~ ZIP(lambda_ij * E_ij, omega), lambda_ij ~ Gamma(nu, sigma)
# ------------------------------------------------------------------
.estimate_zigammapois_mle <- function(n_ij,
                                      E_ij,
                                      method = "BFGS",
                                      ...) {
  expit <- function(x) {
    n <- length(x)
    idx_pos <- x >= 0
    n_idx_pos <- sum(idx_pos)
    out <- rep(NA, n)
    if (n_idx_pos > 0) {
      out[idx_pos] <- 1/(1 + exp(-x[idx_pos]))
    }
    if (n_idx_pos < n) {
      tmp <- exp(x[!idx_pos])
      out[!idx_pos] <- tmp/(1+tmp)
    }
    setNames(out, names(x))
  }

  E_ij_adj <- pmax(E_ij, 1e-20)

  neg_llik <- function(par) {
    nu <- exp(par["log_nu"])
    # nu <- exp(par["nu"])
    # beta <- par["beta1"]

    sigma <- exp(par["log_sigma"])
    # sigma <- par["sigma"]
    sigma_k <- sigma * c(E_ij_adj)

    omega <- expit(par["logit_omega"])
    # pi <- par["pi"]

    tmp <- dzipoisgamma(
      x = c(n_ij),
      nu = nu,
      sigma = sigma_k,
      pi = omega,
      log = TRUE
    )

    -sum(tmp)
  }

  par_init <- c(
    logit_omega = 0,
    log_nu = 0,
    log_sigma = 0
  )

  # start with BFGS
  opt <- tryCatch(
    optim(
      par = par_init,
      fn = neg_llik,
      method = "BFGS",
      ...
    ),
    error = function(e) e
  )

  if (is(opt, "error")) {
    # try Nelder-Mead
    opt <- tryCatch(
      optim(
        par = par_init,
        fn = neg_llik,
        method = "Nelder-Mead",
        ...
      ),
      error = function(e) e
    )

    if (is(opt, "error")) stop(opt)
  }


  est_omega <- opt$par %>%
    .["logit_omega"] %>%
    expit() %>%
    unname()
  est_nu <- opt$par %>%
    .["log_nu"] %>%
    expit() %>%
    unname()
  est_sigma <- opt$par %>%
    .["log_sigma"] %>%
    expit() %>%
    unname()

  out <- c(
    omega = est_omega,
    nu = est_nu,
    sigma = est_sigma,
    llik = -opt$value,
    optim = opt
  )

  out
}

# .estimate_zigammapois_mle <- function(x,
#                                       a = rep(1, length(x)),
#                                       method = "BFGS",
#                                       fast_init = FALSE,
#                                       zero_inflation = TRUE,
#                                       init_pi = 0.5,
#                                       ...) {
#
#   expit <- function(x) {
#     n <- length(x)
#     idx_pos <- x >= 0
#     n_idx_pos <- sum(idx_pos)
#     out <- rep(NA, n)
#     if (n_idx_pos > 0) {
#       out[idx_pos] <- 1/(1 + exp(-x[idx_pos]))
#     }
#     if (n_idx_pos < n) {
#       tmp <- exp(x[!idx_pos])
#       out[!idx_pos] <- tmp/(1+tmp)
#     }
#     setNames(out, names(x))
#   }
#
#   # logit <- function(x) {
#   #   log(x) - log(1 - x)
#   # }
#   if (zero_inflation) {
#     neg_llik <- function(par) {
#       nu <- exp(par["log_nu"])
#       # nu <- exp(par["nu"])
#       beta <- par["beta1"]
#
#       sigma <- exp(par["log_sigma"])
#       # sigma <- par["sigma"]
#       sigma_k <- sigma * a^beta
#
#       pi <- expit(par["logit_pi"])
#       # pi <- par["pi"]
#
#       tmp <- dzipoisgamma(
#         x = x,
#         nu = nu,
#         sigma = sigma_k,
#         pi = pi,
#         log = TRUE
#       )
#
#       -sum(tmp)
#     }
#   } else {
#     neg_llik <- function(par) {
#       nu <- exp(par["log_nu"])
#       # nu <- exp(par["nu"])
#
#       sigma <- exp(par["log_sigma"])
#       # sigma <- par["sigma"]
#
#       beta <- par["beta1"]
#
#       sigma_k <- sigma * a^beta
#
#       tmp <- dzipoisgamma(
#         x = x,
#         nu = nu,
#         sigma = sigma_k,
#         pi = 0,
#         log = TRUE
#       )
#
#       -sum(tmp)
#     }
#   }
#
#   par_init_list <-
#     if (fast_init) {
#       cbind(log_nu = c(-2, 0, 2),
#             # sigma = exp(c(-2, 0, 2)),
#             log_sigma = c(-2, 0, 2),
#             # pi = expit(c(-2, 0, 2))
#             logit_pi = c(-2, 0, 2),
#             beta1 = c(-5, 0, 5))
#     } else {
#       as.matrix(
#         expand.grid(
#           # nu = exp(c(-2, 0, 2)),
#           log_nu = c(-2, 0, 2),
#           # sigma = exp(c(-2, 0, 2)),
#           log_sigma = c(-2, 0, 2),
#           # pi = expit(c(-2, 0, 2))
#           logit_pi = c(-2, 0, 2),
#           beta1 = c(-5, 0, 5)
#         )
#       )
#     }
#
#   if (!zero_inflation) {
#     par_init_list <- par_init_list[, c("log_nu", "log_sigma", "beta1")]
#   }
#
#   # browser()
#
#   # constrained optimization with all 3^3 initial points
#   opt_list <- apply(
#     par_init_list, 1,
#     function(this_init) {
#       out <- tryCatch(
#         optim(
#           par = this_init,
#           fn = neg_llik,
#           method = method,
#           # lower = c(nu = 1e-5, sigma = 1e-7, pi = 0),
#           # upper = c(nu = 1e7, sigma = 1e7, pi = 1-1e-7),
#           ...
#         ),
#         error = function(e) e
#       )
#       if (is(out, "error")) {
#         # use Nelder-Mead (doesn't use gradient)
#         out <- tryCatch(
#           optim(
#             par = this_init,
#             fn = neg_llik,
#             # method = method,
#             # lower = c(nu = 1e-5, sigma = 1e-7, pi = 0),
#             # upper = c(nu = 1e7, sigma = 1e7, pi = 1-1e-7),
#             ...
#           ),
#           error = function(e) e
#         )
#       }
#
#       if (is(out, "error")) browser()
#
#       out
#     }
#   )
#
#
#   # return the optimization instance
#   # that reached minimum across all 27 instances
#   obj_val <- sapply(opt_list, "[[", "value")
#   opt <- opt_list[[which.min(obj_val)[1]]]
#   est_pi <- 0
#   if (zero_inflation) {
#     est_pi <- unname(expit(opt$par["logit_pi"]))
#   }
#   # out <- opt_par
#   out <- c(
#     nu = unname(exp(opt$par["log_nu"])),
#     sigma = unname(exp(opt$par["log_sigma"])),
#     pi = est_pi,
#     beta1 = unname(opt$par["beta1"])
#   )
#   attr(out, "optim") <- opt
#   attr(out, "loglik") <- -opt$value
#
#   out
# }
# # ---------------------------------------------------------------
#
#
#
# # ---------------------------------------------------------------
# # Empirical Bayes estimates of theta
# # given x_k, a_k, nu, sigma, pi
# # ---------------------------------------------------------------
#
# compute_bayes_est_theta <- function(x, a, nu, sigma, pi, beta1 = 1) {
#   sigma_k <- sigma * a^beta1
#   delta_k <- pi / (
#     pi +
#       (1-pi) / ((1 + sigma_k)^nu)
#   )
#   theta_eb <- ifelse(
#     x == 0,
#     # xk = 0
#     delta_k * nu * sigma_k +
#       (1-delta_k) * (x + nu)/(1 + 1/sigma_k),
#     # xk > 0
#     (x + nu)/(1 + 1/sigma_k)
#   )#/a_k
#
#   theta_eb
# }
