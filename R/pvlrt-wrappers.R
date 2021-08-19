#' Pseudo Likelihood Ratio Test under the
#' zero-inflated Poisson model with relative reporting rate parametrization
#' @inheritParams pvlrt
#' @param ... additional arguments passed to pvlrt
#'
#'
#'
#' @note
#' \code{lrt_zi_poisson()} is a wrapper for \code{pvlrt()} with
#' \code{parametrization = "rrr"}.
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
                           ...)
{
  dots <- list(...)
  if (!is.null(dots$parametrization)) {
    msg <- glue::glue(
      "parametrization' is automatically set to 'rrr' in \\
      `lrt_zi_poisson()`. For 'rr' parametrization use \\
      pvlrt() or lrt_poisson()"
    )
    warning(msg)
  }

  dots$parametrization <- "rrr"
  inargs <- list(
    contin_table = contin_table,
    nsim = nsim
  ) %>%
    c(dots)

  do.call(pvlrt, inargs)
}



#' Likelihood Ratio Test under the (vanilla, non-zero-inflated) Poisson model
#' @inheritParams pvlrt
#' @inheritParams lrt_zi_poisson
#'
#' @examples
#'
#' data("lovastatin")
#'
#' # no grouping -- each drug forms its own class
#' test1 <- lrt_poisson(lovastatin)
#' ## extract the observed LRT statistic
#' attr(test1, "lrstat")
#'
#'
#' # grouped drugs --
#' # group1 : drug 1, drug 2
#' # group 2: drug 3
#' drug_groups <- list(c(1, 2), 3)
#' test2 <- lrt_vanilla_poisson(lovastatin, drug_class_idx = drug_groups)

#' @export
lrt_poisson <- function(contin_table,
                        nsim = 1e4,
                        parametrization = "rrr",
                        ...)
{

  dots <- list(...)

  for (nm in c("zi_prob", "omega_vec", "omega_est_vec")) {
    if (!is.null(dots[[nm]])) {
      msg <- glue::glue(
        "{nm}' is automatically set to 0 in `lrt_poisson()`. \\
        For non-zero zi probability use `pvlrt()` or `lrt_zi_poisson()` \\
        and supply as '{nm}'"
      )
      warning(msg)
      dots$omega_vec <- dots$omega_est_vec <- dots$zi_prob <- NULL
    }
  }

  if (parametrization %in% c("rrr", "lambda")) {
    dots$omega_vec <- rep(0, ncol(contin_table))
  }

  inargs <- list(
    contin_table = contin_table,
    nsim = nsim,
    parametrization = parametrization
  ) %>%
    c(dots)

  do.call(pvlrt, inargs)
}

#' @rdname lrt_poisson
#' @export
lrt_vanilla_poisson <- lrt_poisson
