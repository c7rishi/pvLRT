#' Pseudo Likelihood Ratio Test under the
#' zero-inflated Poisson model with relative reporting rate parametrization
#' @inheritParams pvlrt
#' @param ... additional arguments passed to pvlrt
#'
#' @seealso \link{pvlrt}
#'
#' @note
#' \code{lrt_zi_poisson()} is a wrapper for \code{pvlrt()} with
#' \code{parametrization = "rrr"}.
#'
#' @examples
#'
#' data("lovastatin")
#' test1 <- lrt_zi_poisson(lovastatin)
#' test1
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
#' @seealso \link{pvlrt}
#' @note
#' \code{lrt_poisson()} and \code{lrt_vanilla_poisson()}
#' are both wrappers for \code{pvlrt()} with
#' \code{omega_vec = rep(0, ncol(contin_table))}
#'
#'
#' @examples
#'
#' data("lovastatin")
#'
#' # no grouping -- each drug forms its own class
#' test1 <- lrt_poisson(lovastatin)
#'
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
