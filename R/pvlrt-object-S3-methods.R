is.pvlrt <- function(obj, ...) {
  is(obj, "pvlrt")
}


summary_zi_probs <- function(object, ...) {
  all_attr <- attributes(object)

  out <- data.table(
    AE = colnames(object),
    zi = all_attr$omega
  )
  . <- NULL

  if (all_attr$test_omega) {
    zi_lrstat <- zi_p_value <- zi_q_value <- NULL
    out[
      ,
      `:=`(
        lrstat = all_attr$omega_lrstat,
        p_value = all_attr$omega_p_value,
        q_value = all_attr$omega_qvalue
      )
    ]
  }

  out
}

#' Summary method for a pvlrt object
#'
#' @param object a \code{pvlrt} object, which is the output of the function
#' \link{pvlrt} or one of its wrappers such as \link{lrt_zi_poisson},
#' \link{lrt_poisson} and \code{lrt_vanilla_poisson}.
#' @param ... other input parameters. Currently unused.
#' @param show_zi logical. Should summary of the estimates and tests
#' (if performed) of the zero inflation parameters be returned?
#' Defaults to FALSE. If TRUE, then the zero inflation summary is included
#' as an attribute with name "zi". See examples.
#'
#' @return
#' Returns a data.table with rows corresponding to all possible AE/Drug pairs
#' as obtained from the input contingency table, and columns titled "AE", "Drug", "n", "lrstat" (log-likelihood ratio
#' test statistic) and "p_value". Additionally, if `show_zi` is set to `TRUE`, then as
#' an attribute named "zi" a data.table with rows
#' corresponding to Drugs (columns in the input contingency table), and columns titled
#' "AE", "zi", "lrstat" (log-likelihood ratio test statistic for zero-inflation),
#' "p_value" and "q_value" (Benjamini-Hochberg adjusted p-values, as obtained through
#' \link[stats]{p.adjust}) is returned.
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, test_zi = TRUE, nsim = 500)
#' summary(test1)
#' tmp <- summary(test1, show_zi = TRUE)
#' print(tmp)
#' tmp_zi <- attr(tmp, "zi")
#' print(tmp_zi)
#'
#' @seealso
#' \link{pvlrt}
#'
#' @md
#' @export
summary.pvlrt <- function(object, show_zi = FALSE, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  . <- NULL
  lrstat <- NULL

  stopifnot(
    is.logical(show_zi),
    length(show_zi) == 1
  )

  tab_lrstat_pval <- list(
    n = extract_n_matrix,
    lrstat = extract_lrstat_matrix,
    p_value = extract_p_value_matrix
  ) %>%
    mapply(
      function(measure_fn, measure) {
        measure_fn(object) %>%
          data.matrix() %>%
          data.table::as.data.table(
            keep.rownames = TRUE
          ) %>%
          data.table::setnames(
            old = "rn",
            new = "AE"
          ) %>%
          data.table::melt(
            id.vars = "AE",
            variable.name = "Drug",
            value.name = measure
          )
      },
      .,
      names(.),
      SIMPLIFY = FALSE
    ) %>%
    Reduce(
      function(dt1, dt2) {
        merge(dt1, dt2, by = c("AE", "Drug"))
      },
      .
    ) %>%
    # using base pipe instead of magrittr pipe
    # would be faster, I think
    data.table::setDT() %>%
    data.table::setorder(-lrstat)

  if (show_zi) {
    summ_zi <- summary_zi_probs(object)
    attr(tab_lrstat_pval, "zi") <- summ_zi
  }

  tab_lrstat_pval
}

#' Print method for pvlrt objects
#' @inheritParams summary.pvlrt
#' @inheritParams extract_p_value_matrix
#' @param x a \code{pvlrt} object; an output of function \code{pvlrt}().
#' @param topn number of top (with respect to likelihood ratio statistic value)
#' pairs to show at the given significance level.
#' @param digits number of digits to show after the decimal place.
#' @param show_test_summary logical. Should a brief summary showing the top few
#' test results be displayed? defaults to FALSE.
#'
#' @return
#' Invisibly returns the input `pvlrt` object.
#'
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, nsim = 500)
#' print(test1)
#'
#' @seealso
#' \link{pvlrt}
#'
#' @md
#' @export
print.pvlrt <- function(x,
                        significance_level = 0.05,
                        topn = 12,
                        digits = 2,
                        show_test_summary = FALSE,
                        ...) {
  object <- x
  if (!is.pvlrt(object)) {
    stop("x must be a 'pvlrt' object.")
  }
  . <- NULL
  lrstat <- NULL


  all_attr <- attributes(object)
  omega <- all_attr$omega
  do_omega_estimation <- all_attr$do_omega_estimation
  parametrization <- ifelse(
    all_attr$parametrization == "rrr",
    "Relative reporting rate (lambda)",
    "Reporting rate (p-q)"
  )
  stopifnot(
    is.numeric(significance_level),
    is.numeric(topn),
    is.numeric(digits),
    is.logical(show_test_summary)
  )

  run_time_txt <- extract_run_time(object) %>%
    capture.output() %>%
    gsub("Time difference of ", "", .) %>%
    paste("Running time of the original pvlrt call:", .)


  if (!do_omega_estimation & all(omega == 0)) {
    zi_text <- "No zi considered in the model."
  } else {
    zi_est_type <- ifelse(do_omega_estimation, "estimated", "supplied")
    zi_values <- paste0(
      round(omega, digits),
      " (", all_attr$dimnames[[2]], ")"
    ) %>%
      paste(collapse = ", ")

    zi_text <- glue::glue(
      "Drug-specific {zi_est_type} zi probabilities:
      {zi_values}"
    ) %>%
      strwrap(prefix = "\n", initial = "") %>%
      paste(collapse = "")
  }

  is_zi_test <- all_attr$test_omega
  if (is_zi_test) {
    zi_lr_stat <- all_attr$omega_lrstat %>%
      round(digits) %>%
      format(nsmall = digits) %>%
      paste0("LRstat=", ., ", ")

    zi_q_values <- all_attr$omega_qval %>%
      format_pval_(digits = digits) %>%
      paste0(
        " (", all_attr$dimnames[[2]], ")"
      ) %>%
      paste0("q", .) %>%
      gsub("=<", "<", .) %>%
      gsub("qNA", "q=<NA>", .)

    zi_test_text <- paste0(zi_lr_stat, zi_q_values) %>%
      paste(collapse = ", ")

    zi_text <- glue::glue(
      "{zi_text} \n
       LR test of significance for zi probabilities:
      {zi_test_text}"
    ) %>%
      strwrap(prefix = "\n", initial = "") %>%
      paste(collapse = "")
  }

  I <- nrow(object)
  J <- ncol(object)
  n_drug_test_idx <- length(all_attr$test_drug_idx)

  signif_pairs <- extract_significant_pairs(
    object,
    significance_level = significance_level
  )

  n_signif_pair <- nrow(signif_pairs)

  topn <- min(n_signif_pair, topn)

  if (topn >= 1) {
    top_signif_pairs <- signif_pairs %>%
      .[1:topn, ] %>%
      .[, c("AE", "Drug", "lrstat", "p_value"), with = FALSE]

    top_signif_pairs$p_value <- top_signif_pairs$p_value %>%
      format.pval(digits = 3, eps = 1e-5)

    res_signif_pairs <- top_signif_pairs %>%
      print() %>%
      capture.output() %>%
      paste(collapse = "\n")

    signif_pairs_txt <- if (show_test_summary) {
      glue::glue(
        "Total {n_signif_pair} AE-drug \\
      {ifelse(n_signif_pair > 1, 'pairs are', 'pair is')} significant \\
      at level = {significance_level}.
      (FDR controlled) p-values for {ifelse(topn == n_signif_pair, '',  'top')} {topn} \\
      significant pairs:
      {res_signif_pairs}"
      )
    } else {
      glue::glue(
        "Total {n_signif_pair} AE-drug \\
      {ifelse(n_signif_pair > 1, 'pairs are', 'pair is')} significant \\
      at level = {significance_level}."
      )
    }
  } else {
    signif_pairs_txt <- glue::glue(
      "No significant pair at level = {significance_level}."
    )
  }

  lrt_type <- ifelse(
    do_omega_estimation &
      all_attr$parametrization == "rrr",
    "pseudo",
    "ordinary"
  )

  res_tab <- capture.output(print(signif_pairs)) %>%
    paste(collapse = "\n")

  null_boot_type <- attr(object, "null_boot_type")

  top_text <- glue::glue(
    "{parametrization}-based {lrt_type}-LRT on \\
    {I} AE & {J} drugs.
    Hypothesis tests performed on {n_drug_test_idx} \\
    {ifelse(n_drug_test_idx > 1, 'drugs', 'drug')} using \\
    {null_boot_type} bootstrap."
  ) %>%
    strwrap(prefix = "\n", initial = "") %>%
    paste(collapse = "")

  msg <- glue::glue(
    "{top_text}

    {run_time_txt}

    {zi_text}
    {signif_pairs_txt}

    Extract all LR statistics and p-values using `summary()`.
    "
  )
  cat(msg)
  invisible(object)
}

#' Casting a `pvlrt` object as a matrix of log LR statistics
#'
#' @inheritParams summary.pvlrt
#' @inheritParams print.pvlrt
#' @method as.matrix pvlrt
#'
#' @return
#' Returns a matrix with the same dimensions as the input contingency
#' table in the original `pvlrt` call, with each cell providing
#' the corresponding value of the observed log-likelihood ratio
#' test statistic.
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, nsim = 500)
#' as.matrix(test1)
#'
#' @seealso
#' \link{pvlrt}
#'
#' @md
#' @export
as.matrix.pvlrt <- function(x, ...) {
  object <- x
  if (!is.pvlrt(object)) {
    stop("x must be a 'pvlrt' object.")
  }
  attributes(object) <- attributes(
    object
    )[c("dim", "dimnames", "class")]
  class(object) <- "matrix"
  object
}



#' Plotting method for a pvlrt object
#'
#' @inheritParams summary.pvlrt
#' @inheritParams print.pvlrt
#' @param type character string determining the type of plot to show.
#' Available choices are `"bubbleplot"` which calls \link{bubbleplot_pvlrt},
#' `"heatmap"` which calls \link{heatmap_pvlrt}, and
#' `"barplot"` which calls \link{barplot.pvlrt}, with the additional arguments
#' supplied in ...
#' @param ... additional arguments passed to heatmap_pvlrt or barplot.pvlrt
#' depending on \code{type}.
#'
#' @return
#' A \link[ggplot2]{ggplot} object.
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' test1 <- pvlrt(statin46, nsim = 500)
#' plot(test1, type = "bubbleplot")
#' plot(test1, type = "barplot")
#' plot(test1, type = "heatmap")
#'
#' @seealso
#' \link{pvlrt}; \link{bubbleplot_pvlrt}; \link{heatmap_pvlrt}; \link{barplot.pvlrt}
#'
#' @md
#' @export
plot.pvlrt <- function(x, type = "bubbleplot", ...) {
  object <- x
  . <- NULL

  if (!is.pvlrt(object)) {
    stop("x must be a 'pvlrt' object.")
  }

  stopifnot(
    type %in% c("heatmap", "barplot", "bubbleplot")
  )

  out <- if (type == "heatmap") {
    heatmap_pvlrt(object, ...)
  } else if (type == "barplot") {
    barplot(object, ...)
  } else if (type == "bubbleplot") {
    bubbleplot_pvlrt(object, ...)
  }

  out
}



#' Overall Log-likelihood for a pvlrt object
#' @inheritParams summary.pvlrt
#' @param type Type of model and hypothesis combination.
#' Available choices are "full-poisson", "null-poisson", "full-zip" (default), and
#' "null-zip". See details.
#'
#' @note
#' The overall log likelihood must be computed during the original pvlrt run. This is
#' ensured by setting  \code{return_overall_loglik = TRUE}, and
#' \code{parametrization = "lambda"} (or \code{parametrization = "rrr"}) while running pvlrt().
#'
#' @details
#' The function extracts the overall log-likelihood and degrees of freedom
#' (the number of estimated parameters) from a \code{pvlrt} object. There are
#' four possible choices of distribution and hypothesis combinations, supplied
#' through the argument \code{type}, with the default being \code{type = "full-zip"}.
#' In a "full" model the signal parameters lambdas are estimated for all cells
#' in the contingency table from the data (subject to the condition lambda >= 1), whereas under a "null"
#' model each lambda is fixed at 1 for each cell. In a "zip" model
#' (type = "full-zip" and type = "null-zip") the log-likelihood under a zero-inflated
#' Poisson model with estimated or supplied zero inflation parameters (
#' through \code{zi_prob} in \link{pvlrt}) is returned. The degrees of freedom
#' reflects whether the zero-inflation parameters are estimated. Note that if
#' an ordinary Poisson LRT is run either by setting \code{zi_prob = 0} in
#' \link{pvlrt} or equivalently through \link{lrt_poisson} then "full-zip" and
#' "null-zip" refer to zero-inflated poisson models with presepecified
#' zero-inflation probabilities equal to 0. Thus, in such cases the results
#' with type = "full-zip" and type =  "null-zip" are identical to
#' type = "full-poisson" and type = "null-poisson"
#' respectively. See examples.
#'
#' @return
#' An object of class \link[stats]{logLik}. See Details.
#'
#' @examples
#'
#' # 500 bootstrap iterations (nsim) in each example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#'
#' set.seed(100)
#' # estimates zero inflation probabilities
#' test1 <- pvlrt(statin46, nsim = 500)
#' logLik(test1)
#' AIC(test1)
#' BIC(test1)
#'
#' # compare with and without zero inflation
#' BIC(logLik(test1, type = "full-zip"))
#' BIC(logLik(test1, type = "full-poisson"))
#'
#' # ordinary poisson model
#' ## equivalent to pvlrt(statin46, zi_prob = 0, nsim = 500)
#' test2 <- lrt_poisson(statin46, nsim = 500)
#'
#' all.equal(logLik(test2, "full-zip"), logLik(test2, "full-poisson"))
#'
#' @seealso
#' \link{pvlrt}; \link[stats]{AIC}
#'
#' @md
#' @export
logLik.pvlrt <- function(object, type = "full-zip", ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  . <- NULL

  check_loglik <- !is.null(attr(object, "return_overall_loglik"))
  if (check_loglik) {
    check_loglik <- attr(object, "return_overall_loglik")
  }
  if (!check_loglik) {
    msg <- glue::glue(
      "Overall log-likelihood not computed. Rerun `pvlrt()` with
      return_overall_loglik = TRUE"
    )
    stop(msg)
  }

  if (!attr(object, "parametrization") %in% c("rrr", "lambda")) {
    stop("parametrization in pvlrt must be 'rrr' or 'lambda'")
  }

  possible_types <- c(
    "full-poisson", "null-poisson",
    "full-zip", "null-zip"
  )
  if (!type %in% possible_types) {
    msg <- possible_types %>%
      paste0("'", ., "'") %>%
      paste(collapse = ", ") %>%
      paste("type must be one of:", .)
    stop(msg)
  }

  this_type <- type

  dat_loglik <- attr(object, "loglik_df")

  val <- dat_loglik[type == this_type]$logLik
  attr(val, "nobs") <- dat_loglik[type == this_type]$N
  attr(val, "df") <- dat_loglik[type == this_type]$df
  class(val) <- "logLik"

  val
}


#' Test if two pvlrt objects are (Nearly) Equal
#' @param target First pvlrt object (output of \link{pvlrt}).
#' @param current Second pvlrt object.
#' @param ... Arguments passed to \link{all.equal.default}.
#' @details
#' Compares all values and attributes of target and current `pvlrt` objects except running times.
#' See \link{all.equal.default} for details on the generic function.
#' @seealso all.equal.default
#' @method all.equal pvlrt
#' @export
all.equal.pvlrt <- function(target, current, ...) {
  attributes(target)$run_time <- attributes(current)$run_time <- NULL
  all.equal.default(target, current, ...)
}
