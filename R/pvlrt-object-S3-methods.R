is.pvlrt <- function(obj,...) {
  is (obj, "pvlrt")
}

extract_n_matrix <- function(object, ...) {
  attr(object, "contin_table") %>%
    set_dimnames(dimnames(object))
}


summary_zi_probs <- function(object, ...) {
  all_attr <- attributes(object)

  out <- data.table(
    AE = colnames(object),
    zi = all_attr$omega
  )

  if (all_attr$test_omega) {
    zi_lrstat <- zi_p.value <- zi_q.value <- NULL
    out[
      ,
      `:=`(
        lrstat = all_attr$omega_lrstat,
        p.value = all_attr$omega_pvalue,
        q.value = all_attr$omega_qvalue
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
#' @examples
#' \dontrun{
#' test1 <- pvlrt(statin46, test_zi = TRUE)
#' summary(test1)
#' tmp <- summary(test1, show_zi = TRUE)
#' attr(tmp, "zi")
#' }
#' @export
summary.pvlrt <- function(object, show_zi = FALSE, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }

  stopifnot(
    is.logical(show_zi),
    length(show_zi) == 1
  )

  tab_lrstat_pval <- list(
    n = extract_n_matrix,
    lrstat = extract_lrstat_matrix,
    p.value = extract_pvalue_matrix
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
#' @inheritParams extract_pvalue_matrix
#' @param topn number of top (with respect to likelihood ratio statistic value)
#' pairs to show at the given significance level.
#' @param digits number of digits to show after the decimal place.
#' @export
print.pvlrt <- function(object,
                        significance_level = 0.05,
                        topn = 12,
                        digits = 2,
                        ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }

  all_attr <- attributes(object)
  omega <- all_attr$omega
  do_omega_estimation <- all_attr$do_omega_estimation
  parametrization <- ifelse(
    all_attr$parametrization == "rrr",
    "Relative reporting rate (lambda)",
    "Reporting rate (p-q)"
  )

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
    )
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
      Pseudo LR test of significance for zi probabilities:
      {zi_test_text}"
    )
  }

  I <- nrow(all_attr$lrstat)
  J <- ncol(all_attr$lrstat)
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
      .[, c("AE", "Drug", "lrstat", "p.value"), with = FALSE]

    top_signif_pairs$p.value <- top_signif_pairs$p.value %>%
      format.pval(digits = 3, eps = 1e-5)

    res_signif_pairs <- top_signif_pairs %>%
      print() %>%
      capture.output() %>%
      paste(collapse = "\n")

    signif_pairs_txt <- glue::glue(
      "Total {n_signif_pair} AE-drug \\
      {ifelse(n_signif_pair > 1, 'pairs are', 'pair is')} significant \\
      at level = {significance_level}.
      (FDR controlled) p-values for {ifelse(topn == n_signif_pair, '',  'top')} {topn} \\
      significant pairs:
      {res_signif_pairs}"
    )
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

  msg <- glue::glue(
    "{parametrization}-based {lrt_type}-LRT on \\
    {I} AE & {J} drugs.
    Hypothesis tests performed on {n_drug_test_idx} \\
    {ifelse(n_drug_test_idx > 1, 'drugs', 'drug')} using \\
    {null_boot_type} bootstrap.

    {zi_text}

    {signif_pairs_txt}

    Extract all LR statistics and p-values using `summary()`
    "
  )
  cat(msg)
}


#' @inheritParams summary.pvlrt
#' @export
as.matrix.pvlrt <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  class(object) <- "matrix"
  object
}



#' Plotting method for a pvlrt object
#'
#' @inheritParams summary.pvlrt
#' @param type character string determining the type of plot to show.
#' Available choices are "heatmap" (default) which calls \link{heatmap_pvlrt},
#' and "barplot" which calls \link{barplot.pvlrt} with the additional arguments
#' supplied in ...
#' @param ... additional arguments passed to heatmap_pvlrt or barplot.pvlrt
#' depending on \code{type}.
#' @export
plot.pvlrt <- function(object, type = "heatmap", ...) {
  stopifnot(
    type %in% c("heatmap", "barplot")
  )

  out <- if (type == "heatmap") {
    heatmap_pvlrt(object, ...)
  } else if (type == "barplot") {
    barplot(object, ...)
  }

  out
}
