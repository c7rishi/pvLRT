is.pvlrt <- function(obj,...) {
  is (obj, "pvlrt")
}

extract_n_matrix <- function(object, ...) {
  attr(object, "contin_table")
}


#' @export
summary.pvlrt <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }

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

  tab_lrstat_pval
}


#' @export
print.pvlrt <- function(object, significance_level = 0.05,
                        topn = 12, digits = 2, ...) {
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
      paste0("q", .)

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

  if (topn > 1) {
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


#' @export
as.matrix.pvlrt <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }
  class(object) <- "matrix"
  object
}



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
