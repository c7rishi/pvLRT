is.pvlrt <- function(obj,...) {
  is (obj, "pvlrt")
}

`:=` <- data.table::`:=`

#' @export
summary.pvlrt <- function(object, ...) {
  if (!is.pvlrt(object)) {
    stop("object must be a 'pvlrt' object.")
  }

  tab_lrstat_pval <- c("lrstat", "p.value") %>%
    lapply(
      function(measure) {
        object %>%
          {
            if (measure == "p.value") {
              `class<-`(., "matrix")
            } else {
              attr(., "lrstat")
            }
          } %>%
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
      }
    ) %>%
    Reduce(
      function(dt1, dt2) {
        merge(dt1, dt2, by = c("AE", "Drug"))
      },
      .
    ) %>%
    data.table::setDT() %>%
    data.table::setorder(p.value)

  tab_lrstat_pval
}


#' @export
print.pvlrt <- function(object, significance_level = 0.05, topn = 12, digits = 4, ...) {
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
      "Drug-specific {zi_est_type} zi probabilities: {zi_values}."
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
      .[, c("AE", "Drug", "p.value"), with = FALSE]

    top_signif_pairs$p.value <- top_signif_pairs$p.value %>%
      format.pval(digits = 3, eps = 1e-5)

    res_signif_pairs <- top_signif_pairs %>%
      print() %>%
      capture.output() %>%
      paste(collapse = "\n")

    signif_pairs_txt <- glue::glue(
      "Total {n_signif_pair} \\
      {ifelse(n_signif_pair > 1, 'pairs are', 'pair is')} significant\\
      at level = {significance_level}.
      p-values for {ifelse(topn == n_signif_pair, '',  'top')} significant pairs:
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

  msg <- glue::glue(
    "{parametrization}-based {lrt_type}-LRT on \\
    {I} AE & {J} drugs input table \\
    (testing performed on {n_drug_test_idx} \\
    {ifelse(n_drug_test_idx > 1, 'drugs', 'drug')}).

    {zi_text}

    {signif_pairs_txt}

    Extract all LR statistics and p.values using `summary()`
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
