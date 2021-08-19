#' Heatmap showing LR test results
#' @param object pvlrt object
#' @param ... additional arguments passed to pheatmap::pheatmap()
#'
#' @examples
#' test1 <- pvlrt(statin46)
#' heatmap_pvlrt(test1)
#'
#' @export
heatmap_pvlrt <- function(object, nrow = NULL,
                          ncol = NULL,
                          color_by = "p.value",
                          show_n = TRUE,
                          show_pvalue = FALSE,
                          show_lrstat = FALSE,
                          arrange_alphabatical = FALSE,
                          digits = 2,
                          ...) {
  dots <- list(...)

  stopifnot(
    color_by %in% c("p.value", "lrstat")
  )

  # rows are sorted by descending LRstat
  summ <- summary(object)

  AE_all <- unique(summ$AE)
  drug_all <- unique(summ$Drug)

  if (is.null(nrow)) {
    nrow <- min(length(AE_all), 50)
  }

  if (is.null(ncol)) {
    ncol <- min(length(drug_all), 10)
  }

  AE_sub <- AE_all[1:nrow] %>%
    {if (arrange_alphabatical) sort(.) else .}

  drug_sub <- drug_all[1:ncol] %>%
    {if (arrange_alphabatical) sort(.) else .}

  measure_stat <- c("lrstat", "p.value") %>%
    setNames(., .) %>%
    lapply(
      function(measure) {
        out <- summ[
          AE %in% AE_sub &
            Drug %in% drug_sub
        ] %>%
          data.table::dcast(
            AE ~ Drug,
            value.var = measure
          ) %>%
          as.data.frame() %>%
          set_rownames(.$AE)

        out$AE <- NULL

        data.matrix(out)[AE_sub, drug_sub]
      }
    )

  p.value_text <- measure_stat$p.value %>%
    .[AE_sub, drug_sub] %>%
    format_pval_(digits = digits) %>%
    c() %>%
    paste0("p", .) %>%
    matrix(nrow, ncol)

  lrstat_text <- measure_stat$lrstat %>%
    # round(digits = digits) %>%
    .[AE_sub, drug_sub] %>%
    c() %>%
    {ifelse(
      . >= 1e4,
      format(., digits = digits, scientific = TRUE),
      as.character(round(., digits = digits))
    )} %>%
    paste0("LR=", .) %>%
    matrix(nrow, ncol)

  n_text <- attr(object, "contin_table") %>%
    .[AE_sub, drug_sub] %>%
    c() %>%
    {ifelse(
      . >= 1e4,
      format(., digits = digits, scientific = TRUE),
      as.character(round(., digits = digits))
    )} %>%
    paste0("n=", .) %>%
    matrix(nrow, ncol)

  display_text <- matrix("", nrow, ncol) %>%
    {
      if (show_n) paste(., n_text) else .
    } %>%
    {
      if (show_lrstat) paste(., lrstat_text, sep = ", ") else .
    } %>%
    {
      if (show_pvalue) paste(., p.value_text, sep = ", ") else .
    } %>%
    gsub("^\\,", "", .) %>%
    # gsub("^\\ ", "", .) %>%
    # trim white space
    trimws() %>%
    matrix(nrow, ncol)

  unused_args <- c("cluster_rows", "cluster_cols")
  for (nm in unused_args) {
    if (!is.null(dots[[nm]])) {
      msg <- glue::glue(
        "{nm} is set to FALSE when passed to `pheatmap()`"
      )
      warning(msg)
      dots[[nm]] <- NULL
    }
  }

  if (show_n | show_pvalue | show_lrstat) {
    which_show <- "" %>%
      {
        if (show_n) paste(., "show_n = TRUE", sep = ", ") else .
      } %>%
      {
        if (show_lrstat) paste(., "show_lrstat = TRUE", sep = ", ") else .
      } %>%
      {
        if (show_pvalue) paste(., "show_pvalue = TRUE", sep = ", ") else .
      } %>%
      gsub("^\\,", "", .) %>%
      # gsub("^\\ ", "", .) %>%
      # trim white space
      trimws()

    to_show <- "" %>%
      {
        if (show_n) paste(., "n's", sep = ", ") else .
      } %>%
      {
        if (show_pvalue) paste(., "p-value's", sep = ", ") else .
      } %>%
      {
        if (show_lrstat) paste(., "lrstat's", sep = ", ") else .
      } %>%
      gsub("^\\,", "", .) %>%
      # gsub("^\\ ", "", .) %>%
      # trim white space
      trimws()


    if (!is.null(dots$display_numbers)) {
      msg <- glue::glue(
        "'display_numbers' automatically set to display \\
        the {to_show} when {which_show}"
      )
      warning(msg)
      dots$display_numbers <- NULL
    }
  }


  if (is.null(dots$color)) {
    dots$color <- RColorBrewer::brewer.pal(
      n = 9, name = "Blues"
    ) %>%
      rev() %>%
      {grDevices::colorRampPalette(.)(100)}
  }

  if (is.null(dots$number_color)) {
    thresh <- measure_stat[[color_by]] %>%
      c() %>%
      median() %>%
      {./2}
    dots$number_color <- measure_stat[[color_by]] %>%
      {ifelse(
        . >= thresh,
        "black",
        "orange"
      )}
  }


  if (is.null(dots$show_rownames)) {
    dots$show_rownames <- nrow <= 50
  }

  if (is.null(dots$main)) {
    dots$main <- color_by
  }

  pheatmap_args <- list(
    mat = measure_stat[[color_by]],
    cluster_cols = FALSE,
    cluster_rows = FALSE
  ) %>%
    c(dots)

  if (show_pvalue | show_n | show_lrstat) {
    pheatmap_args$display_numbers <- display_text
  }




  do.call(pheatmap::pheatmap, pheatmap_args)

}
