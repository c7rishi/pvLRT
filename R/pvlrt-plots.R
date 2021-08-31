#' Heatmap showing LR test results
#' @param object pvlrt object
#' @param engine heatmap engine to use. Available choices are "ggplot2"
#' (default, can be abbreviated to "ggplot") and "pheatmap".
#' @param ... additional arguments passed to pheatmap::pheatmap()
#' @param fill_measure Measure to govern the filling
#' @examples
#' test1 <- pvlrt(statin46)
#' heatmap_pvlrt(test1)
#' heatmap_pvlrt(test1, arrange_alphabetical = TRUE)
#'
#' @export
heatmap_pvlrt <- function(object,
                          AE = NULL,
                          Drug = NULL,
                          grep = FALSE,
                          fill_measure = "p.value",
                          show_n = TRUE,
                          arrange_alphabetical = FALSE,
                          show_pvalue = FALSE,
                          show_lrstat = FALSE,
                          p.value_lower = 0,
                          p.value_upper = 1,
                          lrstat_lower = 0,
                          lrstat_upper = Inf,
                          n_lower = 0,
                          n_upper = Inf,
                          remove_outside = FALSE,
                          digits = 2,
                          engine = "ggplot2",
                          border_color = "white",
                          ...) {

  dots <- list(...)
  all_inputs <- c(as.list(environment()), dots)

  # processed plotting data from process_plot_data
  dat_pl <- tryCatch(
    do.call(process_plot_data, all_inputs),
    error = function(e) e
  )

  if (is(dat_pl, "error")) {
    stop(dat_pl$message)
  }

  stopifnot(
    engine %in% c("pheatmap", "ggplot", "ggplot2")
  )

  show_text <- show_pvalue | show_n | show_lrstat

  # heatmap engine specific codes

  if (grepl("ggplot", engine)) {
    dat_pl[
      ,
      AE := AE %>%
        factor(levels = rev(levels(.)))
    ]

    darkblue_col <- RColorBrewer::brewer.pal(8, "Blues") %>% tail(1)

    out <- dat_pl %>%
      ggplot2::ggplot(
        ggplot2::aes_string(
          y = "AE",
          x = "Drug",
          fill = fill_measure,
          label = "text"
        )
      ) +
      ggplot2::geom_tile(color = border_color) +
      ggplot2::scale_fill_gradient(
        low = ifelse(fill_measure == "p.value", darkblue_col, "white"),
        high = ifelse(fill_measure == "p.value", "white", darkblue_col)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90, vjust = 0.5, hjust = 1
        ),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = "", y = "")

    if (show_text) {
      out <- out +
        ggplot2::geom_text(
          ggplot2::aes_string(
            color = "text_color"
          ),
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = c(dat_pl$text_color) %>%
            unique() %>%
            rev(),
          guide = "none"
        )
    }
  } else if (engine == "pheatmap") {

    if (!requireNamespace("pheatmap")) {
      stop("package {pheatmap} is required for engine='pheatmap'")
    }

    AE_sub <- levels(dat_pl$AE)
    Drug_sub <- levels(dat_pl$Drug)

    plot_mats <- c(
      "lrstat", "p.value",
      "text", "text_color"
    ) %>%
      setNames(., .) %>%
      lapply(
        function(fill_measure) {
          out <- dat_pl %>%
            data.table::dcast(
              AE ~ Drug,
              value.var = fill_measure
            ) %>%
            as.data.frame() %>%
            set_rownames(.$AE)

          out$AE <- NULL


          as.matrix(out)[AE_sub, Drug_sub, drop = FALSE]
        }
      )

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


    # if show_text then display_number is not used
    if (show_text & !is.null(dots$display_numbers)) {
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

      msg <- glue::glue(
        "'display_numbers' automatically set to display \\
           the {to_show} when {which_show}"
      )

      warning(msg)

      dots$display_numbers <- NULL
    }


    if (show_text) {
      dots$display_numbers <- plot_mats$text
    }


    if (is.null(dots$color)) {
      dots$color <- RColorBrewer::brewer.pal(
        n = 9, name = "Blues"
      ) %>%
        rev() %>%
        {grDevices::colorRampPalette(.)(100)}
    }

    if (is.null(dots$number_color)) {
      dots$number_color <- plot_mats$text_color
    }


    if (is.null(dots$show_rownames)) {
      dots$show_rownames <- length(AE_sub) <= 50
    }

    if (is.null(dots$main)) {
      dots$main <- fill_measure
    }

    pheatmap_args <- list(
      mat = plot_mats[[fill_measure]],
      cluster_cols = FALSE,
      cluster_rows = FALSE
    ) %>%
      c(dots)

    out <- do.call(pheatmap::pheatmap, pheatmap_args)
  }

  out
}


#' @rdname heatmap_pvlrt
#' @export
barplot.pvlrt <- function(object,
                          AE = NULL,
                          Drug = NULL,
                          grep = FALSE,
                          x_axis_measure = "lrstat",
                          fill_measure = "p.value",
                          show_n = FALSE,
                          arrange_alphabetical = FALSE,
                          show_pvalue = FALSE,
                          show_lrstat = FALSE,
                          p.value_lower = 0,
                          p.value_upper = 1,
                          lrstat_lower = 0,
                          lrstat_upper = Inf,
                          n_lower = 0,
                          n_upper = Inf,
                          remove_outside = FALSE,
                          digits = 2,
                          Drug_nrow = 1,
                          x_axis_trans = FALSE,
                          x_axis_trans_fn = "log1p",
                          border_color = "white",
                          ...) {

  dots <- list(...)
  all_inputs <- c(as.list(environment()), dots)

  # processed plotting data from process_plot_data
  dat_pl <- tryCatch(
    do.call(process_plot_data, all_inputs),
    error = function(e) e
  )

  if (is(dat_pl, "error")) {
    stop(dat_pl$message)
  }

  show_text <- show_pvalue | show_n | show_lrstat

  darkblue_col <- RColorBrewer::brewer.pal(8, "Blues") %>% tail(1)


  dat_pl[
    ,
    AE := AE %>%
      factor(levels = rev(levels(.)))
  ]

  fill_range <- RColorBrewer::brewer.pal(n = 10, name = "RdBu") %>%
    rev()

  n_uniq_fill_meas <- dat_pl[[fill_measure]] %>%
    unique() %>%
    length()

  if (n_uniq_fill_meas == 1) {
    fill_range <- fill_range[1]
  }

  out <- dat_pl %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        y = "AE",
        x = x_axis_measure,
        fill = fill_measure,
        label = "text"
      )
    ) +
    ggplot2::geom_bar(stat = "identity", color = border_color) +
    ggplot2::facet_wrap(~Drug, nrow = Drug_nrow) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradientn(colors = fill_range)+
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = "")


  if (show_text) {
    out <- out +
      ggplot2::geom_text(
        hjust = "inward",
        position = ggplot2::position_dodge(width = 1),
        inherit.aes = TRUE,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(dat_pl$text_color) %>%
          unique() %>%
          rev(),
        guide = "none"
      )
  }

  if (x_axis_trans) {
    out <- out +
      ggplot2::coord_trans(x = x_axis_trans_fn)
  }

  out
}
