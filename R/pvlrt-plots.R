#' Heatmap, barplot and bubbleplot displaying liklihood ratio test results in
#' a `pvlrt` object
#' @param object,height pvlrt object; output of \code{pvlrt}()
#' @param ... Other arguments. Currently ignored
#' @param fill_measure Measure to govern the filling color in each cell
#' (in heatmap) or bar (in barplot) or circle/bubble (in bubbleplot) for each
#' drug/AE combination. Defaults to "p_value". Available choices are:
#'  "p.value", "lrstat", and "n".
#' @param AE input parameter determining which adverse effects to show in
#' the plot. This can be a numeric scalar  specifying the
#' number of *top* (in terms of computed LRT values) adverse effects to show.
#' Alternatively, it can be a character vector, specifying the exact
#' adverse effects to show. It can also be a vector of patterns to match
#' (ignores cases) against the full names of all available adverse effects,
#' provided \code{grep} is set to TRUE. Defaults to adverse effects
#' corresponding to the top M pairs where M = max(number of possible pairs, 50).
#' Set AE = Inf to force display of all adverse effects.
#' @param Drug input parameter determining which drugs to show in
#' the plot. This can be a numeric scalar  specifying the
#' number of *top* (in terms of computed LRT values) drugs to show.
#' Alternatively, it can be a character vector, specifying the exact
#' drugs to show. It can also be a vector of patterns to match
#' (ignores cases) against the full names of all available drugs,
#' provided \code{grep} is set to TRUE. Defaults to drugs
#' corresponding to the top M pairs where M = max(number of possible pairs, 50).
#' Set Drug = Inf to force display all drugs.
#' @param grep logical. Match patterns against the supplied AE or Drug names?
#' Ignores if neither AE nor Drug is a character vector.
#' @param show_n logical. show the sample size as inscribed text on each cell?
#' @param show_p_value logical. show the computed p-value as inscribed text on each cell?
#' @param show_lrstat logical. show the computed LRT statistic (on log-scale)
#' inscribed text on each cell?
#' @param arrange_alphabetical logical. should the rows (AEs) and columns (Drugs)
#'  be arranged in alphabetical orders? Defaults to FALSE, in which case
#'  the orderings of the computed LRT values are used.
#' @param p_value_lower,p_value_upper lower and upper limits on the computed p-values to
#' display on the plot.
#' @param lrstat_lower,lrstat_upper lower and upper limits on the computed LRT values to
#' display on the plot.
#' @param n_lower,n_upper lower and upper limits on the computed sample sizes to
#' display on the plot.
#' @param remove_outside logical. Should the values for pairs with p-value, LRT
#' statistics or sample sizes falling outside of the provided ranges through p_value_lower, p_value_upper
#' etc., be replaced with \code{NA}? Defaults to FALSE. Setting this to TRUE may help
#' distinguish drugs or AEs which has some pairs falling within and some pairs falling
#' outside of the provided ranges better.
#' @param digits numeric. Number of decimal places to show on the
#' inscribed texts on the plot.
#' @param border_color character string. Specifies the border color of cells/bars.
#' @param Drug_nrow Number of rows in the panels for Drugs for the barplots.
#' @param x_axis_measure measure to show on the x-axis of the (horizontal) bar
#' plots and bubble plots. Defaults to "lrstat". Available choices are "lrstat", "p_value" and "n".
#' @param size_measure measure to govern sizes of the circles in the bubble plot.
#' Defaults to "n". Available choices are "lrstat", "p_value" and "n".
#' @param x_axis_logscale logical. Should the x axis measure in the bar plot or the bubble plot
#' be log transformed (more precisely, "log1p" transformed with the function
#' f(x) = log(1 + x))? Defaults to TRUE.
#' @param size_logscale logical. Should the circle size measure in the the bubble plot
#' be log transformed (more precisely, "log1p" transformed with the function
#' f(x) = log(1 + x)). Defaults to TRUE.
#' @param fill_gradient_range character vector. Specifies the range of gradient colors used
#' for `fill_measure`. Passed into the \code{colours} argument of
#' `scale_fill_gradientn` from \code{ggplot2}.
#'
#'
#' @examples
#' # 500 bootstrap iterations (nsim) in the example below
#' # are for quick demonstration only --
#' # we recommended setting nsim to 10000 (default) or bigger
#' test1 <- pvlrt(statin46, nsim = 500)
#' bubbleplot_pvlrt(test1)
#' heatmap_pvlrt(test1)
#' barplot(test1)
#'
#'
#' @seealso
#' \link{pvlrt}
#'
#' @return
#' A \link[ggplot2]{ggplot} object.
#'
#' @md
#' @export
heatmap_pvlrt <- function(object,
                          AE = NULL,
                          Drug = NULL,
                          grep = FALSE,
                          fill_measure = "p_value",
                          show_n = FALSE,
                          show_lrstat = FALSE,
                          show_p_value = FALSE,
                          p_value_lower = 0,
                          p_value_upper = 1,
                          lrstat_lower = 0,
                          lrstat_upper = Inf,
                          n_lower = 0,
                          n_upper = Inf,
                          arrange_alphabetical = FALSE,
                          remove_outside = FALSE,
                          digits = 2,
                          border_color = "black",
                          fill_gradient_range = c("darkred", "white"),
                          ...) {
  . <- NULL
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


  show_text <- show_p_value | show_n | show_lrstat

  # heatmap engine specific codes
  # pheatmap support is turned off

  # stopifnot(
  #   engine %in% c("pheatmap", "ggplot", "ggplot2")
  # )


  # if (grepl("ggplot", engine)) {
  dat_pl[
    ,
    AE := AE %>%
      factor(levels = rev(levels(.)))
  ]

  darkblue_col <- RColorBrewer::brewer.pal(8, "Blues") %>% tail(1)

  if (show_text) {
    dat_pl$info <- dat_pl$text
  }

  out <- dat_pl %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        y = "AE",
        x = "Drug",
        fill = fill_measure,
        label = "info"
      )
    ) +
    ggplot2::geom_tile(color = border_color) +
    # ggplot2::scale_fill_gradient(
    #   low = ifelse(fill_measure == "p_value", darkblue_col, "white"),
    #   high = ifelse(fill_measure == "p_value", "white", darkblue_col)
    # ) +
    ggplot2::scale_fill_gradientn(colors = fill_gradient_range) +
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
      ggfittext::geom_fit_text(
        reflow = TRUE,
        contrast = TRUE,
        grow = TRUE
      )
  }

  # } else if (engine == "pheatmap") {
  #
  #   if (!requireNamespace("pheatmap")) {
  #     stop("package {pheatmap} is required for engine='pheatmap'")
  #   }
  #
  #   AE_sub <- levels(dat_pl$AE)
  #   Drug_sub <- levels(dat_pl$Drug)
  #
  #   plot_mats <- c(
  #     "lrstat", "p_value",
  #     "text", "text_color"
  #   ) %>%
  #     setNames(., .) %>%
  #     lapply(
  #       function(fill_measure) {
  #         out <- dat_pl %>%
  #           data.table::dcast(
  #             AE ~ Drug,
  #             value.var = fill_measure
  #           ) %>%
  #           as.data.frame() %>%
  #           set_rownames(.$AE)
  #
  #         out$AE <- NULL
  #
  #
  #         as.matrix(out)[AE_sub, Drug_sub, drop = FALSE]
  #       }
  #     )
  #
  #   unused_args <- c("cluster_rows", "cluster_cols")
  #   for (nm in unused_args) {
  #     if (!is.null(dots[[nm]])) {
  #       msg <- glue::glue(
  #         "{nm} is set to FALSE when passed to `pheatmap()`"
  #       )
  #       warning(msg)
  #       dots[[nm]] <- NULL
  #     }
  #   }
  #
  #
  #   # if show_text then display_number is not used
  #   if (show_text & !is.null(dots$display_numbers)) {
  #     which_show <- "" %>%
  #       {
  #         if (show_n) paste(., "show_n = TRUE", sep = ", ") else .
  #       } %>%
  #       {
  #         if (show_lrstat) paste(., "show_lrstat = TRUE", sep = ", ") else .
  #       } %>%
  #       {
  #         if (show_p_value) paste(., "show_p_value = TRUE", sep = ", ") else .
  #       } %>%
  #       gsub("^\\,", "", .) %>%
  #       # gsub("^\\ ", "", .) %>%
  #       # trim white space
  #       trimws()
  #
  #     to_show <- "" %>%
  #       {
  #         if (show_n) paste(., "n's", sep = ", ") else .
  #       } %>%
  #       {
  #         if (show_p_value) paste(., "p-value's", sep = ", ") else .
  #       } %>%
  #       {
  #         if (show_lrstat) paste(., "lrstat's", sep = ", ") else .
  #       } %>%
  #       gsub("^\\,", "", .) %>%
  #       # gsub("^\\ ", "", .) %>%
  #       # trim white space
  #       trimws()
  #
  #     msg <- glue::glue(
  #       "'display_numbers' automatically set to display \\
  #          the {to_show} when {which_show}"
  #     )
  #
  #     warning(msg)
  #
  #     dots$display_numbers <- NULL
  #   }
  #
  #
  #   if (show_text) {
  #     dots$display_numbers <- plot_mats$text
  #   }
  #
  #
  #   if (is.null(dots$color)) {
  #     dots$color <- RColorBrewer::brewer.pal(
  #       n = 9, name = "Blues"
  #     ) %>%
  #       rev() %>%
  #       {grDevices::colorRampPalette(.)(100)}
  #   }
  #
  #   if (is.null(dots$number_color)) {
  #     dots$number_color <- plot_mats$text_color
  #   }
  #
  #
  #   if (is.null(dots$show_rownames)) {
  #     dots$show_rownames <- length(AE_sub) <= 50
  #   }
  #
  #   if (is.null(dots$main)) {
  #     dots$main <- fill_measure
  #   }
  #
  #   pheatmap_args <- list(
  #     mat = plot_mats[[fill_measure]],
  #     cluster_cols = FALSE,
  #     cluster_rows = FALSE
  #   ) %>%
  #     c(dots)
  #
  #   out <- if (requireNamespace("pheatmap")) {
  #     do.call(pheatmap::pheatmap, pheatmap_args)
  #   } else NULL
  # }

  out
}


#' @rdname heatmap_pvlrt
#' @export
barplot.pvlrt <- function(height,
                          AE = NULL,
                          Drug = NULL,
                          grep = FALSE,
                          x_axis_measure = "lrstat",
                          fill_measure = "p_value",
                          show_n = FALSE,
                          arrange_alphabetical = FALSE,
                          show_p_value = FALSE,
                          show_lrstat = FALSE,
                          p_value_lower = 0,
                          p_value_upper = 1,
                          lrstat_lower = 0,
                          lrstat_upper = Inf,
                          n_lower = 0,
                          n_upper = Inf,
                          remove_outside = FALSE,
                          digits = 2,
                          Drug_nrow = 1,
                          border_color = "black",
                          x_axis_logscale = TRUE,
                          fill_gradient_range = c("darkred", "white"),
                          ...) {
  . <- NULL
  object <- height
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

  show_text <- show_p_value | show_n | show_lrstat

  darkblue_col <- RColorBrewer::brewer.pal(8, "Blues") %>% tail(1)


  dat_pl[
    ,
    AE := AE %>%
      factor(levels = rev(levels(.)))
  ]

  # fill_range <- RColorBrewer::brewer.pal(n = 10, name = "RdBu") %>%
  #   rev()

  n_uniq_fill_meas <- dat_pl[[fill_measure]] %>%
    unique() %>%
    length()

  if (n_uniq_fill_meas == 1) {
    fill_range <- fill_range[1]
  }

  if (show_text) {
    dat_pl$info <- dat_pl$text
  }

  if (x_axis_logscale) {
    x_axis_measure_new <- paste0("log1p_", x_axis_measure)
    # x_axis_label <- bquote('log'[10](1 + .(x_axis_measure)))
    x_axis_label <- glue::glue("log10(1 + {x_axis_measure})") %>% as.character()
    dat_pl[[x_axis_measure_new]] <- log(1 + dat_pl[[x_axis_measure]], base = 10)
  } else {
    x_axis_measure_new <- x_axis_label <- x_axis_measure
    dat_pl[[x_axis_measure_new]] <- dat_pl[[x_axis_measure]]
  }


  out <- dat_pl %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        y = "AE",
        x = x_axis_measure_new,
        fill = fill_measure,
        label = "info"
      )
    ) +
    ggplot2::geom_bar(stat = "identity", color = border_color) +
    ggplot2::facet_wrap(~Drug, nrow = Drug_nrow) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradientn(colors = fill_gradient_range) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = "", x = x_axis_label)


  if (show_text) {
    out <- out +
      ggfittext::geom_bar_text(contrast = TRUE)

  }

  out
}





#' @rdname heatmap_pvlrt
#' @param Drug_nrow Number of rows in the panels for Drugs for the barplots.
#' @param x_axis_measure measure to show on the x-axis of the (horizontal) bar
#' plots. Defaults to "lrstat" available choices are "lrstat", "p_value" and "n".
#' @export
bubbleplot_pvlrt <- function(object,
                             AE = NULL,
                             Drug = NULL,
                             grep = FALSE,
                             x_axis_measure = "lrstat",
                             fill_measure = "p_value",
                             size_measure = "n",
                             show_n = FALSE,
                             arrange_alphabetical = FALSE,
                             show_p_value = FALSE,
                             show_lrstat = FALSE,
                             p_value_lower = 0,
                             p_value_upper = 1,
                             lrstat_lower = 0,
                             lrstat_upper = Inf,
                             n_lower = 0,
                             n_upper = Inf,
                             remove_outside = FALSE,
                             digits = 2,
                             Drug_nrow = 1,
                             border_color = "black",
                             x_axis_logscale = TRUE,
                             size_logscale = TRUE,
                             fill_gradient_range = c("darkred", "white"),
                             ...) {
  . <- NULL
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

  show_text <- show_p_value | show_n | show_lrstat

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


  if (show_text) {
    dat_pl$info <- dat_pl$text
  }


  out <- dat_pl %>%
    ggplot2::ggplot(
      ggplot2::aes_string(
        y = "AE",
        x = x_axis_measure,
        fill = fill_measure,
        size = size_measure,
        label = "info"
      )
    ) +
    ggplot2::geom_point(stat = "identity", color = border_color, shape = 21) +
    ggplot2::facet_wrap(~Drug, nrow = Drug_nrow) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradientn(colors = fill_gradient_range) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1
      ),
      # panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = "")

  if (x_axis_logscale) {
    out <- out +
      ggplot2::scale_x_continuous(trans = "log1p")
  }

  if (size_logscale) {
    out <- out +
      ggplot2::scale_size_continuous(trans = "log1p")
  }


  if (show_text) {
    out <- out + ggplot2::geom_text()
  }

  out
}


