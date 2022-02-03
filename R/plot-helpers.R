process_plot_data <- function(object = object,
                              AE = NULL,
                              Drug = NULL,
                              grep = FALSE,
                              fill_measure = "p_value",
                              show_n = FALSE,
                              show_p_value = FALSE,
                              show_lrstat = FALSE,
                              arrange_alphabetical = FALSE,
                              p_value_lower = p_value_lower,
                              p_value_upper = p_value_upper,
                              lrstat_lower = lrstat_lower,
                              lrstat_upper = lrstat_upper,
                              n_lower = 0,
                              n_upper = Inf,
                              remove_outside = FALSE,
                              digits = 2,
                              fill_gradient_range = c("red", "white"),
                              ...) {
  stopifnot(
    is(object, "pvlrt"),
    fill_measure %in% c("p_value", "lrstat", "n"),
    all(is.character(fill_gradient_range))
  )


  . <- NULL
  n_text <- p_value_text <- lrstat_text <- text_color <- NULL

  dots <- list(...)
  all_inputs <- c(as.list(environment()), dots)

  old_args <- c("show_pvalue", "show_p.value")
  for (this_arg in old_args) {
    if (!is.null(dots[[this_arg]])) {
      msg <- glue::glue(
        "The argument '{this_arg}' is deprecated. Use 'show_p_value' instead."
      )
      show_p_value <- all_inputs$show_p_value <- dots[[this_arg]]
    }
  }

  processed <- list()

  # defined inside data.table using NSE
  lrstat <- p_value <- NULL

  # get summary statistics, sort by decreasing lr statistic
  summ <- summary(
    object
  )[
    ,
    `:=`(
      AE = as.character(AE),
      Drug = as.character(Drug)
    )
  ][
    !is.na(lrstat)
  ] %>%
    data.table::setorder(-lrstat)

  all_names <- list(
    AE = unique(summ$AE),
    Drug = unique(summ$Drug)
  )


  summ_AE_Drug <- summ # data.table::copy(summ)


  # first subset on AEs and Drugs #
  # keeps the above decreasing LRstat orders
  # if the provided AEs and Drugs are numeric or null
  # but uses the provided order if the provided
  # AEs and Drugs are characters
  measure_order <- if (all(is.numeric(AE))) {
    c("Drug", "AE")
  } else if (all(is.numeric(Drug))) {
    c("AE", "Drug")
  } else {
    c("Drug", "AE")
  }


  for (nm in measure_order) {
    tmp <- all_inputs[[nm]]
    if (is.null(tmp)) {
      n_elem <- all_names[[nm]] %>%
        length(.) %>%
        min(50)
      processed[[nm]] <- all_names[[nm]][1:n_elem]
    } else if (all(is.numeric(tmp))) {
      if (length(tmp) > 1) {
        msg <- glue::glue(
          "{nm}, if numeric, must have length 1"
        )
        stop(msg)
      }
      n_elem <- all_names[[nm]] %>%
        length(.) %>%
        min(tmp)
      processed[[nm]] <- summ_AE_Drug[[nm]] %>%
        unique() %>%
        .[1:n_elem]
    } else if (all(is.character(tmp))) {
      if (!grep) {
        extra <- setdiff(tmp, all_names[[nm]])
        if (length(extra) > 1) {
          msg <- glue::glue(
            "some of the requested {nm}s are \\
            not present in the pvlrt object.
            To match by pattern set grep = TRUE.
            Use `extract_{nm}_names()` to extract \\
            the available {nm}-names in the \\
            pvlrt object."
          )
          warning(msg)
        }
        matches <- tmp %>%
          intersect(all_names[[nm]]) %>%
          unique(.)
        if (length(matches) < 1) {
          msg <- glue::glue(
            "No {nm}s found in pvlrt object \\
            with the requested names.
            To match by pattern set grep = TRUE.
            Use `extract_{nm}_names()` to extract \\
            the available {nm}-names in the \\
            pvlrt object.
            "
          )
          stop(msg)
        }
      } else {
        pattern <- tmp %>%
          paste(collapse = "|") %>%
          gsub("\\|$", "", .) %>%
          trimws()
        # if (length(pattern) > 1) browser()
        matches <- grep(
          pattern = pattern,
          x = all_names[[nm]],
          ignore.case = TRUE,
          value = TRUE
        )
        if (length(matches) < 1) {
          msg <- glue::glue(
            "no {nm}s found in pvlrt object \\
            with the requested pattern.
            Use `extract_{nm}_names()` to extract \\
            the available {nm}-names in the \\
            pvlrt object.
            "
          )
          warning(msg)
        }
      }

      processed[[nm]] <- matches
    } else {
      msg <- glue::glue("invalid {nm} provided")
      stop(msg)
    }



    # arrange AE & Drug categories alphabetically,
    # if requested
    if (arrange_alphabetical) {
      processed[[nm]] <- sort(processed[[nm]])
    }

    # subset on the current set of nm,
    # rearrange the rows based on lrstat
    summ_AE_Drug <- glue::glue(
      "summ_AE_Drug[
        {nm} %in% unique(processed${nm})
      ] %>%
      data.table::setorder(-lrstat)"
    ) %>%
      parse(text = .) %>%
      eval(.)
  }


  # summ_AE_Drug <- summ[
  #   AE %in% unique(processed$AE) &
  #     Drug %in% unique(processed$Drug)
  # ]

  if (nrow(summ_AE_Drug) < 1) {
    msg <- glue::glue(
      "no result at the provided AE and Drug combination"
    )
    stop(msg)
  }

  # next threshold by measure -- extract all AE's and
  # Drugs that satisfy the thresholds in at least one
  # comparison
  meas_all <- c("p_value", "lrstat", "n")

  filter_char <-
    # specify the ranges for p_value, lrstat and n
    # as character
    glue::glue(
      "{meas_all} %between% c({meas_all}_lower, {meas_all}_upper)"
    ) %>%
    paste(collapse = " & ")

  summ_pval_lrt <- filter_char %>%
    # data table filtering statement as character
    {
      glue::glue("summ_AE_Drug[{.}]")
    } %>%
    # parse and eval text
    parse(text = .) %>%
    eval()


  if (nrow(summ_pval_lrt) < 1) {
    msg <- glue::glue(
      "no result at the provided p_value_lower, p_value_upper, \\
      lrstat_lower, and lrstat_upper combination"
    )
    stop(msg)
  }

  # subset on both (AE, Drug) subset
  # and measure value based subset
  dat_pl <- summ_AE_Drug[
    AE %in% unique(summ_pval_lrt$AE) &
      Drug %in% unique(summ_pval_lrt$Drug)
  ]


  # defined inside data.table using NSE
  p_value <- n <- .N <- .SD <-
    threshold <- display_text_color <- NULL

  dat_pl <- dat_pl[
    ,
    `:=`(
      AE = AE %>%
        factor(
          levels = intersect(
            processed$AE, unique(.)
          )
        ),
      Drug = Drug %>%
        factor(
          levels = intersect(
            processed$Drug, unique(.)
          )
        )
    )
  ][
    ,
    `:=`(
      p_value_text = p_value %>%
        format_pval_(digits = digits) %>%
        paste0("p", .),
      lrstat_text = lrstat %>%
        {
          ifelse(
            . >= 1e4,
            format(., digits = digits, scientific = TRUE),
            as.character(round(., digits = digits))
          )
        } %>%
        paste0("logLR=", .),
      n_text = n %>%
        {
          ifelse(
            . >= 1e4,
            format(., digits = digits, scientific = TRUE),
            as.character(round(., digits = digits))
          )
        } %>%
        paste0("n=", .)
    )
  ][
    ,
    `:=`(
      info = paste(
        n_text,
        lrstat_text,
        p_value_text,
        sep = "; "
      ),
      text = rep("", .N) %>%
        {
          if (show_n) paste(., n_text) else .
        } %>%
        {
          if (show_lrstat) paste(., lrstat_text, sep = "; ") else .
        } %>%
        {
          if (show_p_value) paste(., p_value_text, sep = "; ") else .
        } %>%
        gsub("^\\,", "", .) %>%
        # gsub("^\\ ", "", .) %>%
        # trim white space
        trimws()
    )
  ][
    ,
    threshold := .SD[[fill_measure]] %>%
      unlist() %>%
      c() %>%
      # median() %>%
      # {./2}
      range(.) %>%
      mean(.)
  ][
    ,
    text_color := .SD[[fill_measure]] %>%
      {
        ifelse(
          . >= threshold,
          "black",
          "orange"
        )
      }
  ][,
    c(
      "AE", "Drug", "n", "lrstat",
      "p_value", "info", "text", "text_color"
    ),
    with = FALSE
  ]


  if (remove_outside) {
    assign_na_char <- glue::glue(
      "{meas_all} = NA"
    ) %>%
      paste(collapse = ", ") %>%
      {
        glue::glue("`:=`({.})")
      }

    glue::glue(
      "dat_pl[
     !({filter_char}),
     {assign_na_char}
     ]"
    ) %>%
      parse(text = .) %>%
      eval()
  }


  dat_pl
}
