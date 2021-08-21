process_plot_data <- function(object = object,
                              AE = NULL,
                              Drug = NULL,
                              grep = FALSE,
                              measure = "p.value",
                              show_n = FALSE,
                              show_pvalue = FALSE,
                              show_lrstat = FALSE,
                              arrange_alphabetical = FALSE,
                              p.value_lower = p.value_lower,
                              p.value_upper = p.value_upper,
                              lrstat_lower = lrstat_lower,
                              lrstat_upper = lrstat_upper,
                              n_lower = 0,
                              n_upper = Inf,
                              remove_outside = FALSE,
                              digits = 2,
                              ...) {


  stopifnot(
    is(object, "pvlrt"),
    measure %in% c("p.value", "lrstat", "n")
  )

  dots <- list(...)
  all_inputs <- c(as.list(environment()), dots)
  processed <- list()

  # defined inside data.table using NSE
  lrstat <- p.value <- AE <- Drug <-  NULL

  # First threshold by the provided limists of p.values and the LRstats
  summ_full <- summary(object) %>%
    .[,
      `:=`(
        AE = as.character(AE),
        Drug = as.character(Drug)
      )
    ]



  meas_all <- c("p.value", "lrstat", "n")

  filter_char <-
    # specify the ranges for p.value, lrstat and n
    # as character
    glue::glue(
      "{meas_all} %between% c({meas_all}_lower, {meas_all}_upper)"
    ) %>%
    paste(collapse = " & ")

  summ_sub_pval_lrt <- filter_char %>%
    # data table filtering statement as character
    {glue::glue("summ_full[{.}]")} %>%
    # parse and eval text
    parse(text = .) %>%
    eval()


  if (nrow(summ_sub_pval_lrt) < 1) {
    msg <- glue::glue(
      "no result at the provided p.value_lower, p.value_upper, \\
      lrstat_lower, and lrstat_upper combination"
    )
    stop(msg)
  }


  summ_sub <- summ_full[
    AE %in% unique(summ_sub_pval_lrt$AE) &
      Drug %in% unique(summ_sub_pval_lrt$Drug)
  ] %>%
    data.table::setorder(-lrstat)

  all_names <- list(
    AE = unique(summ_sub$AE),
    Drug = unique(summ_sub$Drug)
  )


  # keeps the above decreasing LRstat orders
  # if the provided AEs and Drugs are numeric or null
  # but uses the provided order if the provided
  # AEs and Drugs are characters
  for (nm in c("AE", "Drug")) {
    tmp <- all_inputs[[nm]]
    if (is.null(tmp)) {
      n_elem <- all_names[[nm]] %>% length(.) %>% min(50)
      processed[[nm]] <- all_names[[nm]][1:n_elem]
    } else if (all(is.numeric(tmp))) {
      if (length(tmp) > 1) {
        msg <- glue::glue(
          "{nm}, if numeric, must have length 1"
        )
        stop(msg)
      }
      n_elem <- all_names[[nm]] %>% length(.) %>% min(tmp)
      processed[[nm]] <- all_names[[nm]][1:n_elem]
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
        matches <- tmp %>% intersect(all_names[[nm]]) %>% unique(.)
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

  }

  # arrange AE & Drug categories alphabeticallys
  if (arrange_alphabetical) {
    processed <- lapply(processed, sort)
  }


  dat_pl <- summ_sub[
    AE %in% unique(processed$AE) &
      Drug %in% unique(processed$Drug)
  ]

  if (nrow(dat_pl) < 1) {
    msg <- glue::glue(
      "no result at the provided AE, Drug, p.value_lower, p.value_upper, \\
      lrstat_lower, and lrstat_upper combination"
    )
    stop(msg)
  }

  # defined inside data.table using NSE
  p.value <- n <- .N <- .SD <-
    threshold <- display_text_color <-  NULL

  dat_pl <- dat_pl[
    ,
    `:=`(
      AE = AE %>% factor(levels = processed$AE),
      Drug = Drug %>% factor(levels = processed$Drug)
    )
  ][
    ,
    `:=`(
      p.value_text = p.value %>%
        format_pval_(digits = digits) %>%
        paste0("p", .),

      lrstat_text = lrstat %>%
        {ifelse(
          . >= 1e4,
          format(., digits = digits, scientific = TRUE),
          as.character(round(., digits = digits))
        )} %>%
        paste0("LR=", .),

      n_text = n %>%
        {ifelse(
          . >= 1e4,
          format(., digits = digits, scientific = TRUE),
          as.character(round(., digits = digits))
        )} %>%
        paste0("n=", .)
    )
  ][,
    `:=`(
      text = rep("", .N) %>%
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
        trimws()
    )
  ][,
    threshold := .SD[[measure]] %>%
      unlist() %>%
      c() %>%
      # median() %>%
      # {./2}
      range(.) %>%
      mean(.)
  ][,
    text_color := .SD[[measure]] %>%
      {ifelse(
        . >= threshold,
        "black",
        "orange"
      )}
  ][,
    c("AE", "Drug", "n", "lrstat",
      "p.value", "text", "text_color"),
    with = FALSE
  ]


  if (remove_outside) {
    assign_na_char <- glue::glue(
      "{meas_all} = NA"
    ) %>%
      paste(collapse = ", ") %>%
      {glue::glue("`:=`({.})")}

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
