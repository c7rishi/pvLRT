#' Convert raw AE-Drug incidence data into a contingency table
#'
#'
#' @param rawdata a `data.frame` or an object that can be converted to a `data.frame`.
#' Must contain 3 columns (i) DRUG: the drug names/labels, (ii) AE: the AE names, and
#' either (iii) CASEID: case ids corresponding to each combination of AE and DRUG,
#' (if `aggregated` is `FALSE`), or (iii') COUNT: the total number of incidences of
#' each AE and DRUG combination (if `aggregated` is `TRUE`). If these columns are named
#' differently in `rawdata`, supply the appropriate column names through `Drug_col_name`,
#' `AE_col_name`, `id_col_name` and `count_col_name`.
#'
#' @param Drug_col_name,AE_col_name,id_col_name,count_col_name Drug, AE, case id and count
#' column names in `rawdata`. Defaults to `DRUG`, `AE`, `CASEID` and `COUNT`.
#'
#' @param aggregated logical. Has the incidences been already aggregated/summarized into counts
#' in `rawdata`? If `TRUE` then then `COUNT` column is used to produce the output
#' contingency table. If `FALSE` (default) incidences are first aggregated into counts
#' before converting to contingency tables.
#'
#' @param create_other_AE_row logical. Add a row in the contingency table for "Other AEs"? This can aid computation
#' in situations where there are certain AEs of primary interest. See `other_AE_excludes` for details on
#' how to specify the "Other AE" row.
#' @param other_AE_excludes character vector cataloging AEs that are NOT to be included in the
#' row for Other AEs. If NULL (default) then then no AEs are included in Other AEs (i.e.,
#' `other_AE_excludes` contains all AEs in the raw data). Ignored if  `create_other_AE_row = FALSE`.
#' @param other_AE_rowname character. Row name for the "Other AE" row created. Defaults to "Other AE". Ignored if
#' `create_other_AE_row = FALSE`.
#'
#' @param create_other_Drug_col logical. Add a column in the contingency table for "Other Drugs"? This
#' column plays the role of a "baseline" group of drugs that typically do not indicate an adverse event
#' association with the signal of interest. Care should be taken while determining which drugs to
#' include in this group; See Ding et al (2020) for guidance.
#' @param other_Drug_excludes character vector cataloging Drugs that are NOT to be included in the
#' column for Other Drugs. If NULL (default) then then no Drugs are included in Other Drugs (i.e.,
#' `other_Drug_excludes` contains all Drugs in the raw data). Ignored if
#' `create_other_Drug_col = FALSE`.
#' @param other_Drug_colname character. Row name for the "Other Drug" column created. Ignored if
#' `create_other_Drug_col = FALSE`.
#'
#' @param ... unused.
#'
#' @details
#'
#' This is a convenience function that creates a contingency table cataloging counts of
#' AE-Drug incidences from a raw Drug/AE incidence data frame.
#' It accepts both raw incidence data (each row is one incidence of a Drug-AE combination,
#' indexed by case ids) and summarized count data (each row catalogs the total counts of incidences
#' of each Drug-AE pair). The output is a matrix (contingency table) enumerating total count of cases for
#' each pair of AE (along the rows) and drug (along the columns) with appropriately
#' specified row and column names, and can be passed to a pvlrt() call. See the examples for more details.
#'
#'
#' The output can be fed into \link{pvlrt} or its wrappers as `contin_table`
#'
#' @references
#' Ding, Y., Markatou, M. and Ball, R., 2020. An evaluation of statistical approaches to postmarketing surveillance. Statistics in Medicine, 39(7), pp.845-874.
#'
#' Chakraborty, S., Liu, A., Ball, R. and Markatou, M., 2022. On the use of the likelihood ratio test methodology in pharmacovigilance. Statistics in Medicine, 41(27), pp.5395-5420.
#'
#' @examples
#'
#' # convert to contingency table form incidence (non-aggregated) raw data
#' # AE subset = AEs in statin46
#' # Durg subset = union of statin46 and gbca drugs
#' tab1 <- convert_raw_to_contin_table(
#'   rawdata = faers22q3raw,
#'   Drug_col_name = "DRUG",
#'   AE_col_name = "AE",
#'   id_col_name = "CASEID",
#'   aggregated = FALSE,
#'   other_AE_excludes = rownames(statin46),
#'   other_Drug_excludes = union(colnames(gbca), colnames(statin)),
#'   create_other_Drug_col = TRUE,
#'   create_other_AE_row = FALSE
#' )
#'
#' # convert to contingency table AFTER aggregating and counting
#' # the total number of incidences of each (AE, Drug) pair
#' ## Same AE and Drug subsets as before
#' ## aggregation (counting) done using data.table dt[i, j, by] syntax
#' ## uses magrittr %>% pipe
#' tab2 <- data.table::as.data.table(
#'   faers22q3raw
#' )[
#'   ,
#'   .(COUNT = length(unique(CASEID))),
#'   by = .(DRUG, AE)
#' ] %>%
#'   convert_raw_to_contin_table(
#'     Drug_col_name = "DRUG",
#'     AE_col_name = "AE",
#'     count_col_name = "COUNT",
#'     aggregated = TRUE,
#'     other_AE_excludes = rownames(statin46),
#'     other_Drug_excludes = union(colnames(gbca), colnames(statin)),
#'     create_other_Drug_col = TRUE,
#'     create_other_AE_row = FALSE
#'   )
#'
#' all.equal(tab1, tab2)
#'
#' # use the contingency table produced above in pvlrt()
#' ## 500 bootstrap iterations (nsim) in the example below
#' ## is for quick demonstration only --
#' ## we recommended setting nsim to 10000 (default) or bigger
#' test1 <- pvlrt(tab1, nsim = 500)
#'
#' @md
#' @export
convert_raw_to_contin_table <- function(rawdata,
                                        Drug_col_name = "DRUG",
                                        AE_col_name = "AE",
                                        id_col_name = "CASEID",
                                        count_col_name = "COUNT",
                                        aggregated = FALSE,
                                        create_other_Drug_col = FALSE,
                                        other_Drug_excludes = NULL,
                                        other_Drug_colname = "Other_Drug",
                                        create_other_AE_row = FALSE,
                                        other_AE_excludes = NULL,
                                        other_AE_rowname = "Other_AE",
                                        ...) {

  tst_2_df <- tryCatch(
    {
      rawdata <- as.data.frame(rawdata)
    },
    error = function(e) e
  )

  if (is(tst_2_df, "error")) stop(tst_2_df$message)

  stopifnot(
    is(rawdata, "data.frame"),
    Drug_col_name %in% names(rawdata),
    AE_col_name %in% names(rawdata),
    is.logical(aggregated),
    is.logical(create_other_Drug_col),
    is.logical(create_other_AE_row)
  )

  if (aggregated) {
    if (! count_col_name %in% names(rawdata)) {
      msg <- glue::glue(
        "aggregated = TRUE but column '{count_col_name}' \\
         (count_col_name) is not a column of rawdata"
      )
      stop(msg)
    }
  } else {
    if (! id_col_name %in% names(rawdata)) {
      msg <- glue::glue(
        "aggregated = FALSE but column '{id_col_name}' \\
         (id_col_name) is not a column of rawdata"
      )
      stop(msg)
    }
  }


  # for safer handling with data.table variables in package
  . <- AE <- DRUG <-
    AE_mod <- DRUG_mod <-
    COUNT <- CASEID <-
    ncases <- NULL

  dat_in <- as.data.table(rawdata) %>%
    data.table::setnames(
      old = c(Drug_col_name, AE_col_name),
      new = c("DRUG", "AE")
    ) %>%
    {
      if (!aggregated) {
        data.table::setnames(
          .,
          old = c(id_col_name),
          new = c("CASEID")
        ) %>%
          .[, .(CASEID, AE, DRUG)]
      } else {
        data.table::setnames(
          .,
          old = c(count_col_name),
          new = c("COUNT")
        ) %>%
          .[, .(AE, DRUG, COUNT)]
      }
    }

  if (is.null(other_AE_excludes)) other_AE_excludes <- unique(dat_in$AE)
  if (is.null(other_Drug_excludes)) other_Drug_excludes <- unique(dat_in$DRUG)

  dat_in <- dat_in[
    ,
    `:=`(
      AE_mod = AE,
      DRUG_mod = DRUG
    )
  ] %>%
    {
      if (create_other_AE_row)
        .[
          !(AE %in% other_AE_excludes),
          AE := other_AE_rowname
        ]
      else .
    } %>%
    {
      if (create_other_Drug_col)
        .[
          !(DRUG %in% other_Drug_excludes),
          DRUG_mod := other_Drug_colname
        ]
      else .
    }


  levels_AE <- unique(dat_in$AE_mod) %>%
    sort(decreasing = TRUE) %>%
    {
      if (create_other_AE_row) c(other_AE_rowname, .) %>% unique()
      else .
    } %>%
    rev()


  levels_Drug <- unique(dat_in$DRUG_mod) %>%
    sort(decreasing = TRUE) %>%
    {
      if (create_other_Drug_col) c(other_Drug_colname, .) %>% unique()
      else .
    } %>%
    rev()


  dat_count <- if (!aggregated) {
    dat_in[
      ,
      .(ncases = length(unique(CASEID))),
      keyby = .(AE_mod, DRUG_mod)
    ]
  } else {
    dat_in[
      ,
      .(ncases = sum(COUNT)),
      keyby = .(AE_mod, DRUG_mod)
    ]
  }


  out_mat <- dat_count[
    ,
    `:=`(
      AE_mod = factor(AE_mod, levels = levels_AE),
      DRUG_mod = factor(DRUG_mod, levels = levels_Drug)
    )
  ] %>%
    data.table::dcast(
      AE_mod ~ DRUG_mod,
      value.var = "ncases",
      fill = 0
    ) %>%
    {
      rn <- .$AE_mod
      .[, AE_mod := NULL] %>%
        as.data.frame() %>%
        set_rownames(rn)
    } %>%
    data.matrix()

  out_mat
}


