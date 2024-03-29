# uses  |> pipe, so, requires R 4.1.0 or newer
library(tidyverse)
# library(vroom)
# library(future)
library(pvLRT)
library(data.table)
# ncore <- availableCores()
# plan(multicore) #uncomment for parallel
set.seed(42)

important_drug <- colnames(statin) %>%
  c(colnames(gbca)) %>%
  unique() %>%
  setdiff("Other")

important_pt <- rownames(statin) %>%
  c(rownames(gbca)) %>%
  unique()

# from the supplementary document,
# contains list of important drugs and AEs
# tab_s4 <- readxl::read_excel("articles/table_s4_old.xlsx")

# # important drugs
# important_drug <- names(tab_s4) %>%
#   setdiff(
#     c("PTID", "PT",
#       "Other Drugs",
#       "PT Marginal Total",
#       "Estimated Probability")
#   )


library(pvLRT)

# important AE
# important_pt <- tab_s4$PT %>% na.omit()


process_data_dir <- "data-raw/fda-processed-data"
if (!dir.exists(process_data_dir)) {
  dir.create(process_data_dir)
}

all_files <-
  # list all zip files in the data/ascii directory
  list.files(
    "data-raw/fda-faers-data/ascii/",
    pattern = "\\.zip",
    full.names = TRUE
  ) %>%
  str_subset("2022") %>%
  str_subset("(?i)Q3")

# separate out faers and aers files
files_faers <-  grep(
  "faers",
  all_files,
  value = TRUE,
  ignore.case = TRUE
)
files_aers <- setdiff(all_files, files_faers)

# extract and process drug and AE data from
# each data file
read_data_faers_files <- function(this_zip, ...) {

  # browser()

  # full path
  full_zip_path <- this_zip

  # make a list of all elements inside the zip file
  this_zip_contents <- glue::glue(
    'unzip -l {this_zip}'
  ) %>%
    system(intern = TRUE) %>%
    strsplit(" ") %>%
    map_chr(tail, 1)


  # tmp <- unz(filename = full_zip_path,
  #            glue::glue(
  #   "ascii/ther12q4.txt"
  # ))
  #
  # tmp <- unz(
  #   full_zip_path,
  #   glue::glue("ascii/outc12q4.txt")
  #   # "ascii/DRUG04Q1.TXT"
  # ) %>%
  #   vroom::vroom(delim = "$")
  #
  #
  # tmp

  # extract the quarter name from the zip path
  this_quarter_name <- this_zip %>%
    strsplit("_") %>%
    map_chr(tail, 1) %>%
    strsplit("\\.") %>%
    map_chr(1) %>%
    toupper()


  # remove the "decade" 20 from the beginning of the
  # string
  this_quarter_name_short <- this_quarter_name %>%
    str_sub(3)


  # ascii drug file

  dat_files <- list(
    "DRUG" = "DRUG",
    "REACT" = "REAC"
  ) %>%
    lapply(
      function(this_type) {

        # browser()
        # if (this_type == "REAC") browser()
        # DRUG or REAC file
        this_filename_string_txt <- glue::glue(
          "{this_type}{this_quarter_name_short}.txt"
        )

        # elements in the zip file that matches
        # this_filename_string_txt
        this_zipfile_element <- this_zip_contents %>%
          grep(
            this_filename_string_txt,
            .,
            value = TRUE,
            ignore.case = TRUE
          )

        cat("** Now reading", this_filename_string_txt, "\n\n")

        # read_cmd <- glue::glue(
        #   "unzip -p {full_zip_path} {this_zipfile_element} | tr \\$ ' '"
        # )

        # extract the specific file from the zip,
        # then paste on the linux command line
        # then replace multiple $ delimeters by a single $
        read_cmd <- glue::glue(
          "unzip -p <<full_zip_path>> <<this_zipfile_element>> |\\
          sed 's/\\$\\{2,\\}/\\$/g'",
          .open = "<<",
          .close = ">>"
        )


        dat <- tryCatch(
          data.table::fread(
            cmd = read_cmd,
            fill = TRUE,
            sep = "$"
          ),
          error = function(e) e
        )

        if (is(dat, "error")) browser()

        data.table::setnames(
          dat,
          names(dat),
          names(dat) %>% toupper()
        )



        # glue::glue(
        #   "unzip -l {full_zip_path}"
        # ) %>%
        #   system()
        #
        #
        # # create a connection with the data file
        # # inside the zip file, and then read using vroom
        # dat <- tryCatch(
        #   unz(
        #     full_zip_path,
        #     glue::glue("{this_zipfile_element}")
        #     # "ascii/DRUG04Q1.TXT"
        #   ) %>%
        #     read_lines() %>%
        #     str_replace_all("\\$", "\t") %>%
        #     cat(file = tmpfile)
        #
        #   pp <- data.table::fread(tmpfile)
        #
        #
        #   vroom::vroom(delim = "$"), #%>%
        #   # mutate(QUARTER = this_quarter_name),
        #   error = function(e) e
        # )

        if (is(dat, "error")) browser()

        target_varnames <- c(
          "PRIMARYID", "CASEID",
          "ROLE_COD", "DRUGNAME",
          "PROD_AI", "PT"
        )

        current_varnames <- sapply(
          target_varnames,
          function(xx) {
            grep(
              xx, names(dat),
              ignore.case = TRUE, value = TRUE
            )
          }
        ) %>%
          unlist()


        out <- dat %>%
          # select only the variables of interest
          select(
            all_of(unname(current_varnames))
          ) %>%
          # clean variable names
          data.table::setDT() %>%
          unique() %>%
          data.table::setnames(
            old = unname(current_varnames),
            new = names(current_varnames)
          ) %>%
          data.table::setkey(PRIMARYID, CASEID) %>%
          unique()
        # # separate filtering for DRUG & REAC data
        # {
        #   if (this_type == "DRUG") {
        #     # subset rows on  ROLE_COD %in% c('PS', 'SS')
        #     # .[
        #     #   grepl("PS|SS", ROLE_COD, ignore.case = TRUE)
        #     # ][
        #     #   ,
        #     #   ROLE_COD := NULL
        #     # ] %>%
        #     #   unique()
        #     unique(.)
        #
        #   } else {
        #
        #     browser()
        #     tmp <- .
        #     # filter out all AEs outside the 46 in the list
        #     tmp[
        #       str_detect(
        #         PT,
        #         paste(important_pt, collapse = "|") %>%
        #           # ignore cases
        #           paste0("(?i)", .),
        #       )
        #       # grepl(
        #       #   paste(important_pt, collapse = "|"),
        #       #   PT,
        #       #   ignore.case = TRUE
        #       # )
        #     ][,
        #       PT := toupper(PT)
        #     ] %>%
        #       unique()
        #
        #   }
        # }

        out
      }
    )

  # browser()

  #
  #   dat1 <- dat_files$DRUG[
  #     ,
  #     .(sub = .(.SD)),
  #     by = .(PRIMARYID, CASEID)
  #   ][, nrow := sapply(sub, nrow)]
  #
  #
  #   dat2 <- dat_files$REACT[
  #     ,
  #     .(sub = .(.SD)),
  #     by = .(PRIMARYID, CASEID)
  #   ][, nrow := sapply(sub, nrow)]

  out <- merge(
    dat_files$DRUG,
    dat_files$REACT,
    by = c("PRIMARYID", "CASEID"),
    allow.cartesian = TRUE
  )[,
    QUARTER := this_quarter_name
  ] %>%
    {
      # add PROD_AI = NA if doesn't exist
      if (is.null(.$PROD_AI)) .[, PROD_AI := NA]
      else .
    }

  out
}


all_data_list <- files_faers %>%
  # replace lapply by future.lapply for parallel
  lapply(read_data_faers_files, future.seed = TRUE)


all_data_combined <- data.table::rbindlist(
  all_data_list,
  use.names = TRUE
)
#[
#   toupper(PT) %in% toupper(important_pt)
# ]

faers22q3raw <- all_data_combined[
  # ROLE_COD %in% c("PS", "SS")
  ROLE_COD == "PS"
][,
  .(CASEID, DRUGNAME, PT)
][,
  DRUGNAME := str_to_sentence(DRUGNAME)
][
  PT %in% important_pt
] |>
  setnames(
    old = c("PT", "DRUGNAME"),
    new = c("AE", "DRUG")
  )

usethis::use_data(faers22q3raw, compress = "xz", overwrite = TRUE)

# data_final %>%
#   object.size()

# data_final[
#   ,
#   .(ncases = length(unique(CASEID))),
#   keyby = .(DRUGNAME, PT)
# ]

# # all_data_combined %>% nrow()
# # map DRUGNAMEs by PROD_AI
# uniq_drugname_prod_ai <- all_data_combined[
#   !is.na(PROD_AI)
# ][,
#   .(DRUGNAME_list = .(unique(DRUGNAME))),
#   by = PROD_AI
# ][
#   ,
#   DRUGNAME_processed := map2_chr(
#     DRUGNAME_list,
#     PROD_AI,
#     function(this_names, this_prod_ai) {
#       # browser()
#       is_important <- sapply(
#         important_drug,
#         function(this_important_drug) {
#           # browser()
#           # this_names %>%
#           #   str_subset(this_important_drug)
#           #
#           grepl(
#             this_important_drug,
#             this_names,
#             ignore.case = TRUE
#           ) %>%
#             any()
#         },
#         USE.NAMES = TRUE
#       )
#
#       out <-  if (any(is_important)) {
#         is_important %>%
#           which() %>%
#           names() %>%
#           .[1]
#       } else {
#         this_prod_ai
#       }
#       # if (is.na(out)) out <- "Other_drug"
#
#       toupper(out)
#     }
#   )
# ][
#   ,
#   .(DRUGNAME_list, DRUGNAME_processed)
# ] %>%
# # [,
# #   DRUGNAME_list := unlist(DRUGNAME_list)
# # ] |>
# #   data.table::setnames(
# #     "DRUGNAME_list",
# #       "DRUGNAME"
# #     ) |>
# #     unique()
# #   #
# #   # #%>%
# unnest(DRUGNAME_list) %>%
#   as_tibble() %>%
#   rename(DRUGNAME = DRUGNAME_list) %>%
#   distinct(DRUGNAME, .keep_all = TRUE) %>%
#   data.table::setDT()
#
# # drug_name_vec <- drugname_by_prod_ai$DRUGNAME_processed %>%
# #   setNames(drugname_by_prod_ai$DRUGNAME)
# #
# # all_data_combined_mapped <- data.table::as.data.table(
# #   all_data_combined
# # )[, DRUGNAME_processed := drug_name_vec[DRUGNAME]]
#
# gc()
#
# DT <- `[`
#
# tmp <- left_join(
#   all_data_combined,
#   drugname_by_prod_ai,
#   by = "DRUGNAME"
# )
#
# all_data_combined_mapped <- merge(
#   x = all_data_combined,
#   y = uniq_drugname_prod_ai,
#   by = c("DRUGNAME"),
#   all.x = TRUE
# )[,
#   PROD_AI := NULL
# ][
#   !is.na(DRUGNAME_processed),
#   # DRUGNAME_processed := "Other_drug"
# ][
#   PT %in% toupper(important_pt)
# ] |>
#   unique() |>
#   DT(
#     (DRUGNAME_processed == "Other_drug") |
#       (
#         grepl(
#           "statin",
#           DRUGNAME_processed,
#           ignore.case = TRUE
#         ) &
#           ROLE_COD %in% c("PS", "SS")
#       )
#   )

# save the mapped data as R data object
save_data_name <- glue::glue(
  "{process_data_dir}/all_faers_files_combined.RDS"
)
saveRDS(all_data_combined_mapped, save_data_name, compress = "xz")


# # new table : with all 2012 - 2020 quarters
# # recreated table : 2014Q3 - 2017Q1
# old_version_quarters <- expand_grid(
#   year = 2014:2017,
#   quarter = paste0("Q", 1:4)
# ) %>%
#   mutate(year_quarter = paste0(year, quarter)) %>%
#   pull(year_quarter) %>%
#   setdiff(
#     c(paste0("2014Q", 1:2), paste0("2017Q", 2:4))
#   )
#
# # create contingency table and write as excel
# cont_mat_list <- list(
#   new = all_data_combined_mapped,
#   recreated = all_data_combined_mapped[
#     QUARTER %in% old_version_quarters
#   ]
# ) %>%
#   lapply(
#     function(this_dat) {
#       this_dat %>%
#         Matrix.utils::dMcast(
#           PT ~ DRUGNAME_processed,
#           fun.aggregate = "length"
#         ) %>%
#         magrittr::set_colnames(
#           colnames(.) %>%
#             str_remove("^DRUGNAME_processed")
#         ) %>%
#         data.matrix() %>%
#         as_tibble(rownames = "PT")
#     }
#   )
#
# cont_mat_list$old <- tab_s4
#
# cont_mat_list %>%
#   writexl::write_xlsx(
#     glue::glue("{process_data_dir}/table_s4_new.xlsx")
#   )
