library(tidyverse)

all_data <- here::here("data-raw/") %>%
  list.files(pattern = "\\.xlsx$", full.names = TRUE) %>%
  setNames(
    strsplit(., "/") %>%
      map_chr(tail, 1) %>%
      str_remove("\\.xlsx") %>%
      strsplit(" ") %>%
      map_chr(head, 1) %>%
      str_remove("_all") %>%
      str_remove("_") %>%
      tolower()
  ) %>%
  lapply(readxl::read_excel) %>%
  lapply(
    . %>%
      select(-pt_marginal) %>%
      filter(!toupper(pt) == "TOTAL") %>%
      # filter(!grepl("Total", pt, ignore.case = TRUE)) %>%
      as.data.frame() %>%
      magrittr::set_rownames(.$pt) %>%
      select(-pt) %>%
      data.matrix() %>%
      magrittr::set_rownames(
        rownames(.) %>%
          str_to_title()
      ) %>%
      magrittr::set_colnames(
        colnames(.) %>%
          str_to_title()
      )
  )


## code to prepare `DATASET` dataset goes here
for (xx in names(all_data)) assign(xx, all_data[[xx]])


usethis::use_data(statin, overwrite = TRUE, compress = "xz")
usethis::use_data(statin1491, overwrite = TRUE, compress = "xz")
usethis::use_data(statin46, overwrite = TRUE, compress = "xz")
usethis::use_data(gbca, overwrite = TRUE, compress = "xz")
usethis::use_data(rv, overwrite = TRUE, compress = "xz")
usethis::use_data(rvyoung, overwrite = TRUE, compress = "xz")
usethis::use_data(rvold, overwrite = TRUE, compress = "xz")
