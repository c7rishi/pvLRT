#' FDA Statin dataset with 46 adverse events
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FAERS databases for the
#' quarters 2014Q3 - 2020Q4
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 46 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that drug/AE pair and detected
#' in the FDA FAERS database during 2014Q3 - 2020Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin
#'
#' The 46 adverse events presented across the rows are considered
#' significant by FDA.
#'
#' @seealso
#' \link{statin}, \link{statin1491}, \link{gbca}
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin46"


#' FDA Statin dataset with 1491 adverse events
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FAERS databases for the
#' quarters 2014Q3 - 2020Q4
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 1491 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that drug/AE pair and detected
#' in the FDA FAERS database during 2014Q3 - 2020Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin
#'
#' The 1491 AEs stored in the dataset represent the intersection
#' of adverse events of the statin class of drugs and the GBCA drugs
#'
#' @seealso
#' \link{statin46}, \link{statin}, \link{gbca}
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin1491"


#' FDA Statin dataset with 6039 adverse events
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FAERS databases for the
#' quarters 2014Q3 - 2020Q4
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 6039 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that drug/AE pair and detected
#' in the FDA FAERS database during 2014Q3 - 2020Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin
#'
#' Corresponding to all 6039 observed adverse events (PTs) observed in statins
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin"


#' FDA lovastatin dataset
#'
#'
#' Lovastatin Drug-AE dataset (contingency table) obtained from the FAERS
#' databases for the quarters 2014Q3 - 2020Q4
#'
#' Drugs are stored across columns, and the 46 AEs are stored across the rows
#'
#'
#' Includes 1 column for lovastatin, and one combined column for all other drugs.
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"lovastatin"


#' FDA GBCA dataset with all observed 1707 adverse events
#'
#'
#' Drug-AE dataset (contingency table) obtained from the FAERS
#' databases for the quarters 2014Q3 - 2020Q4
#'
#' Drugs are stored across columns, and the 1707 AEs are stored across the rows
#'
#' the dataset includes 6 GBCA drugs:
#'
#' gadobenate, gadobutrol, gadodiamide, gadofosveset, gadopentetate, gadoterate,
#' gadoteridol, gadoversetamide, gadoxetate
#'
#' Corresponding to all 1707 observed adverse events (PTs) as curated in FAERS database.
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"gbca"
