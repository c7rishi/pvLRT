#' FDA Statin dataset with 46 adverse events
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA FAERS database for the
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
#' obtained from the FDA FAERS database for the
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
#' obtained from the FDA FAERS database for the
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
#' Corresponding to all 6039 observed adverse events (AEs) observed in statins
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin"


#' FDA lovastatin dataset
#'
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA FAERS database for the
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
#' The dataset contains 1 column for the lovastatin drug,
#' and one column for all other drugs combined.
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"lovastatin"


#' FDA GBCA dataset with all observed 1707 adverse events
#'
#'
#' A Drug-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2014Q3 - 2020Q4
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 1707 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that drug/AE pair and detected
#' in the FDA FAERS database during 2014Q3 - 2020Q4.
#'
#' The dataset contains 6 Gadolinium-Based Contrast Agents (GBCAs) as drugs:
#'
#' gadobenate, gadobutrol, gadodiamide, gadofosveset, gadopentetate, gadoterate,
#' gadoteridol, gadoversetamide, gadoxetate
#'
#' Corresponding to all 1707 observed adverse events (AEs) as curated in FAERS database.
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"gbca"




#' FDA rotavirus vaccine dataset with 794 adverse events observed among
#' combined old (age >= 1 year) and  young (age < 1 year) individuals
#'
#'
#' A vaccine-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA VAERS database for the year 1999
#'
#'
#' @details
#' Data are stored in the forms of contingency table, with the
#' vaccines listed across the columns and the 794 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that vaccine/AE pair and detected
#' in the FDA VAERS database for the year 1999.
#'
#' The dataset contains two columns -- one for the rotavirus vaccine, and another
#' for other (37 vaccines combined).
#'
#' @source \url{https://vaers.hhs.gov/data/datasets.html}
"rv"


#' FDA rotavirus vaccine dataset with 727 adverse events observed among
#' "old" (non-infant; age >= 1 year) individuals
#'
#'
#' A vaccine-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA VAERS database for the year 1999
#'
#'
#' @details
#' Data are stored in the forms of contingency table, with the
#' vaccines listed across the columns and the 727 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts (from "old" individuals with age >= 1 year) associated
#' with that vaccine/AE pair and detected
#' in the FDA VAERS database for the year 1999.
#'
#' The dataset contains two columns -- one for the rotavirus vaccine, and another
#' for other (37 vaccines combined).
#'
#' @source \url{https://vaers.hhs.gov/data/datasets.html}
"rvold"


#' FDA rotavirus vaccine dataset with 346 adverse events observed among
#' young (infant -- 1 year) individuals
#'
#'
#' A vaccine-Adverse event (AE) count dataset (contingency table)
#' obtained from the FDA VAERS database for the year 1999
#'
#'
#' @details
#' Data are stored in the forms of contingency table, with the
#' vaccines listed across the columns and the 346 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts from young individuals with age < 1 year associated
#' with that vaccine/AE pair and detected
#' in the FDA VAERS database for the year 1999.
#'
#' The dataset contains two columns -- one for the rotavirus vaccine, and another
#' for other (37 vaccines combined).
#'
#' @source \url{https://vaers.hhs.gov/data/datasets.html}
"rvyoung"


