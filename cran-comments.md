# pvLRT v0.5

## Resubmission

This is a re-submission of the package pvLRT. In this version we have made the following changes:

-   Add a new function 'convert_raw_to_contin_table' to create contingency table from raw Drug-AE incidence data frames.
-   fix typos in the docuementation texts

## Test environments

-   local Windows 10 install, R 4.2.2
-   Ubuntu 20.04.3 LTS, R 4.2.2
-   win-builder (devel, release, oldrelease)

## R CMD Check Results

-   There were no ERRORs or WARNINGs or NOTEs.

# pvLRT v0.4

## Resubmission

This is a re-submission of the package pvLRT. In this version we have made the following changes:

-   Add a new function 'extract_run_time' for pvlrt objects. Run times are also shown when a pvlrt object is printed.
-   Add argument 'fill_gradient_range' to all plotting functions.
-   use ggfittext to fit texts inside pvlrt plots
-   add rotavirus vaccine data
-   Default colors are unified across all plots (bubble, bar, heatmap)

## Test environments

-   local Windows 10 install, R 4.1.2
-   Ubuntu 20.04.3 LTS, R 4.1.2
-   win-builder (devel, release, oldrelease)

## R CMD Check Results

-   There were no ERRORs or WARNINGs or NOTEs.

# pvLRT v0.3

This is a CRAN submission for the (new) R package pvLRT v0.3, which provides a suite of functions implementing likelihood ratio test-based approaches to pharmacovigilance. We have ensured that all examples take \< 10s time to pass the CRAN check in this version. (One example failed this test during the initial CRAN submission of v0.2.)

In addition, we have taken into account all the helpful comments and suggestions provided by Gregor Seyer on behalf of CRAN during the submission of pvLRT v0.2.1. We greatly appreciate the comments; we have made all the necessary changes in v0.3, and our point-by-point responses are provided below.

## Point-by-point responses addressing Greg Seyer's comments on pvLRT v0.2.1

-   Please consider the use of the [Authors\@R](mailto:Authors@R){.email} field where you can declare Maintainer, Authors and Contributors with their appropriate roles with person() calls. e.g. something like: [Authors\@R](mailto:Authors@R){.email}: c(person("Alice", "Developer", role = c("aut", "cre","cph"), email = "[alice.developer\@some.domain.net](mailto:alice.developer@some.domain.net){.email}"), person("Bob", "Dev", role = "aut") )

--\> Response: Thanks for the suggestion. In v0.3, we have updated the DESCRIPTION with an appropriately populated [Authors\@R](mailto:Authors@R){.email} filed.

-   Please reduce the length of the title to less than 65 characters.

--\> Response: We have updated the title to: "Likelihood Ratio Test-Based Approaches to Pharmacovigilance" (59 characters)

\*Please omit the redundant "An R package" from your description.

--\> Response: We have made this change.

-   If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) \<doi:...\> authors (year) \<arXiv:...\> authors (year, ISBN:...) or if those are not available: \<https:...\> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

--\> Response: Thanks, we don't yet have references describing the methods in our package, so we have not made any change on this point.

-   Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar) Missing Rd-tags in up to 11 .Rd files, e.g.: as.matrix.pvlrt.Rd: \value extract_AE_names.Rd: \value extract_lrstat_matrix.Rd: \value heatmap_pvlrt.Rd: \value logLik.pvlrt.Rd: \value lrt_poisson.Rd: \value ...

--\> Response: Thank you for the suggestion. In pvLRT v0.3, we have added a \value tag to each exported method, which describes the output structure and what the output means.

\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("\# Not run:") as a warning for the user. Does not seem necessary.

Please unwrap the examples if they are executable in \< 5 sec, or replace \dontrun{} with \donttest{}.

--\> Response: Thank you for pointing this out. We have removed the \dontrun{} wrappers throughout. We have used fewer iterations in the examples wherever appropriate (with added comments) to ensure that all but one example take \<5 secs to run. One example under the function `pvlrt` still takes \>5 secs; we have wrapped it under \donttest{}.

## Test environments

-   local Windows 10 install, R 4.1.2
-   Ubuntu 20.04.3 LTS, R 4.1.2
-   win-builder (devel, release, oldrelease)

## R CMD Check Results

-   There were no ERRORs or WARNINGs.

-   There were no NOTEs in local Windows 10 and Ubuntu 20.04.3 checks

-   There was one NOTE in win-builder:

-   checking CRAN incoming feasibility ... NOTE Maintainer: 'Saptarshi Chakraborty [chakra.saptarshi\@gmail.com](mailto:chakra.saptarshi@gmail.com){.email}'

New submission

Possibly misspelled words in DESCRIPTION: Pharmacovigilance (2:50) pharmacovigilance (11:71)

We believe this NOTE is a false alarm.

# pvLRT v0.2.1

This is a new CRAN submission for R package pvLRT, which provides a suite of functions implementing likelihood ratio test-based approaches to pharmacovigilance. In this version we have made sure that:

-   all examples take \< 10s time to pass the CRAN check. One example failed this test during initial CRAN submission.

## Test environments

-   local Windows 10 install, R 4.1.2
-   Ubuntu 20.04.3 LTS, R 4.1.2
-   win-builder (devel, release, oldrelease)

## R CMD Check Results

-   There were no ERRORs or WARNINGs.

-   There were no NOTEs in local Windows 10 and Ubuntu 20.04.3 checks

-   There was one note in win-builder:

-   checking CRAN incoming feasibility ... NOTE Maintainer: 'Saptarshi Chakraborty [chakra.saptarshi\@gmail.com](mailto:chakra.saptarshi@gmail.com){.email}'

New submission

Possibly mis-spelled words in DESCRIPTION: Pharmacovigilance (3:23) pharmacovigilance (8:97)

We believe this NOTE is a false alarm, as 'pharmacovigilance' refers to the field relating to the collection, detection, assessment, monitoring, and prevention of adverse effects with pharmaceutical products (Wikipedia; <https://en.wikipedia.org/wiki/Pharmacovigilance>).
