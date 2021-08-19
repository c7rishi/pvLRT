`:=` <- data.table::`:=`

set_attr <- `attr<-`

set_dimnames <- `dimnames<-`

set_rownames <- magrittr::set_rownames

set_colnames <- magrittr::set_colnames

format_pval_ <- function(pv, digits = 2, eps = 1/(10^(digits+1)), ...) {
  pv %>%
    {ifelse(
      . >= 0.90,
      ">0.90",
      paste0("=", format.pval(., digits = digits, eps = eps))
    )}
}

scientific_10 <- function(n, digits = 2, ...) {
out <-
  #Transforms the number into scientific notation even if small
  format(n, scientific = TRUE, digits = digits, ...) %>%
  #Replace e with 10^
  sub("e", "x10^", .) %>%
  #Remove + symbol and leading zeros on expoent, if > 1
  sub("\\+0?", "", .) %>%
  #Leaves - symbol but removes leading zeros on expoent, if < 1
  sub("-0?", "-", .)

out
}
