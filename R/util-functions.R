`:=` <- data.table::`:=`

between <- data.table::between
`%between%` <- data.table::`%between%`

set_attr <- `attr<-`

set_dimnames <- `dimnames<-`

set_rownames <- magrittr::set_rownames

set_colnames <- magrittr::set_colnames

format_pval_ <- function(pv, digits = 2, eps = 1 / (10^(digits + 1)), ...) {
  . <- NULL

  pv %>%
    {
      ifelse(
        . >= 0.90,
        ">0.90",
        paste0("=", format.pval(., digits = digits, eps = eps))
      )
    } %>%
    gsub("=<", "<", .) %>%
    gsub("=>", ">", .)
}

scientific_10 <- function(n, digits = 2, ...) {
  . <- NULL

  out <-
    # Transforms the number into scientific notation even if small
    format(n, scientific = TRUE, digits = digits, ...) %>%
    # Replace e with 10^
    sub("e", "x10^", .) %>%
    # Remove + symbol and leading zeros on expoent, if > 1
    sub("\\+0?", "", .) %>%
    # Leaves - symbol but removes leading zeros on expoent, if < 1
    sub("-0?", "-", .)

  out
}


# process zero inflation such that certain
# places never get zero inflated
# zmat: observed 0-1 tables indicating zero inflation
# 1-means zero-inflation, 0-means no zero-inflation
.process_zero_inflation <- function(zmat,
                                    no_zero_infl_idx = NULL,
                                    n_row = nrow(zmat),
                                    n_col = ncol(zmat)) {
  if (!is.null(no_zero_infl_idx)) {
    for (kk in no_zero_infl_idx) {
      ii <- if (kk[1] != 0) kk[1] else 1:n_row
      jj <- if (kk[2] != 0) kk[2] else 1:n_col
      zmat[ii, jj] <- 0
    }
  }
  zmat
}

# returns 0 if the numerator is zero,
# other wise returns the ratio
.safe_divide <- function(x, y) {
  ifelse(
    x == 0,
    0,
    x / y
  )
}



# orders w.r.t. given names,
# then remove names
.order_then_unname <- function(obj, names, dimnames) {
  . <- NULL
  out <- unname(obj)
  if (is.matrix(obj)) {
    if (!is.null(dimnames(obj))) {
      out <- obj[dimnames[[1]], dimnames[[2]]] %>%
        unname(.)
    }
  }
}
