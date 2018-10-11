##############
# General utility functions
#
# Joshua Schwab
# jschwab77@berkeley.edu
#---------------------------

library(dplyr)
library(data.table)
library(openxlsx)
library(googlesheets)
library(matrixStats)

#avoid comparisons between character and logical
"%in%" <- function(x, table) {
  if ((is.character(x) & is.logical(table)) || (is.character(table) & is.logical(x))) {
    dput(x)
    dput(table)
    stop("bad classes using %in%") 
  } else {
    match(x, table, nomatch = 0) > 0
  }
}

NA_Date <- as.Date(NA)

#halfway between, regardless of which is later
HalfwayDate <- function(d1, d2) d1 + (d2 - d1) / 2

#d1 (strictly) before d2 and d1 is not NA 
DateBefore <- function(d1, d2, strict = T) { 
  stopifnot(class(d1) == "Date" && class(d2) == "Date")
  if (strict) {
    is.less <- !is.na(d1) & (d1 < d2)
  } else {
    is.less <- !is.na(d1) & (d1 <= d2)
  }
  stopifnot(!anyNA(is.less))
  return(is.less)
}

#make sure setnames does not result in multiple columns with the same name
setnames <- function(x, old, new) {
  data.table::setnames(x, old, new)
  if (any(duplicated(names(x)))) {
    stop("after setnames, x has multiple columns with the same name")
  }
  invisible(x)
}

