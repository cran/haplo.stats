## $Id: sr.strsplit.q,v 1.1 2005/11/21 20:10:11 lunde Exp $

## sr.strsplit splits a vector of strings on the character string provided in
## sep.  By default, since sr.strsplit has no idea about the vector lengths
## after splitting, the return value is a list.  sr.strsplit was intended to
## separate strings by one character (space, /, #, etc.).  Note: using a sep
## of "" produces different results in S-PLUS and R.

## If collapse is TRUE, then sr.strsplit will attempt to 'collapse' the
## results into a matrix.  This will only succeed, when all of the input
## strings split into the same number of sub-strings.

## If drop is FALSE and the length of strings is 1, sr.strplit will keep the
## result in either list or matrix format (depending on the value of
## collapse).  By default, if strings as length 1, sr.strsplit's result is a
## simple vector.

## Differences between S-PLUS and R:
## - multicharacter sep values (e.g. "   " and "abc") are honored in R, but
##   only the first character is honored in S-PLUS
## - After splitting, S-PLUS removes leading and trailing whitespace
##   characters (space, tab, and newlines)

## Eric Lunde, 2005-11-21
## Mayo Clinic, Division of Biostatistics
## srlocal package development

sr.strsplit <- function(strings, sep = " ", collapse = FALSE, drop = TRUE) {
  result <- as.character(strings)

  if(is.R()) {
    result <- strsplit(result, split = sep, fixed = TRUE)
  } else {
    result <- sapply(result, splitString, sep = sep, trim = FALSE,
                     simplify = FALSE)
  }

  if(collapse) {
    lengths <- unlist(lapply(result, length))
    if(all(lengths == lengths[1])) {
      result <- matrix(unlist(result), nrow = length(result),
                       ncol = lengths[1], byrow = TRUE)
    } else {
      warning("Lengths are not the same, cannot collapse")
    }
  }

  if(drop) {
    if(is.list(result)) {
      if(length(result) == 1) {
        result <- result[[1]]
      }      
    } else if(is.matrix(result)) {
      if(nrow(result) == 1) {
        result <- result[1, ]
      }
    }
  }

  return (result)
}
