\name{sr.strsplit}
\alias{sr.strsplit}
\title{
  Split a vector of strings by a single character
}
\description{
Split a character vector using the sep argument.  Results will be a
vector, matrix, or list, depending on the value of collapse and drop.
}
\usage{
sr.strsplit(strings, sep=" ", collapse=FALSE, drop=TRUE)
}
\arguments{
\item{strings}{
character vector
}
\item{sep}{
single character to separate the strings
}
\item{collapse}{
If collapse is TRUE, then sr.strsplit will attempt to 'collapse' the
results into a matrix.  This will only succeed, when all of the input
strings split into the same number of sub-strings. 
}
\item{drop}{
If drop is FALSE and the length of strings is 1, sr.strplit will keep
the result in either list or matrix format (depending on the value of
collapse).  By default, if strings as length 1, sr.strsplit's result
is a simple vector.
}
}
\value{
a list, matrix, or vector of the split results of the strings supplied
}
\details{
Differences between S-PLUS and R:

- multicharacter sep values (e.g. "   " and "abc") are honored in R, but
  only the first character is honored in S-PLUS

- After splitting, S-PLUS removes leading and trailing whitespace
  characters (space, tab, and newlines)
}
\seealso{
See splitString for S-PLUS and strsplit for R
}
\examples{
x <- c("str:split" , "new:string")
split <- sr.strsplit(x, sep = ":")
## split has the value:
##[[1]]:
##[1] "str"   "split"
##
##[[2]]:
##[1] "new"    "string"

split <- sr.strsplit(x, sep = ":", collapse = TRUE)
## split has the value:
##      [,1]     [,2] 
##[1,] "str" "split" 
##[2,] "new" "string"

y <- "str:split"
split <- sr.strsplit(y, sep = ":")
## split has the value:
##[1] "str"   "split"
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
