\name{print.haplo.score}
\alias{print.haplo.score}
\title{
  Print a haplo.score object
}
\description{
Method function to print a class of type haplo.score
}
\usage{
print.haplo.score(x, digits, ...)
}
\arguments{
\item{x}{
The object returned from haplo.score (which has class haplo.score).
}
\item{digits}{
Number of digits to round the numeric output.
}
\item{\dots }{
Dynamic parameter for the values of additional parameters for 
the print method.  
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
This is a print method function used to print information from
haplo.score class, with haplotype-specific information given in a
table. Because haplo.score is a class, the generic print function 
can be used, which in turn calls this print.haplo.score function.
}
\seealso{
haplo.score
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
