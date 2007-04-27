\name{print.haplo.group}
\alias{print.haplo.group}
\title{
  Print a haplo.group object
}
\description{
Method function to print a class of type haplo.group
}
\usage{
print.haplo.group(x, digits=max(options()$digits-2, 5), nlines=NULL, ...)
}
\arguments{
\item{x}{
The object returned from haplo.group (which has old class haplo.group).
}
\item{digits}{
Set the number of significant digits to print for haplotype probabilities.
}
\item{nlines}{
For shorter output, print first 1:nlines rows of the large data frame
}
\item{\dots}{
Optional arguments for the print method
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
This is a print method function used to print information from the
haplo.group class, with haplotype-specific information given in a
table. Because haplo.group is a class, the generic print function 
can be used, which in turn calls this print.haplo.group function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland  GA.
Expected haplotype frequencies for association of traits with 
haplotypes when linkage phase is ambiguous. 
Submitted to Amer J Hum Genet.
}
\seealso{
haplo.score, haplo.group, haplo.em
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
