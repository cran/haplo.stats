\name{print.haplo.score.merge}
\alias{print.haplo.score.merge}
\title{
  Print a haplo.score.merge object
}
\description{
Method function to print a class of type haplo.score.merge
}
\usage{
print.haplo.score.merge(x, order.by="score", all.haps=FALSE, 
                        digits=max(options()$digits-2, 5), ...)
}
\arguments{
\item{x}{
The object returned from haplo.score.merge (which has old class \{S\} 
haplo.score.merge).
}
\item{order.by}{
Column of the haplo.score.merge object by which to order the results.
}
\item{all.haps}{
Logical, if (T)rue prints a row for all haplotypes.  If (F)alse, the
default, only prints the haplotypes kept in haplo.score for modelling.
}
\item{digits}{
Set the number of significant digits to print for the numeric output.
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
This is a print method function used to print information from the
haplo.score.merge class. Because haplo.score.merge is a class, 
the generic print function can be used, which in turn calls this 
print.haplo.score.merge function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland  GA.
Expected haplotype frequencies for association of traits with 
haplotypes when linkage phase is ambiguous. 
Submitted to Amer J Hum Genet.
}
\seealso{
haplo.score.merge, haplo.score, haplo.group
}
\examples{
#see example for haplo.score.merge
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
