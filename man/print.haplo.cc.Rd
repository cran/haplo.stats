\name{print.haplo.cc}
\alias{print.haplo.cc}
\title{
  Print a haplo.cc object
}
\description{
Display results for a haplotype analysis on a case-control study.  
}
\usage{
print.haplo.cc(x, order.by="score", digits=max(options()$digits-2, 5),
    nlines=NULL ...)
}
\arguments{
\item{x}{
A haplo.cc object, made by the haplo.cc function.
}
\item{order.by}{
Order the printed data frame by haplotype score (score), haplotype alleles
(haplotype), or haplotype frequency (freq).
}
\item{digits}{
Number of digits to display for the numeric columns of the data frame.
}
\item{nlines}{
Print the first nlines of the cc.df data frame of the haplo.cc object,
keeps output short if desired.
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{

}
\seealso{
\code{\link{haplo.cc}},
}
\examples{
## for a haplo.cc object named cc.test, 
##    order results by haplotype 
# print(cc.test, order.by="haplotype")
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
