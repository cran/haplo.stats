\name{plot.haplo.score}
\alias{plot.haplo.score}
\title{
  Plot Haplotype Frequencies versus Haplotype Score Statistics
}
\description{
Method function to plot a class of type haplo.score
}
\usage{
plot.haplo.score(x, ...)
}
\arguments{
\item{x}{
The object returned from haplo.score (which has class haplo.score).
}
\item{\dots }{
Dynamic parameter for the values of additional parameters for 
the plot method.  
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
This is a plot method function used to plot haplotype frequencies on
the x-axis and haplotype-specific scores on the y-axis. Because
haplo.score is a class, the generic plot function 
can be used, which in turn calls this plot.haplo.score function.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
"Score tests for association of traits with haplotypes when
linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}
\seealso{
haplo.score
}
\examples{
  setupData(hla.demo)
  geno <- as.matrix(hla.demo[,c(17,18,21:24)])
  keep <- !apply(is.na(geno) | geno==0, 1, any)
  hla.demo <- hla.demo[keep,]
  geno <- geno[keep,]
  attach(hla.demo)
  label <- c("DQB","DRB","B")
 
# For quantitative, normally distributed trait:

  score.gaus <- haplo.score(resp, geno, locus.label=label, 
                            trait.type = "gaussian")

  plot.haplo.score(score.gaus)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
