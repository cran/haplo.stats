\name{haplo.scan.obs}
\alias{haplo.scan.obs}
\title{
  For observed data, slide a fixed-width window over each marker locus
and scan all possible haplotypes within the window
}
\description{
Slide a fixed-width window over each marker locus and scan all
possible haplotypes within the window.  For all locus-subsets, 
calculate score test statistics on observed data, and save the necessary
information for calculating the test statistics in haplo.scan.sim.
}
\usage{
haplo.scan.obs(y, em.obj, width)
}
\arguments{
\item{y}{
Vector of trait values, must be 1 for cases and 0 for controls.
}
\item{em.obj }{
A haplo.em object containing necessary haplotype information
}
\item{width}{
Width of the slide window for scanning
}

}
\value{

\item{svec}{
Vector of locus-specific maximum test statistics
}
\item{save.lst}{
A list object containing haplotype information necessary for calculating
score test statistics on the simulated data sets in haplo.scan.sim
}
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{haplo.scan}},
\code{\link{haplo.scan.sim}},
\code{\link{haplo.chistat}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
