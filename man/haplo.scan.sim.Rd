\name{haplo.scan.sim}
\alias{haplo.scan.sim}
\title{
  For simulated data, slide a fixed-width window over each marker locus
and scan all possible haplotypes within the window
}
\description{
For simulated data, slide a fixed-width window over each marker locus
and scan all possible haplotypes within the window
}
\usage{
haplo.scan.sim(y.reord, save.lst, nloci)
}
\arguments{
\item{y.reord }{
Re-ordered vector of case/control status, must be 1 for cases and 0 for controls.
}
\item{save.lst }{
List of haplotype information saved from haplo.scan.obs
}
\item{nloci}{
Number of loci in the whole region
}

}
\value{
Vector of the locus-specific test statistics based on
simulated case/control status
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{haplo.scan}},
\code{\link{haplo.scan.obs}},
\code{\link{haplo.chistat}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
