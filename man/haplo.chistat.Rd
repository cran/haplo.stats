\name{haplo.chistat}
\alias{haplo.chistat}
\title{
  Calculate a score test statistic for haplotypes
}
\description{
Calculate a score test statistic for haplotypes, given case-control
status.  For use in haplo.scan
}
\usage{
haplo.chistat(h1, h2, post, y, nrep)
}
\arguments{

\item{h1}{
Vector of codes for the first of a haplotype pair per subject
}
\item{h2}{
Vector of codes for the second of a haplotype pair per subject
}
\item{post}{
Posterior probability of the haplotype pair
}
\item{y}{
The case-control status (1 for cases; 0 for controls)
}
\item{nrep}{
Number of haplotype pairs for each subject
}
}
\value{
The score test statistic for the set of haplotypes
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{

}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
