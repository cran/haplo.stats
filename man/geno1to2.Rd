\name{geno1to2}
\alias{geno1to2}
\title{
  convert genotype matrix from 1-column 2-column 
}
\description{
convert 1-column genotype matrix to 2-column genotype matrix,
converting from a minor allele count (0,1,2) to (1/1, 1/2, 2/2) where 2 is the
minor allele. (not supported for x-linked markers)
}
\usage{
geno1to2(geno, locus.label=NULL)
}
\arguments{
\item{geno}{
1-column representation of genotype matrix for 2-allele loci.  Values
are 0, 1, or 2, usually the count of minor alleles
}
\item{locus.label }{
Vector of labels for loci, If a locus name is "A", its columns
will be "A.1" and "A.2"
}
}
\value{
a 2-column genotype matrix
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
geno1 <- matrix(c(0,0,1,
                  1,0,2,
                  2,1,0), ncol=3, byrow=TRUE)
geno1to2(geno1, locus.label=c("A", "B", "C"))

## demonstrate how NA and 3 will be coded
geno1[1,3] <- NA
geno1[1,1] <- 3
geno1to2(geno1)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.