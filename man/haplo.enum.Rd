\name{haplo.enum}
\alias{haplo.enum}
\title{
  Enumerate all possible pairs of haplotypes that are consistent with 
a set of un-phased multilocus markers
}
\description{
Given subject un-phased genotype hmat, enumerate all possible pairs 
of haplotypes, and return enumerated pairs in matrices h1 and h2.
}
\usage{
haplo.enum(hmat)
}
\arguments{

\item{hmat}{
A genotype vector of length 2*K (K = number of loci).  When used in 
haplo.em, it is a single row of a genotype matrix.
}
}
\value{
List with components:

\item{h1}{
A matrix of enumerated haplotypes.  If there are N
enumerations, h1 will have dimension N x K.
}
\item{h2}{
Similar to h1, a matrix of enumerated haplotypes for the second 
members of the pairs of haplotypes.  Haplotype pairs in h1 and h2
match by the same row number.
}
}
\section{Side Effects}{

}
\details{
For a pair of haplotypes, if there are H sites that are
heterozygous, then there are 2 raised to (H-1) possible pairs to
enumerate. To achieve this, the algorithm moves across the
loci that are heterozygous  (after the 1st heterozygous
locus), flipping alleles at heterozygous locations to
enumerate all possible pairs of haplotpes, and appending
results as rows of the output matrices h1, and h2.
}
\section{References}{

}
\seealso{
\code{\link{haplo.em}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
