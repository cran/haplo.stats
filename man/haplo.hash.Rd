\name{haplo.hash}
\alias{haplo.hash}
\title{
  Integer Rank Codes for Haplotypes
}
\description{
Create a vector of integer codes for the input  matrix  of
haplotypes.   The haplotypes in the input matrix are
converted to character strings, and if there are C unique
strings, the integer codes for the haplotypes will be 1,
2, ..., C.
}
\usage{
haplo.hash(hap)
}
\arguments{

\item{hap}{
A matrix of haplotypes. If there are N haplotypes for K
loci, hap have dimensions N x K.
}
}
\value{
List with elements:

\item{hash}{
Vector of integer codes for the input data  (hap).  The
value of hash is the row number of the unique haplotypes
given in the returned matrix hap.mtx.
}
\item{hap.mtx}{
Matrix of unique haplotypes.
}
}
\section{Side Effects}{

}
\details{
The alleles that make  up  each  row  in  hap  are  pasted
together  as character strings, and the unique strings are
sorted so that the rank order of  the  sorted  strings  is
used as the integer code for the unique haplotypes.
}
\section{References}{

}
\seealso{
haplo.em
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
