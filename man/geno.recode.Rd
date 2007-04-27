\name{geno.recode}
\alias{geno.recode}
\title{
  Recode Genotypes
}
\description{
For all loci as pairs of columns in a matrix, recode
alleles
}
\usage{
geno.recode(geno, miss.val=0)
}
\arguments{
\item{geno}{
Matrix of alleles, such that each locus has a  pair  of
adjacent  columns  of  alleles.  If there are K loci, then
ncol(geno) = 2*K. Rows represent alleles for each subject.
}
\item{miss.val }{
Vector of codes for missing values of alleles.
}
}
\value{
List with components:

\item{grec}{
Matrix of recoded alleles - see allele.recode
}
\item{alist}{
List of allele labels.  For K loci, there are K
components in the list, and the kth  component is a vector
of sorted unique allele labels for the kth locus.
}
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
allele.recode
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
