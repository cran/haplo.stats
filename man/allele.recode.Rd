\name{allele.recode}
\alias{allele.recode}
\title{
  Recode allele values to integer ranks
}
\description{
Genotypes for subjects represented by a pair  of  vectors,
with the vectors containing allele values (either numeric,
factor, or character), are recoded to the  rank  order  of
allele values.
}
\usage{
allele.recode(a1, a2, miss.val=NA)
}
\arguments{
\item{a1}{
Vector of "first" alleles.
}
\item{a2}{
Vector of "second" alleles.
}
\item{miss.val }{
Vector of missing value codes for alleles.
}
}
\value{
List with components:

\item{a1}{
Vector of recoded "first" alleles.
}
\item{a2}{
Recode of recoded "second" alleles.
}
\item{allele.label}{
Vector of labels for unique alleles.
}
}
\section{Side Effects}{

}
\details{
If alleles are numeric, they are recoded to the rank order
of the alleles.  If the alleles are factor or character,
they are recoded to interger values that correspond to the
indices of the sorted values of the unique alleles, but
sorted as character values.
}
\section{References}{

}
\seealso{
geno.recode
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
