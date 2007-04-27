\name{haplo.design}
\alias{haplo.design}
\title{
  Build a design matrix for haplotypes
}
\description{
Build a design matrix for haplotypes estimated from a haplo.em
object. 
}
\usage{
haplo.design(obj, haplo.effect="additive", hapcodes=NA, min.count=5, haplo.base=NA)
}
\arguments{
\item{obj}{
an object created from haplo.em
}
\item{haplo.effect }{
The "effect" pattern of haplotypes on the response. This
parameter determines the coding for scoring the haplotypes. 
Valid coding options for heterozygous and homozygous carriers
of a haplotype are "additive" (1, 2, respectively),
"dominant" (1,1, respectively), and "recessive" (0, 1,
respectively). 
}
\item{hapcodes}{
codes assigned in haplo.em, corresponding to the row numbers in the
obj\$haplotypes matrix
}
\item{min.count }{
The minimum number of estimated counts of the haplotype in the sample
in order for a haplotype to be included in the design matrix.  
}
\item{haplo.base }{
code for which haplotype will be the reference group, or to be
considered the baseline of a model.  The code is the row number of the
obj\$haplotypes matrix. This haplotype is removed from the design matrix.
}
}
\value{
Matrix of columns for haplotype effects. Column names are "hap.k"
where k is the row number of the unique haplotypes within the haplo.em
object's "haplotypes" item. 
}
\section{Side Effects}{

}
\details{
First a matrix is made for the possible haplotypes for each
person, coded for the haplo.effect, weighted by the posterior
probability of those possible haplotypes per person, and then 
collapsed back to a single row per person.
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
% Converted by Sd2Rd version 37351.