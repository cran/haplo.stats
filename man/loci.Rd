\name{loci}
\alias{loci}
\title{
  Create a group of locus objects from a genotype matrix, assign to
'model.matrix' class.  
}
\description{
The function makes each pair of columns a locus object, which
recodes alleles to numeric and saves the original alleles as an
attribute of the model.matrix. 
}
\usage{
loci(geno, locus.names, chrom.label=NULL, x.linked=FALSE, sex=NULL, 
     male.code="M", female.code="F", miss.val=NA, map=NA)
}
\arguments{
\item{geno}{
Matrix of alleles, such that each locus has a pair of adjacent
columns of alleles, and the order of columns corresponds to the
order of loci on a chromosome. If there are K loci, then
ncol(geno) = 2*K. Rows represent alleles for each subject.
}
\item{locus.names }{
A vector containing the locus name for each locus.
}
\item{chrom.label }{
Chromosome Label 
}
\item{x.linked }{
A logical value denoting whether the chromosome is X-linked.
}
\item{sex}{
A vector containing the sex of each individual. If x.linked=F then argum
ent sex is not required and may be left as the default value of NULL.
}
\item{male.code }{
The code denoting a male in the sex vector.
}
\item{female.code }{
The code denoting a female in the sex vector.
}
\item{miss.val }{
A vector of codes denoting missing values for the allele labels.  
Note that NA will always be treated as a missing value, and alleles
matching miss.val are assigned NA.  Also note that the original
missing value code for a specific individual can not be retrieved 
from the returned object.
}
\item{map}{
An optional chromosome map of class "cmap"
}
}
\value{
An object of class "model.matrix", with all alleles recoded to a
numeric value.  It contains the following attributes:

\item{locus.names}{
A vector of labels for the loci, of length nloci.
}
\item{map}{
Will be better defined later.
}
\item{x.linked}{
A logical value denoting whether the chromosome is X-linked.
}
\item{unique.alleles}{
The original allele labels are stored in the 'unique.alleles' attribute.
The ith item of the unique.alleles list is a vector of unique alleles
for the ith locus.
}
\item{male.code}{
The code denoting a male in the sex vector.
}
\item{female.code}{
The code denoting a female in the sex vector.
}
\item{chrom.label}{
Chromosome Label
}
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\note{
A matrix that contains all elements of mode character will be sorted
in alphabetic order.
}
\seealso{
locus, setupGeno
}
\examples{
# Create some loci to work with
a1 <- 1:6
a2 <- 7:12

b1 <- c("A","A","B","C","E","D")
b2 <-c("A","A","C","E","F","G")

c1 <- c("101","10","115","132","21","112")
c2 <- c("100","101","0","100","21","110")

myloci <- data.frame(a1,a2,b1,b2,c1,c2)
myloci <- loci(myloci, locus.names=c("A","B","C"),miss.val=c(0,NA))
myloci

attributes(myloci)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
