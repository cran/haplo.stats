\name{haplo.group}
\alias{haplo.group}
\title{
  Frequencies for Haplotypes by Grouping Variable
}
\description{
Calculate maximum likelihood estimates of haplotype probabilities for 
the entire dataset and separately for each subset defined by the
levels of a group variable.  Only autosomal loci are considered.
}
\usage{
haplo.group(group, geno, locus.label=NA, miss.val=0, 
            control=haplo.em.control())
}
\arguments{
\item{group}{
Group can be of logical, numeric, character, or factor class type.
}
\item{geno}{
Matrix of alleles, such that each locus has a pair of
adjacent columns of alleles, and the order of columns
corresponds to the order of loci on a chromosome.  If
there are K loci, then geno has 2*K columns. Rows represent all
observed alleles for each subject.  
}
\item{locus.label }{
Vector of labels for loci, of length K (see definition of geno matrix).
}
\item{miss.val }{
Vector of codes for allele missing values.
}
\item{control}{
list of control parameters for haplo.em (see haplo.em.control).
}
}
\value{

\item{ }{
A list as an object of the haplo.group class.  The three
elements of the list are described below.
}
\item{group.df}{
A data frame with the columns described as follows.
-haplotype:
Names for the K columns for the K alleles in the haplotypes.
-total:
Estimated frequencies for haplotypes from the total sample.
-group.name.i:
Estimated haplotype frequencies for the haplotype
if it occurs in the group referenced by 'i'.  Frequency is NA 
if it doesn't occur for the group.  The column name is the actual 
variable name joined with the ith level of that variable.
}
\item{group.count}{
Vector containing the number of subjects for each level of the 
grouping variable.
}
\item{n.loci}{
Number of loci occuring in the geno matrix.
}
}
\section{Side Effects}{

}
\details{
Haplo.em is used to compute the maximum likelihood estimates of 
the haplotype frequencies for the total sample, then for each of 
the groups separately.  
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
"Score tests for association of traits with haplotypes when
linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}
\seealso{
print.haplo.group, haplo.em
}
\examples{
  setupData(hla.demo)
  geno <- as.matrix(hla.demo[,c(17,18,21:24)])

# remove any subjects with missing alleles for faster examples, 
# but you may keep them in practice
  keep <- !apply(is.na(geno) | geno==0, 1, any)
  hla.demo <- hla.demo[keep,]
  geno <- geno[keep,]
  attach(hla.demo)
  
  y.ord <- as.numeric(resp.cat)
  y.bin <-ifelse(y.ord==1,1,0)
  group.bin <- haplo.group(y.bin, geno, miss.val=0)
  print.haplo.group(group.bin)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
