\name{geno.count.pairs}
\alias{geno.count.pairs}
\title{
  Counts of Total Haplotype Pairs Produced by Genotypes 
}
\description{
Provide a count of all possible haplotype pairs for each subject,
according to the phenotypes in the rows of the geno matrix.  
The count for each row includes the count for complete phenotypes, as 
well as possible haplotype pairs for phenotypes where there are
missing alleles at any of the loci.
}
\usage{
geno.count.pairs(geno)
}
\arguments{
\item{geno}{
Matrix of alleles, such that each locus has a  pair  of adjacent 
columns of alleles,  and the order of columns corresponds to the 
order of loci  on a chromosome.   If there are K loci, then geno 
has  2*K  columns.  Rows represent all observed alleles for each 
subject, their phenotype.
}

}
\value{
Vector where each element gives a count of the number haplotype pairs 
that are consistent with a  subject's phenotype, where a phenotype may
include 0, 1, or 2 missing alleles at any locus.
}
\section{Side Effects}{

}
\details{
When a subject has no missing alleles, and has h heterozygous sites, 
there are 2**(h-1) haplotype pairs that are possible ('**'=power).
For loci with missing alleles, we consider all possible pairs of alleles 
at those loci.  Suppose that there are M loci with missing alleles, and 
let the vector V have values 1 or 0 acccording to whether these loci 
are imputed to be heterozygous or homozygous, respectively. The length
of V is M.  The total number of  possible states of V is 
2**M. Suppose that  the vector W, also of length M, provides a count
of the number of possible heterozygous/homozygous states at the loci
with missing data. For example, if one allele  is missing, and there
are K possible alleles at that locus, then there can be one homozygous
and (K-1) heterozygous genotypes. If two alleles are missing, there
can be K homozygous and K(K-1)/2 heterozygous genotypes. Suppose the
function H(h+V) counts the total number of heterozygous sites among
the loci without missing data (of which h are heterozygous) and the
imputed loci (represented by the  vector V). Then, the total number of
possible pairs of haplotypes can be respresented as SUM(W*H(h+V)),
where the sum is over all possible values for the vector V.
}
\seealso{
\code{\link{haplo.em}},
\code{\link{summaryGeno}}
}
\examples{
setupData(hla.demo)
geno <- hla.demo[,c(17,18,21:24)]
geno <- geno.recode(geno)$grec
count.geno <- geno.count.pairs(geno)
print(count.geno)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
