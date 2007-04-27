\name{seqhap}
\alias{seqhap}
\title{
  Sequential Haplotype Scan Association Analysis for Case-Control Data
}
\description{
Seqhap implements sequential haplotype scan methods to perform
association analyses for case-control data.  When evaluating each
locus, loci that contribute additional information to haplotype 
associations with disease status will be added sequentially. 
This conditional evaluation is based on the 
Mantel-Haenszel (MH) test.  Two sequential methods are provided, a
sequential haplotype method and a sequential summary method, as well as
results based on the traditional single-locus method.  Currently, seqhap
only works with bialleleic loci (single nucleotide polymorphisms, 
or SNPs) and binary traits.
}
\usage{
seqhap(y, geno, pos, n.sim=1000, weight=NULL, 
       mh.threshold=3.84, r2.threshold=0.95, haplo.freq.min=0.005, 
       miss.val=c(0, NA), control=haplo.em.control())
}
\arguments{
\item{y}{
vector of binary response (1=case, 0=control). The length is equal 
to the number of rows in geno. 
}
\item{geno}{
matrix of alleles, such that each locus has a pair of adjacent columns
of alleles, and the order of columns corresponds to the order of loci
on a chromosome. If there are K loci, then ncol(geno)=2*K. Rows
represent the alleles for each subject. Currently, only bi-allelic loci
(SNPs) are allowed. 
}
\item{pos}{
vector of physical positions (or relative physical positions) for
loci. If there are K loci, length(pos)=K. The scale (in kb, bp, or
etc.) doesn't affect the results.
}
\item{n.sim }{
number of permutations that will be used to calculate permutation p-values. The
default is 1000. If it is set to zero, permuted regional and
pointwise p-values will not be available, but the pointwise p-values
based on asymptotic chi-square distributions will still be provided. 
}
\item{weight}{
weights for observations (rows of geno matrix).
}
\item{mh.threshold }{
threshold for the Mantel-Haenszel statistic that evaluates whether a
locus contributes additional information of haplotype association to
disease, conditional on current haplotypes. The default is 3.84, which
is the 95th percentile of the chi-square distribution with 1 degree 
of freedom. 
}
\item{r2.threshold }{
threshold for a locus to be skipped. When scanning locus k, loci
with correlations r-squared (the square of the Pearson's correlation)
greater than r2.threshold with locus k will be ignored, so that the haplotype
growing process continues for markers that are further away from locus k.
}
\item{haplo.freq.min }{
the minimum haplotype frequency for a haplotype to be included in the
association tests. The haplotype frequency is based on the EM
algorithm that estimates haplotype frequencies independent of trait. 
}
\item{miss.val }{
 
vector of values that represent missing alleles.
}
\item{control}{
A list of parameters that control the EM algorithm for estimating 
haplotype frequencies when phase is unknown.  The list is created by 
the function haplo.em.control - see this function for more details. 
}
}
\value{
list with components:

\item{converge}{
indicator of convergence of the EM algorithm (see haplo.em); 
1 = converge, 0=failed
}
\item{locus.label}{
vector of labels for loci, of length K.
}
\item{inlist}{
matrix that shows which loci are combined for association analysis in
the sequential scan. The non-zero values of the kth row of inlist are
the indices of the loci combined when scanning locus k. 
}
\item{chi.stat}{
chi-square statistics of single-locus analysis.
}
\item{chi.p.point}{
permuted pointwise p-values of single-locus analysis.
}
\item{chi.p.region}{
permuted regional p-value of single-locus analysis.
}
\item{hap.stat}{
chi-square statistics of sequential haplotype analysis.
}
\item{hap.df}{
degrees of freedom of sequential haplotype analysis.
}
\item{hap.p.point}{
permuted pointwise p-values of sequential haplotype analysis.
}
\item{hap.p.region}{
permuted region p-value of sequential haplotype analysis.
}
\item{sum.stat}{
chi-square statistics of sequential summary analysis.
}
\item{sum.df}{
degrees of freedom of sequential summary analysis.
}
\item{sum.p.point}{
permuted pointwise p-values of sequential summary analysis.
}
\item{sum.p.region}{
permuted regional p-value of sequential summary analysis.
}
}
\section{Side Effects}{

}
\details{

}
\section{References}{
Yu Z, Schaid DJ. (2007) Sequential haplotype scan methods for
association analysis. Genet Epidemiol, in print.
}
\seealso{
\code{\link{haplo.em}},
\code{\link{print.seqhap}}
}
\examples{
# load example data with response and genotypes. 
setupData(seqhap.dat)
mydata.y <- seqhap.dat[,1]
mydata.x <- seqhap.dat[,-1]
# load positions
setupData(seqhap.pos)
pos=seqhap.pos$pos
# run seqhap with default settings
myobj <- seqhap(y=mydata.y, geno=mydata.x, pos=pos)
print.seqhap(myobj)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
