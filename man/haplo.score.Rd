\name{haplo.score}
\alias{haplo.score}
\title{
  Score Statistics for Association of Traits with Haplotypes
}
\description{
Compute score statistics to evaluate the association of a trait
with haplotypes, when linkage phase is unknown and diploid marker 
phenotypes are observed among unrelated subjects. For now, only 
autosomal loci are considered.  
}
\usage{
haplo.score(y, geno, trait.type="gaussian", offset = NA,
            x.adj = NA, skip.haplo=.005, locus.label=NA,
            miss.val=c(0,NA), simulate=FALSE, sim.control=score.sim.control(),
            em.control=haplo.em.control())
}
\arguments{
\item{y}{
Vector of trait values. For trait.type = "binomial", y must 
have values of 1 for event, 0 for no event.
}
\item{geno}{
Matrix of alleles, such that each locus has a pair of
adjacent columns of alleles, and the order of columns
corresponds to the order of loci on a chromosome.  If
there are K loci, then ncol(geno) = 2*K. Rows represent
alleles for each subject.
}
\item{trait.type }{
Character string defining type of trait, with values of 
"gaussian", "binomial", "poisson", "ordinal".
}
\item{offset}{
Vector of offset when trait.type = "poisson"
}
\item{x.adj }{
Matrix of non-genetic covariates used to adjust the score 
statistics.  Note that intercept should not be included, 
as it will be added in this function.
}
\item{skip.haplo }{
Skip score statistics for haplotypes with frequencies < skip.haplo
}
\item{locus.label }{
Vector of labels for loci, of length K (see definition of geno matrix)
}
\item{miss.val }{
Vector of codes for missing values of alleles
}
\item{simulate}{
Logical:  if [F]alse, no empirical p-values are computed; if [T]rue,
simulations are performed. Specific simulation parameters
can be controlled in the sim.control parameter list.  
}
\item{sim.control }{
A list of control parameters to determine how simulations are performed for simulated
p-values.  The list is created by the function
score.sim.control and the default values of this function can be
changed as desired.  See score.sim.control for details.
}
\item{em.control }{
A list of control parameters to determine how to perform the EM algorithm for
estimating haplotype frequencies when phase is unknown.  The list is
created by the function haplo.em.control - see this function for more details.
}
}
\value{
List with the following components:

\item{score.global}{
Global statistic to test association of trait with 
haplotypes that have frequencies >= skip.haplo.
}
\item{df}{
Degrees of freedom for score.global.
}
\item{score.global.p}{
P-value of score.global based on chi-square distribution, 
with degrees of freedom equal to df.
}
\item{score.global.p.sim}{
P-value of score.global based on simulations
(set equal to NA when simulate=F).
}
\item{score.haplo}{
Vector of score statistics for individual
haplotypes that have frequencies >= skip.haplo.
}
\item{score.haplo.p}{
Vector of p-values for score.haplo, based on a
chi-square distribution with 1 df.
}
\item{score.haplo.p.sim}{
Vector of p-values for score.haplo, based on
simulations (set equal to NA when simulate=F).
}
\item{score.max.p.sim}{
P-value  of  maximum  score.haplo, based on
simulations  (set equal to NA when simulate=F).
}
\item{haplotype}{
Matrix of hapoltypes  analyzed.  The ith row of
haplotype corresponds to the ith item of score.haplo,
score.haplo.p, and score.haplo.p.sim.
}
\item{hap.prob}{
Vector of haplotype probabilies, corresponding to the
haplotypes in the matrix haplotype.
}
\item{locus.label}{
Vector of labels for loci, of length  K  (same  as
input argument).
}
\item{simulate}{
Same as function input  parameter.  If [T]rue, simulation results are
included in the haplo.score object.
}
\item{n.val.global}{
Vector containing the number of valid simulations used in the
global score statistic simulation. The number of valid simulations can
be less than the number of simulations requested (by sim.control)
if simulated data sets produce unstable variances of the score
statistics.  
}
\item{n.val.haplo}{
Vector containing the number of valid simulations used in the p-value
simulations for maximum-score statistic and scores for the individual
haplotypes.
}
}
\section{Side Effects}{

}
\details{
Compute the maximum likelihood estimates of the haplotype frequencies
and the posterior probabilities of the pairs of haplotypes for each
subject using an EM algorithm. The algorithm begins with haplotypes
from a subset of the loci and progressively discards those with low
frequency before inserting more loci.  The process is repeated until
haplotypes for all loci are established.  The posterior probabilities
are used to compute the score statistics for the association of
(ambiguous) haplotypes with traits. The glm function is used to
compute residuals of the regression of the trait on the non-genetic
covariates.  
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
"Score tests for association of traits with haplotypes when
linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}
\seealso{
\code{\link{haplo.em}},
\code{\link{plot.haplo.score}},
\code{\link{print.haplo.score}},
\code{\link{haplo.em.control}},
\code{\link{score.sim.control}}
}
\examples{
  # establish all hla.demo data, remove genotypes with missing alleles
  # so haplo.score runs faster

  setupData(hla.demo)
  geno <- as.matrix(hla.demo[,c(17,18,21:24)])
  keep <- !apply(is.na(geno) | geno==0, 1, any)
  hla.demo <- hla.demo[keep,]
  geno <- geno[keep,]
  attach(hla.demo)
  label <- c("DQB","DRB","B")
 
# For quantitative, normally distributed trait:

  score.gaus <- haplo.score(resp, geno, locus.label=label, 
                            trait.type = "gaussian")
  print(score.gaus)

# For ordinal trait:
  y.ord <- as.numeric(resp.cat)
  score.ord <- haplo.score(y.ord, geno, locus.label=label,
                           trait.type="ordinal")
  print(score.ord)

# For a  binary trait and simulations,
# limit simulations to 500 in score.sim.control, default is 20000
  y.bin <-ifelse(y.ord==1,1,0)
  score.bin.sim <- haplo.score(y.bin, geno, trait.type = "binomial",
                     locus.label=label, simulate=TRUE, sim.control=
                     score.sim.control(min.sim=200,max.sim=500))

  print(score.bin.sim)

# For a binary trait, adjusted for sex and age:

  x <- cbind(male, age)
  score.bin.adj <- haplo.score(y.bin, geno, trait.type = "binomial", 
                               locus.label=label, x.adj=x)
  print(score.bin.adj)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
