\name{haplo.cc}
\alias{haplo.cc}
\title{
  Haplotype Association Analysis in a Case-Control design
}
\description{
Combine results from haplo.score, haplo.group, and haplo.glm for
case-control study designs.  Analyze the association between the 
binary (case-control) trait and the haplotypes relevant to the
unrelated individuals' genotypes.
}
\usage{
haplo.cc(y, geno, haplo.min.count=5, locus.label=NA, ci.prob=0.95, miss.val=c(0,NA),
                     simulate=FALSE, sim.control=score.sim.control(),
                     control=haplo.glm.control())
}
\arguments{
\item{y}{
Vector of trait values, must be 1 for cases and 0 for controls.
}
\item{geno}{
Matrix of alleles, such that each locus has a pair of
adjacent columns of alleles, and the order of columns
corresponds to the order of loci on a chromosome.  If
there are K loci, then ncol(geno) = 2*K. Rows represent
alleles for each subject.
}
\item{haplo.min.count }{
The minimum number of expected counts for a haplotype in the sample
based on estimated frequencies. 
}
\item{ci.prob }{
Probability level for confidence interval on the Odds Ratios of
each haplotype to span the true value. 
}
\item{locus.label }{
Vector of labels for loci, of length K (see definition of geno matrix)
}
\item{miss.val }{
Vector of codes for missing values of alleles
}
\item{simulate}{
Logical:  if [F]alse, no empirical p-values are computed; if [T]rue,
simulations are performed within haplo.score. Specific simulation parameters
can be controlled in the sim.control parameter list.  
}
\item{sim.control }{
A list of control parameters to determine how simulations are
performed for simulated p-values.  The list is created by the function
score.sim.control and the default values of this function can be
changed as desired.  See score.sim.control for details.
}
\item{control}{
A list of control parameters for managing the execution of
haplo.cc. The list is created by the function haplo.glm.control, which
also manages control parameters for the execution of haplo.em.
}
}
\value{
A list including the haplo.score object (score.lst), vector of subject
counts by case and control group (group.count), haplo.glm object (fit.lst), 
confidence interval probability (ci.prob), and a data frame (cc.df) with the 
following components:

\item{haplotypes}{
The first K columns contain the haplotypes used in the analysis.
}
\item{Hap-Score}{
Score statistic for association of haplotype with the binary trait.
}
\item{p-val}{
P-value for the haplotype score statistic, based on a chi-square
distribution with 1 degree of freedom.
}
\item{sim.p.val}{
Vector of p-values for score.haplo, based on simulations in
haplo.score (omitted when simulations not performed).
P-value of score.global based on simulations
(set equal to NA when simulate=F).
}
\item{pool.hf}{
Estimated haplotype frequency for cases and controls pooled together.
}
\item{control.hf}{
Estimated haplotype frequency for control group subjects.
}
\item{case.hf}{
Estimated haplotype frequency for case group subjects.
}
\item{glm.eff}{
The haplo.glm function modeled the haplotype effects as: baseline (Base),
additive haplotype effect (Eff), or rare haplotypes pooled into a
single group (R).
}
\item{OR.lower}{
Lower limit of the Odds Ratio Confidence Interval.
}
\item{ OR }{
Odds Ratio based on haplo.glm model estimated coefficient for the haplotype.
}
\item{OR.upper}{
Upper limit of the Odds Ratio Confidence Interval.
}
}
\section{Side Effects}{

}
\details{
All function calls within haplo.cc are for the analysis of association
between haplotypes and the case-control status (binomial trait).  No 
additional covariates may be modeled with this function.  Odd Ratios
are in reference to the baseline haplotype.  Odds Ratios will change
if a different baseline is chosen using haplo.glm.control.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
"Score tests for association of traits with haplotypes when
linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}
\section{}{
Lake S, LH, Silverman E, Weiss S, Laird N, Schaid DJ.  "Estimation and
tests of haplotype-environment interaction when linkage phase is
ambiguous. Human Heredity. 55 (2003):  56-65
}
\seealso{
\code{\link{haplo.em}},
\code{\link{haplo.score}},
\code{\link{haplo.group}},
\code{\link{haplo.score.merge}},
\code{\link{haplo.glm}}
\code{\link{print.haplo.cc}}
}
\examples{
#  For a genotype matrix geno.test, case/control vector y.test
#  The function call will be like this
#  cc.test <- haplo.cc(y.test, geno.test, locus.label=locus.label, haplo.min.count=3, ci.prob=0.95)
#
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
