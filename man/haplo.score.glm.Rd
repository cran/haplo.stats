\name{haplo.score.glm}
\alias{haplo.score.glm}
\title{
  Compute haplotype score statistics for GLM
}
\description{
This function is used by haplo.score when analyzing traits
by a GLM score.
}
\usage{
haplo.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
}
\arguments{
\item{y}{
Vector of trait values.
}
\item{mu}{
Expected value of y.
}
\item{a}{
scale parameter
}
\item{v}{
v= b''/a for a GLM.
}
\item{x.adj }{
Matrix of non-genetic covariates used  to  adjust  the
score  statistics. Note that intercept should  be included
in this matrix.
}
\item{nreps}{
Vector for the count of haplotype  pairs  that  map  to
each subject's marker genotypes (see haplo.em).
}
\item{x.post }{
Matrix for posterior mean of x per subject.
}
\item{post}{
Vector of posterior probabilities of pairs of
haplotypes for a person, given thier marker phenotypes
(see haplo.em).
}
\item{x}{
Matrix  of  scores  for  enumerated  haplotypes  for  each
subject,  with  elements  0,  1,  2  (counts  of  specific
haplotypes).
}
None.
}
\value{
List with components:

\item{u.score}{
Vector of scores for the chosen haplotypes
}
\item{v.score}{
Covariance matrix for u.score
}
}
\section{Side Effects}{

}
\details{
Using posterior probabilities of pairs of haplotypes,  the
"design" matrix for the haplotype effects, and the GLM
residuals, compute the score vector and its variance
matrix, adjusted for the non-genetic covariates.
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
Score tests for association of traits with haplotypes when
linkage phase is ambiguous. Submitted to Amer J Hum Genet.
}
\seealso{
haplo.score
}
\examples{

}
\keyword{}
% docclass is functionhaplo
% Converted by Sd2Rd version 37351.
