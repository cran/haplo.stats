\name{haplo.glm.control}
\alias{haplo.glm.control}
\title{
  Create list of control parameters for haplo.glm
}
\description{
Create a list of control pararameters for haplo.glm. If no parameters
are passed to this function, then all default values are used.
}
\usage{
haplo.glm.control(haplo.effect="add", haplo.base=NULL,haplo.freq.min=0.001, 
                  sum.rare.min=0.001, haplo.min.info=0.001, 
                  keep.rare.haplo=TRUE, glm.c=glm.control(maxit=500), 
                                                  em.c=haplo.em.control())
}
\arguments{
\item{haplo.effect }{
  the "effect" of a haplotypes, which determines the covariate (x)
  coding of haplotypes. Valid options are "additive" (causing x = 0,
  1, or 2, the count of a  particular haplotype), "dominant" (causing
 x = 1 if heterozygous or homozygous carrier of a particular
  haplotype; x = 0 otherwise), and "recessive" (causing x = 1 if
  homozygous for a particular haplotype; x = 0 otherwise).
}
\item{haplo.base }{
  the index for the haplotype to be used as the base-line for
  regression. By default, haplo.base=NULL, so that the most frequent
  haplotype is chosen as the base-line.
}
\item{haplo.freq.min }{
  the minimum haplotype frequency for a haplotype to be included in
  the regression model as its own effect. The haplotype frequency is
  based on the EM algorithm that estimates haplotype frequencies
  independent of trait.
}
\item{sum.rare.min }{
  the sum of the "rare" haplotype frequencies must be larger than
  sum.rare.min in order for the pool of rare haplotypes to be included
  in the regression model as a separate term. If this condition
  is not met, then the rare haplotypes are pooled with the base-line
  haplotype (see keep.rare.haplo below).
}
\item{haplo.min.info }{
  the minimum haplotype frequency for determining the contribution of 
  a haplotype to the observed
  information matrix. Haplotypes with less frequency are
  dropped from the observed information matrix. The haplotype
  frequency is that from the final EM that iteratively updates 
  haplotype frequencies and regression coefficients. 
}
\item{keep.rare.haplo }{
  TRUE/FALSE to determine if the pool of rare haplotype should be kept as a
  separate term in the regression model (when keep.rare.haplo=TRUE), or
  pooled with the base-line haplotype (when keep.rare.haplo=FALSE).
}
\item{glm.c }{
  list of control parameters for the usual glm.control (see
glm.control).
}
\item{em.c }{
  list of control parameters for the EM algorithm to estimate
  haplotype frequencies, independent of trait (see haplo.em.control). 
}
}
\value{
  the list of above components
}
\seealso{
  glm.control, haplo.em.control
}
\examples{
# using the data set up in the example for haplo.glm,
# the control function is used in haplo.glm as follows
#  > fit <- haplo.glm(y ~ male + geno, family = gaussian,  
#  >          na.action="na.geno.keep",
#  >          data=my.data, locus.label=locus.label,
#  >          control = haplo.glm.control(haplo.freq.min =
#  >          0.02,em.c=haplo.em.control(n.try=1)))
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
