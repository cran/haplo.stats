\name{haplo.em.control}
\alias{haplo.em.control}
\title{
  Create the Control Parameters for the EM Computation of Haplotype 
Probabilities, with Progressive Insertion of Loci
}
\description{
Create a list of parameters that control the EM algorithm for 
estimating haplotype frequencies, based on progressive insertion 
of loci. Non-default parameters for the EM algorithm can be
set as parameters passed to haplo.em.control.
}
\usage{
haplo.em.control(loci.insert.order=NULL, insert.batch.size = 6,
                             min.posterior = 1e-07, tol = 1e-05,
                             max.iter=5000, random.start=0, n.try = 10,
                             iseed=NULL, verbose=0)
}
\arguments{

\item{loci.insert.order }{
Numeric vector with specific order to insert the loci. If this value
is NULL, the insert oder will be in sequential order (1, 2, ..., No. Loci).
}
\item{insert.batch.size }{
Number of loci to be inserted in a single batch.
}
\item{min.posterior }{
Minimum posterior probability of a haplotype pair, conditional on
observed marker genotypes. Posteriors below this minimum value will
have their pair of haplotypes "trimmed" off the list of possible
pairs.
}
\item{tol}{
If the change in log-likelihood value between EM steps is less 
than the tolerance (tol), it has converged.
}
\item{max.iter }{
Maximum number of iterations allowed for the EM algorithm before it
stops and prints an error.  If the error is printed, double max.iter.
}
\item{random.start }{
If random.start = 0, then the inititial starting values of the
posteriors for the first  EM attempt will be based on assuming equal
posterior probabilities (conditional on genotypes). If random.start = 1,
then the initial starting values of the first EM attempt will be based on assuming a uniform distribution for
the initial posterior probabilities.
}
\item{n.try }{
Number of times to try to maximize the lnlike by the EM algorithm. The
first try uses, as initial starting values for the posteriors,
either equal values or uniform random variables, as determined by
random.start. All subsequent tries will use random uniform values as
initial starting values for the posterior probabilities.
}
\item{iseed}{
An integer or a saved copy of .Random.seed. This allows simulations to
be reproduced by using the same initial seed.
}
\item{verbose}{
Logical, if TRUE, print procedural messages to the screen.  If FALSE, 
do not print any messages. 
}
}
\value{
A list of the parameters passed to the function.
}
\section{Side Effects}{

}
\details{
The default is to use n.try = 10. If this takes too much time, it may
be worthwhile to decrease n.try. Other tips for computing haplotype
frequencies for a large number of loci, particularly if some have many
alleles, is to decrease the batch size (insert.batch.size), increase
the memory (max.haps.limit).
}
\section{References}{

}
\seealso{
\code{\link{haplo.em}},
\code{\link{haplo.score}}
}
\examples{
# This is how it is used within haplo.score
#    > score.gauss <- haplo.score(resp, geno, trait.type="gaussian", 
#    >           em.control=haplo.em.control(insert.batch.size = 2, n.try=1))
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
