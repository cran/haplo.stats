\name{score.sim.control}
\alias{score.sim.control}
\title{
  Create the list of control parameters for simulations in haplo.score
}
\description{
In the call to haplo.score, the sim.control
parameter is a list of parameters that control the simulations. This
list is created by this function, score.sim.control, making it easy to
change the default values.
}
\usage{
score.sim.control(p.threshold=0.25, min.sim=1000, max.sim=20000.,verbose=FALSE)
}
\arguments{

\item{p.threshold }{
A paremeter used to determine p-value precision from Besag and
Clifford (1991).  For a p-value calculated after min.sim simulations,
continue doing simulations until the p-value's sample standard error is
less than p.threshold * p-value.  The dafault value for p.threshold = 1/4
corresponds approximately to having a two-sided 95\% confidence
interval for the p-value with a width as wide as the p-value itself.
Therefore, simulations are more precise for smaller p-values.
Additionally, since simulations are stopped as soon as this criteria
is met, p-values may be biased high.  
}
\item{min.sim }{
The minimum number of simulations to run.
}
\item{max.sim }{
The upper limit of simulations allowed.  When the number of simulations reaches
max.sim, p-values are approximated based on simulation results at that
time.
}
\item{verbose}{
Logical, if (T)rue, print updates from every simulation to the
screen.  If (F)alse, do not print these details.
}
}
\value{
A list of the control parameters:

\item{p.threshold}{
As described above
}
\item{min.sim}{
As described above.  If run-time is an issue, a
lower minimum (e.g. 500) may be useful.
}
\item{max.sim}{
As described above
}
\item{verbose}{
As described above
}
}
\section{Side Effects}{

}
\details{
In simulations for haplo.score, employ the simulation p-value
precision criteria of Besag and Clifford (1991).  The criteria
ensures both the global and the maximum score statistic simulated
p-values be precise for small p-values.  First, perform min.sim
simulations to guarantee sufficient precision for the score
statistics on individual haplotypes.  Then continue simulations as
needed until simulated p-values for both the global and max score
statistics meet precision requirements set by p.threshold.
}
\section{References}{
Besag, J and Clifford, P.  "Sequential Monte Carlo p-values."
Biometrika.  78, no. 2 (1991): 301-304.
}
\seealso{
\code{\link{haplo.score}}
}
\examples{
# it would be used in haplo.score as appears below
#
# score.sim.500 <- haplo.score(y, geno, trait.type="gaussian", simulate=T, 
#                sim.control=score.sim.control(min.sim=500, max.sim=2000)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
