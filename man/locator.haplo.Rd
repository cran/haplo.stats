\name{locator.haplo}
\alias{locator.haplo}
\title{
  Find Location from Mouse Clicks and Print Haplotypes on Plot
}
\description{
Much like the R/Splus locator function is used to find x-y coordinates
on a plot.  Find all x-y coordinates that are chosen by the user's
mouse clicks.  Then print haplotype labels at the chosen positions.
}
\usage{
locator.haplo(obj)
}
\arguments{
\item{obj}{
An object (of class haplo.score) that is returned from haplo.score.
}
}
\value{
List with the following components:

\item{x.coord}{
Vector of x-coordinates.
}
\item{y.coord}{
Vector of y-coordinates.
}
\item{hap.txt}{
Vector of character strings for haplotypes.
}
}
\details{
After plotting the results in obj, as from plot(obj), the function
locator.haplo is used to place on the plot the text strings for
haplotypes of interest. After the function call (e.g.,
locator.haplo(obj)), the user can click, with the left mouse button,
on as many points in the plot as desired. Then, clicking with the
middle mouse button will cause the haplotypes to be printed on the
plot. The format of a haplotype is "a:b:c", where a, b, and c are
alleles, and the separator ":" is used to separate alleles on a
haplotype. The algorithm chooses the closest point that the user
clicks on, and prints the haplotype either above the point (for points
on the lower-half of the plot) or below the point (for points in the
upper-half of the plot).
}
\seealso{
\code{\link{haplo.score}}
}
\examples{
# follow the pseudo-code
#  score.out <-  haplo.score(y, geno, trait.type = "gaussian")

#  plot(score.out)

#  locator.haplo(score.out)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
