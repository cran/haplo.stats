\name{plot.haplo.score.slide}
\alias{plot.haplo.score.slide}
\title{
  Plot a haplo.score.slide Object 
}
\description{
Method function to plot an object of class haplo.score.slide.  The
p-values from haplo.score.slide are for sub-haplotypes of a larger
chromosomal region, and these are plotted to visualize the change in
p-values as the sub-haplotype "slides" over a chromosome. Plot
-log10(p-value) on the y-axis vs. the loci over which it was computed
on the x-axis.
}
\usage{
plot.haplo.score.slide(x, pval="global", dist.vec=1:x$n.loci, 
                       cex=.8, srt=270, ...)
}
\arguments{
\item{x}{
The object returned from haplo.score.slide
}
\item{pval}{
Character string for the choice of p-value to plot. Options are:
"global" (the global score statistic p-value based on an asymptotic
chi-square distribution), "global.sim" (the global score statistic
simulated p-value), and "max.sim" (the simulated p-value for the
maximum score statistic).
}
\item{dist.vec }{
Numeric vector for position (i.e. in cM) of the loci along a
chromosome. Distances on x-axis will correspond to these positions. 
}
\item{cex}{
Character expansion size.
}
\item{srt}{
String rotation in degrees measured counterclockwise from
horizontal.  Applies to x-axis (locus) labels.  
}
\item{\dots }{
Dynamic parameter for the values of additional parameters for 
the plot method.  
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
The x-axis has tick marks for all loci.  The y-axis is the -log10() of
the selected p-value.  For each haplo.score result, plot a horizontal
line at the height of -log10(p-value) drawn across the loci over which
it was calculated.  Therefore a p-value of 0.001 for the first 3 loci
will plot as a horizontal line plotted at y=3 covering the first three
tick marks.  
}
\section{References}{
Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
"Score tests for association of traits with haplotypes when
linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}
\seealso{
\code{\link{haplo.score.slide}}
}
\examples{
#This example is run completely in the haplo.score.slide

#   setupData(hla.demo)
#   attach(hla.demo)
#   geno.11 <- hla.demo[,-c(1:4)]
#   label.11 <- c("DPB","DPA","DMA","DMB","TAP1","TAP2","DQB","DQA","DRB","B","A")

#For an ordinal trait, slide by 3 loci, and simulate p-values:
#   y.ord <- as.numeric(resp.cat)
#   slide.ord.sim <-  haplo.score.slide(y.ord, geno.11, trait.type = "ordinal",
#                              n.slide=3, locus.label=label.11, simulate=TRUE,
#                              sim.control=score.sim.control(min.sim=500))

#   print(slide.ord.sim)
#   plot(slide.ord.sim)
#   plot(slide.ord.sim, pval="global.sim")
#   plot(slide.ord.sim, pval="max.sim")
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
