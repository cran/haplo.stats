\name{plot.seqhap}
\alias{plot.seqhap}
\title{
  Plot a seqhap object
}
\description{
Method to plot an object of class seqhap. The p-values at each
locus are based on sequentially combined loci, and they
are plotted to visualize the p-values when scanning each locus using
seqhap methods. Plots -log10(p-value) on the y-axis vs. the
loci over which it was computed on the x-axis.
}
\usage{
plot.seqhap(x, pval="hap", single=TRUE, ...)
}
\arguments{
\item{x}{
The object returned from seqhap
}
\item{pval}{
Character string for the choice of p-value to plot. Options are:
"hap" (sequential haplotype asymptotic p-value), "hap.sim" (sequential
haplotype simulated p-value), "sum" (sequential summary asymptotic
p-value), and "sum.sim" (sequential summary simulated p-value).  
}
\item{single}{
Logical, indicating whether to plot p-values for single-locus
association tests. If TRUE, the pointwise p-values from the
single-locus will be plotted using a dotted line.
}
\item{\dots }{
Dynamic parameter for the values of additional parameters for 
the plot method.  Accept values that will apply to plot() and axis().
Recommended values to make locus labels vertical on the x-axis:
for R:  las=2,  cex.axis=1.2
for S+: srt=90, cex.axis=1.2, adj=1
}
}
\value{
Nothing is returned.
}
\section{Side Effects}{

}
\details{
The x-axis has tick marks for all loci.  The y-axis is the -log10() of
the selected p-value.  For the sequential result for each locus, a
horizontal line at the height of -log10(p-value) is drawn across the
loci combined. The start locus is indicated by a filled triangle and
other loci combined with the start locus are indicated by an asterisk
or circle.

If the simulated p-value for a set of loci is zero, it is set to .5 /
n.sim.
}
\section{References}{
Yu Z, Schaid DJ. (2007) Sequential haplotype scan methods for
association analysis. Genet Epidemiol, in print.
}
\seealso{
\code{\link{seqhap}},
\code{\link{print.seqhap}}
}
\examples{
setupData(seqhap.dat)
mydata.y <- seqhap.dat[,1]
mydata.x <- seqhap.dat[,-1]
setupData(seqhap.pos)
myobj <- seqhap(y=mydata.y, geno=mydata.x, pos=seqhap.pos$pos, n.sim=100)
plot(myobj)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
