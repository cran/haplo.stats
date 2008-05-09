% $Author: sinnwell $ 
% $Date: 2008/01/08 22:23:35 $ 
% $Header: /people/biostat3/sinnwell/Haplo/Make/RCS/summary.haplo.em.Rd,v 1.2 2008/01/08 22:23:35 sinnwell Exp $ -->
% $Locker:  $ 
% $Log: summary.haplo.em.Rd,v $
% Revision 1.2  2008/01/08 22:23:35  sinnwell
% change comments
%
% Revision 1.1  2007/11/07 21:33:46  sinnwell
% Initial revision
%
% Revision 1.2  2004/04/06 20:46:27  sinnwell
% add nlines
%
% Revision 1.1  2004/02/26 22:36:18  sinnwell
% Initial revision

\name{summary.haplo.em}
\alias{summary.haplo.em}
\title{
  Summarize contents of a haplo.em object
}
\description{
Display haplotype pairs and their posterior probabilities by subject.
Also display a table with number of max haplotype pairs for a
subject versus how many were kept (max vs. used).  
}
\usage{
\method{summary}{haplo.em}(object, show.haplo=FALSE, digits=max(options()$digits-2, 5), nlines=NULL, ...)
}
\arguments{
\item{object}{
A haplo.em object
}
\item{show.haplo }{
Logical. If TRUE, show the alleles of the haplotype pairs, otherwise show
only the recoded values.
}
\item{digits}{
number of significant digits to be printed for numeric values
}
\item{nlines}{
To shorten output, print the first 1:nlines rows of the large data frame.  
}
\item{\dots }{
Optional arguments for the summary method
}
}
\value{
}
\section{Side Effects}{

}
\details{

}
\section{References}{

}
\seealso{
\code{\link{haplo.em}}
}
\examples{

}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
