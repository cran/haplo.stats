% $Author: sinnwell $
% $Date: 2008/01/08 20:39:02 $
% $Header: /people/biostat3/sinnwell/Haplo/Make/RCS/print.haplo.em.Rd,v 1.3 2008/01/08 20:39:02 sinnwell Exp $ 
% $Locker:  $
% $Log: print.haplo.em.Rd,v $
% Revision 1.3  2008/01/08 20:39:02  sinnwell
% fix comments
%
% Revision 1.2  2008/01/08 15:54:14  sinnwell
% minor change to desc
%
% Revision 1.1  2007/11/07 21:27:38  sinnwell
% Initial revision
%
% Revision 1.2  2004/02/27 15:29:37  sinnwell
% move haplo.em from keyword to 'see'
%
% Revision 1.1  2004/02/26 22:36:48  sinnwell
% Initial revision

\name{print.haplo.em}
\alias{print.haplo.em}
\title{
  Print contents of a haplo.em object
}
\description{
  Print a data frame with haplotypes and their frequencies.  Likelihood
  information is also printed.
}
\usage{
\method{print}{haplo.em}(x, digits=max(options()$digits-2, 5), nlines=NULL, ...)
}
\arguments{

\item{x}{
A haplo.em object
}
\item{digits}{
number of significant digits to print for numeric values
}
\item{nlines}{
To shorten output, print the first 1:nlines rows of the large data frame.  
}
\item{\dots }{
optional arguments for print
}
}
\value{
Nothing is returned
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
