%  $Author: sinnwell $ 
%  $Date: 2004/03/12 17:56:54 $
%  $Header: /people/biostat3/sinnwell/Rdir/Make/RCS/glm.fit.nowarn.Rd,v 1.3 2004/03/12 17:56:54 sinnwell Exp $
%  $Locker:  $ 

\name{glm.fit.nowarn}
\alias{glm.fit.nowarn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Modified from glm.fit function to not warn users for binomial
  non-integer weights. }
\description{
 An internal function for the haplo.stats library
}
\usage{
glm.fit.nowarn(x, y, weights = rep(1, nobs), start = NULL,
             etastart = NULL, mustart = NULL, offset = rep(0,nobs),
             family=gaussian(), control=glm.control(), intercept=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  \code{x} }
  \item{y}{  \code{y} }
  \item{weights}{ \code{weights} }
  \item{start}{  \code{start} }
  \item{etastart}{  \code{etastart} }
  \item{mustart}{  \code{mustart} }
  \item{offset}{  \code{offset} }
  \item{family}{  \code{family} }
  \item{control}{ \code{control} }
  \item{intercept}{ \code{intercept} }
}
\details{
 
}
\value{ }
\references{  }
\author{ Sinnwell JP }
\note{ }

% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{ \code{\link{haplo.glm}} }
\examples{ }
\keyword{ }% at least one, from doc/KEYWORDS

