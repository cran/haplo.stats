%  $Author: sinnwell $ 
%  $Date: 2011/11/10 15:29:41 $
%  $Header: /projects/genetics/cvs/cvsroot/haplo.stats/man/glm.fit.nowarn.Rd,v 1.6 2011/11/10 15:29:41 sinnwell Exp $
%  $Locker:  $ 

\name{glm.fit.nowarn}
\alias{glm.fit.nowarn}
\title{ Modified from glm.fit function to not warn users for binomial
  non-integer weights. }
\description{
 An internal function for the haplo.stats library
}
\usage{
glm.fit.nowarn(x, y, weights = rep(1, nobs), start = NULL,
             etastart = NULL, mustart = NULL, offset = rep(0,nobs),
             family=gaussian(), control=list(), intercept=TRUE)
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
\author{ Sinnwell JP }

% ~Make other sections like Warning with \section{Warning }{....} ~ 
\seealso{ \code{\link{haplo.glm}} }

