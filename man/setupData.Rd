\name{setupData}
\alias{setupData}
\title{
  Set up an example dataset provided within the library.
}
\description{
This function defines an alias function to run exactly as data() in R
and does nothing in Splus.  R keeps a data set within the working data
frame, so we only want to load data it when calling an example.  Splus
keeps it in the background, so it is already loaded upon library(mypkg).
}
\usage{
setupData(...)
}
\arguments{

\item{\dots }{
The name of a dataset provided within the Splus/R library.
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

}
\examples{
## for a data set named my.data load it by
# setupData(my.data)

## check the names of my.data to see if it is loaded
# names(my.data)
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
