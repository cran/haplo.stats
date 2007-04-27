\name{sr.class}
\alias{sr.class}
\alias{sr.class<-}
\title{
  Return or set the class attribute of an object, compatible for both
R and S-PLUS
}
\description{
Return or set the class attribute of an object, compatible for both
R and S-PLUS.  See details.
}
\usage{
sr.class(x)
sr.class(x) <- value
}
\arguments{
\item{x}{
Any object
}
\item{value}{
character vector with the class(es) to assign to x
}
}
\value{
character vector giving the S3 class(es) of x
}
\section{Side Effects}{

}
\details{
Assigns the S3 class inheritence vector for S-PLUS objects by
oldClass() and R objects by class()
}
\section{References}{

}
\seealso{

}
\examples{
# create a simple object
myobj <- "This is an object"
sr.class(myobj) <- c("temp.class", "character")

# observe the characteristics of myobj
myobj
sr.class(myobj)  # returns: "temp.class" "character"
inherits(myobj, "character") # returns TRUE
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
