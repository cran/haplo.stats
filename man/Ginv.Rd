\name{Ginv}
\alias{Ginv}
\title{
  Compute Generalized Inverse of Input Matrix
}
\description{
Singular value decomposition (svd) is used to compute a
generalized inverse of input matrix.
}
\usage{
Ginv(x)
}
\arguments{

\item{x}{
A matrix.
}
}
\value{
List with components:

\item{Ginv}{
Generalized inverse of x.
}
\item{rank}{
Rank of matrix x.
}
}
\section{Side Effects}{

}
\details{
The function svd is used to compute the singular values of
the input matrix, and the rank of the matrix is determined
by the number of singular values  that  are  at  least  as
large   as  max(svd)*eps,  where  eps  is  a  small  value
(currently eps = .000001).
}
\section{References}{
Press  WH,  Teukolsky  SA,  Vetterling  WT,  Flannery  BP.
Numerical  recipes  in C. The art of scientific computing.
2nd ed. Cambridge University Press,  Cambridge.1992.  page
61.
}
\seealso{
svd
}
\examples{
# for matrix x, extract the generalized inverse and 
# rank of x as follows
#   > save <- Ginv(x)
#   > ginv.x <- save$Ginv
#   > rank.x <- save$rank
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
