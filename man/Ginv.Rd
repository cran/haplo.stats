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
Ginv(x, eps=1e-6)
}
\arguments{

\item{x}{
A matrix.
}
\item{eps}{
minimum cutoff for singular values in svd of x
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
The svd function uses the LAPACK standard library to compute the 
singular values of the input matrix, and the rank of the matrix is 
determined by the number of singular values  that  are  at  least  as
large   as  max(svd)*eps,  where  eps  is  a  small  value.  
For S-PLUS, the Matrix library is required (Ginv loads Matrix if not already
done so).
}
\section{References}{
Press  WH,  Teukolsky  SA,  Vetterling  WT,  Flannery  BP.
Numerical  recipes  in C. The art of scientific computing.
2nd ed. Cambridge University Press,  Cambridge.1992.  page
61.

Anderson, E., et al. (1994). LAPACK User's Guide, 
2nd edition, SIAM, Philadelphia.
}
\seealso{
svd
}
\examples{
# for matrix x, extract the generalized inverse and 
# rank of x as follows
    x <- matrix(c(1,2,1,2,3,2),ncol=3)
    save <- Ginv(x)
    ginv.x <- save$Ginv
    rank.x <- save$rank
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 37351.
