\name{printBanner}
\alias{printBanner}
\title{
  Print a nice banner
}
\description{

}
\usage{
printBanner(str, banner.width=options()$width, char.perline=.75*banner.width, border="=")
}
\arguments{
\item{str}{
character string - a title within the banner
}
\item{banner.width }{
width of banner, the default is set to fit current options
}
\item{char.perline }{
number of characters per line for the title, the default is 75\% of 
the banner.width parameter 
}
\item{border}{
type of character for the border
}
}
\value{
}
\section{Side Effects}{

}
\details{
This function prints a nice banner in both R and S-PLUS
}
\section{References}{

}
\seealso{
options
}
\examples{
printBanner("This is a pretty banner", banner.width=40, char.perline=30)

# the output looks like this:
# ========================================
#         This is a pretty banner
# ========================================
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
