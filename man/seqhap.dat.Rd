%  $Author: sinnwell $ 
%  $Date: 2007/04/17 22:03:56 $ 
%  $Header: /people/biostat3/sinnwell/Haplo/Make/RCS/seqhap.dat.Rd,v 1.1 2007/04/17 22:03:56 sinnwell Exp $
%  $Locker:  $ 


\name{seqhap.dat}
\alias{seqhap.dat}
\alias{seqhap.pos}
\title{
  Simulated data for seqhap examples
}
\description{
 Simulated data set for the demonstration of seqhap functionality.
 Contains one column for disease status and columns representing 10 SNP
 loci with a known association.  seqhap.pos contains a column for
 chromosome position, as required by seqhap.
}
\usage{
data(seqhap.dat)
data(seqhap.pos)
}
\format{
\item{ }{
seqhap.dat contains disease status and 20 columns for a 10-locus
genotype matrix of SNPs.  seqhap.pos is a 1-column data.frame containing
the chromosome positions of the 10 loci.  
}
}
\source{
  
}
\section{References}{
Yu Z, Schaid DJ (2007) Sequantial haplotype scan methods for
association analysis.  To appear in {\em Gen Epi}.
}

\keyword{datasets}
% docclass is function