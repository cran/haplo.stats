%  $Author: sinnwell $ 
%  $Date: 2004/03/12 21:51:18 $ 
%  $Header: /people/biostat3/sinnwell/Rdir/Make/RCS/hla.demo.Rd,v 1.3 2004/03/12 21:51:18 sinnwell Exp $
%  $Locker:  $ 


\name{hla.demo}
\alias{hla.demo}
\title{
  HLA Loci and Serologic Response to Measles Vaccination.
}
\description{
  Eleven HLA-region loci genotyped for 220 subjects, phase not
  known. Contains measles vaccination response with covariate data.
}
\usage{
data(hla.demo)
}
\format{
\item{ }{
Data Frame with the following columns:
}
\item{ resp }{
Quantitative response to Measles Vaccination 
}
\item{ resp.cat }{
Category of response as low, normal, or high; based on 'resp'
}
\item{ male }{
Binary indicator of gender, 1=male, 0=female
}
\item{ age }{
Age of the subject
}
\item{ allele columns 5 - 26 }{
  Pairs of columns represent the allele pairs for each subject at the
  locus.
}
}

\source{
  Data set kindly provided by Gregory A. Poland, M.D. and the Mayo
  Clinic Vaccine Research Group for illustration only, and my not be
  used for publication.
}
\section{References}{
  Schaid DJ, Rowland CM, Tines DE, Jacobson RM,  Poland  GA.
  "Score tests for association of traits with haplotypes when
  linkage phase is ambiguous." Amer J Hum Genet. 70 (2002):  425-434.
}

\keyword{datasets}
% docclass is function
