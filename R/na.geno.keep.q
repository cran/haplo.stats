#$Author: schaid $
#$Date: 2004/02/26 17:31:02 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/na.geno.keep.q,v 1.3 2004/02/26 17:31:02 schaid Exp $
#$Locker:  $
#$Log: na.geno.keep.q,v $
#Revision 1.3  2004/02/26 17:31:02  schaid
#changed F to FALSE
#
#Revision 1.2  2003/12/03 15:38:36  schaid
#fixed subsetting to response & covariates to not drop dim, to retain matrix class
#
#Revision 1.1  2003/09/16 16:03:08  schaid
#Initial revision
#
na.geno.keep <- function(m) {
  # determine which item in a model.frame is the genotype matrix
  gindx <- mf.gindx(m)

  # ignore genotype matrix when determining missing values for all
  # other variables (response and other covaraites)
  miss <- apply(is.na(m[, -gindx, drop=FALSE]),1,any)

  return(m[!miss,])

}

