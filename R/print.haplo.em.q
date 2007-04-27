#$Author: sinnwell $
#$Date: 2004/04/06 20:39:39 $
#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/print.haplo.em.q,v 1.4 2004/04/06 20:39:39 sinnwell Exp $
#$Locker:  $
#$Log: print.haplo.em.q,v $
#Revision 1.4  2004/04/06 20:39:39  sinnwell
#use nlines to limit printouts in vignettes
#
#Revision 1.3  2004/02/26 23:04:28  sinnwell
#print.banner to printBanner
#
#Revision 1.2  2003/08/26 22:09:38  sinnwell
#added GPL License
#
#Revision 1.1  2003/08/26 20:59:40  schaid
#Initial revision
#
# License: 
# 
# Copyright 2003 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; either 
# version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
# more details.
# 
# You should have received a copy of the GNU General Public License along with this 
# program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA
# 
# For other licensing arrangements, please contact Daniel J. Schaid.
# 
# Daniel J. Schaid, Ph.D.
# Division of Biostatistics
# Harwick Building – Room 775
# Mayo Clinic
# 200 First St., SW
# Rochester, MN 55905
# 
# phone: 507-284-0639
# fax:      507-284-9542
# email: schaid@mayo.edu
# 
print.haplo.em <- function(x, nlines=NULL, ...){

  printBanner("Haplotypes")
  df <- data.frame(x$haplotype,round(x$hap.prob,5))
  names(df) <- c(x$locus.label, "hap.freq")
  if(is.null(nlines)) print(df)
  else print(df[1:nlines,])
  invisible()

  if(x$converge==0)
     warning("EM failed to converge")

  printBanner("Details")

  pval <- NA
  if(x$df.lr > 0) pval = 1-pchisq(x$lr, x$df.lr)
  
  cat("lnlike = ",round(x$lnlike,6),"\n")
  cat("lr stat for no LD = ",round(x$lr,6),", df = ",x$df.lr,", p-val = ",round(pval,6),"\n")
  n.rem <- length(x$rows.rem)
  if(n.rem > 0) {
     cat("\nRow number of subjects removed because max number of pairs of haplotypes > enum.limit\n")
     tbl <- data.frame(row.num=x$rows.rem, max.pairs=x$max.pairs[x$rows.rem])
    
   }

  invisible()
}
