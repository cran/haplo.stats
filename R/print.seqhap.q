#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/print.seqhap.q,v 1.3 2007/05/25 15:38:20 sinnwell Exp $
#$Locker:  $
#$Log: print.seqhap.q,v $
#Revision 1.3  2007/05/25 15:38:20  sinnwell
#change inlist to scanned.loci
#
#Revision 1.2  2007/04/16 20:11:29  sinnwell
#add ...
#
#Revision 1.1  2007/04/06 19:31:34  sinnwell
#Initial revision
#
#$Author: sinnwell $
#$Date: 2007/05/25 15:38:20 $

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

print.seqhap <- function(x, digits=max(options()$digits-2, 5), ...)


{
  
  if(x$converge==0)
    warning("EM for haplotype probabilities failed to converge")
  else
    {
      printBanner("Single-locus Chi-square Test")
      cat("Regional permuted P-value based on single-locus test is ", 
          x$chi.p.region,"\n")
      chi.test <- data.frame(chi.stat=round(x$chi.stat,digits), 
                             perm.point.p=x$chi.p.point, 
                             asym.point.p=round(1-pchisq(x$chi.stat,1),digits))
      row.names(chi.test) <- c(x$locus.label)
      print(chi.test)
      cat("\n\n")
      
      printBanner("Sequential Scan")     
      cat("Loci Combined in Sequential Analysis\n")
      scanned.loci <- x$scanned.loci
      for(i in 1:dim(scanned.loci)[1]){
        cat('seq-')
        cat(x$locus.label[i], scanned.loci[i,scanned.loci[i,]!=0],sep=" ","\n")}
      cat("\n\n")

      printBanner("Sequential Haplotype Test")
      cat("Regional permuted P-value based on sequential haplotype test is ", 
          x$hap.p.region,"\n")
      hap.test <- data.frame(hap.stat=round(x$hap.stat,digits),
                             df=x$hap.df,
                             perm.point.p=x$hap.p.point,
                             asym.point.p=round(1-pchisq(x$hap.stat,x$hap.df),digits))
      row.names(hap.test) <- paste(rep('seq-',length(x$locus.label)),
                                   x$locus.label,sep='')
      print(hap.test)
      cat("\n\n")
      
      printBanner("Sequential Summary Test")
      cat("Regional permuted P-value based on sequential summary test is ", 
            x$sum.p.region,"\n")
      sum.test <- data.frame(sum.stat=round(x$sum.stat,digits),
                             df=x$sum.df,
                             perm.point.p=x$sum.p.point,
                             asym.point.p=round(1-pchisq(x$sum.stat,x$sum.df),digits))
      row.names(sum.test) <- paste(rep('seq-',length(x$locus.label)),
                                   x$locus.label,sep='')
      print(sum.test)
    }

  invisible()
}

