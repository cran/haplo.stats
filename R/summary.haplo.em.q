#$Author: sinnwell $
#$Date: 2004/02/26 23:10:39 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/summary.haplo.em.q,v 1.6 2004/02/26 23:10:39 sinnwell Exp $
#$Locker:  $
#$Log: summary.haplo.em.q,v $
#Revision 1.6  2004/02/26 23:10:39  sinnwell
#print.banner to printBanner
#
#Revision 1.5  2004/02/16 19:41:01  sinnwell
#add '...' to make method for haplo.em object
#
#Revision 1.4  2003/12/08 20:24:09  sinnwell
# changed T,F to TRUE,FALSE
#
#Revision 1.3  2003/10/15 15:41:18  schaid
#added show.haplo option for showing haplotypes instead of their codes
#
#Revision 1.2  2003/08/26 22:09:38  sinnwell
#added GPL License
#
#Revision 1.1  2003/08/26 21:00:25  schaid
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
summary.haplo.em <- function(object, show.haplo=FALSE, ...){

  
  printBanner("Subjects: Haplotype Codes and Posterior Probabilities")

  if(show.haplo){
     hap1 <- object$haplotype[object$hap1code,]
     hap2 <- object$haplotype[object$hap2code,]
     df <- data.frame(subj.id=object$subj.id, hap1=hap1, hap2=hap2,
                    posterior=round(object$post,5))
  }  else{
     df <- data.frame(subj.id=object$subj.id, hap1code=object$hap1code, hap2code=object$hap2code,
                    posterior=round(object$post,5))
  }

  print(df)
 
  printBanner("Number of haplotype pairs: max vs used")
 

  x <- object$max.pairs
  if(length(object$rows.rem) > 0 ){
    x <- x[-object$rows.rem]
  } 
  print(table(x, object$nreps))


  invisible()
}
