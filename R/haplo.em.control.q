#$Author: schaid $
#$Date: 2003/11/17 23:28:10 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.em.control.q,v 1.3 2003/11/17 23:28:10 schaid Exp $
#$Locker:  $
#$Log: haplo.em.control.q,v $
#Revision 1.3  2003/11/17 23:28:10  schaid
#made compatible with R
#
#Revision 1.2  2003/08/26 22:09:38  sinnwell
#added GPL License
#
#Revision 1.1  2003/08/26 21:02:09  schaid
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
haplo.em.control <- function(loci.insert.order=NULL, insert.batch.size = 6,
                              min.posterior=0.0000001, tol=0.00001, max.iter=500, random.start=0,
                              n.try = 10, iseed=NULL, max.haps.limit = 2e6, verbose=0){

  if(min.posterior < 0 | min.posterior > .9) {
    warning("The value of min.posterior is out of range, the devault value of 0.0001 is used instead")
    min.posterior <-  0.0001
  }

  if(tol < 0 | tol > .5) {
    warning("The value of tol is out of range, the default value of 0.00001 is used instead")
    tol <- 0.00001
  }


  if(max.iter < 0 | max.iter > 5000) { 
    warning("The value of max.iter is out of range, the default value of 500 is used instead")
    max.iter <- 500
  }

  if(random.start !=0  & random.start != 1){
    warning("The value of random.start is not valid, the default value of 0 is used instead")
    random.start <- 0
  }

  if(n.try < 1){
    warning("The value of n.try is not valid, the default value of 10 is used instead")
    n.try <- 10
  }

  if(max.haps.limit < 1000 | max.haps.limit > 1e12){
    warning("The value of max.haps.limit is not valid, the default value of 2e6 is used instead")
      max.haps.limit = 2e6
  }

  return(list(loci.insert.order=loci.insert.order, 
              insert.batch.size = insert.batch.size,
              min.posterior=min.posterior, 
              tol=tol, 
              max.iter=max.iter, 
              random.start=random.start,
              n.try=n.try,
              iseed=iseed,
              max.haps.limit = max.haps.limit, 
              verbose=verbose))
}

