#$Author: sinnwell $
#
#$Date: 2005/03/31 15:18:08 $
#
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/Ginv.q,v 1.7 2005/03/31 15:18:08 sinnwell Exp $
#
#$Id: Ginv.q,v 1.7 2005/03/31 15:18:08 sinnwell Exp $
#
#$Locker:  $
#
#$Log: Ginv.q,v $
#Revision 1.7  2005/03/31 15:18:08  sinnwell
#re-new Matrix class Ginv, problem was in haplo.score
#
#
#Revision 1.6  2005/03/30 16:35:02  sinnwell
#revert to version 1.4, errors with Matrix class
#
#Revision 1.5  2005/03/17 22:35:02  sinnwell
#use LAPACK version of svd Splus: svd.Matrix; R it is svd(LINPACK=FALSE)
#
#Revision 1.4  2003/08/26 16:39:04  sinnwell
#change license statement
#
#Revision 1.3  2003/03/26 17:12:10  sinnwell
#fix assign operator from '_' to <-
#
#Revision 1.2  2003/03/06 23:02:34  sinnwell
#add license text
#
#Revision 1.1  2002/09/09 19:53:18  sinnwell
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
# fax:   507-284-9542
# email: schaid@mayo.edu
#

Ginv<-function(x){
  if(exists("is.R") && is.function(is.R) && is.R()) {
    savesvd <- svd(x, LINPACK=FALSE)
    U.svd<-savesvd$u
    V.svd<-savesvd$v
    d.svd<-savesvd$d

  } else{
    library(Matrix)
    savesvd<-svd.Matrix(x)
    U.svd<-savesvd$vectors$left
    V.svd<-savesvd$vectors$right
    d.svd<-savesvd$values
  }
  
  eps<- 1e-6
  maxd<-max(d.svd)                                  
  w<-ifelse((d.svd/maxd) < eps, rep(0,length(d.svd)), 1/d.svd)
  df<-sum(d.svd/maxd >= eps)

  Ginv <- V.svd %*% diag(w) %*% t(U.svd)

  list(Ginv=Ginv,rank=df)
}
