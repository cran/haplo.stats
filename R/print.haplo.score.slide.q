#$Author: sinnwell $
#$Date: 2003/09/19 16:48:38 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/print.haplo.score.slide.q,v 1.2 2003/09/19 16:48:38 sinnwell Exp $
#$Locker:  $
#$Log: print.haplo.score.slide.q,v $
#Revision 1.2  2003/09/19 16:48:38  sinnwell
#fix digits syntax error
#
#Revision 1.1  2003/08/22 21:09:41  sinnwell
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

print.haplo.score.slide <- function(x, digits=max(options()$digits-2, 5), ...)
  # Developed by Schaid, DJ; Sinnwell, JP 2003
  # Mayo Clinic Rochester, Div of Biostatistics
{
  
# print method for haplo.score.slide object.  Only print the data frame portion
# of the object.  It has the start and end locus,
#  and the p-values of that call to haplo.score.
  
  print(as.data.frame(x$df), digits=digits, ...)
  invisible()
  
}
