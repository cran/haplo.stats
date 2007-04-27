#$Author: sinnwell $
#$Date: 2004/07/02 14:18:11 $
#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/zzz.haplo.stats.R,v 1.1 2004/07/02 14:18:11 sinnwell Exp $
#$Locker:  $
#$Log: zzz.haplo.stats.R,v $
#Revision 1.1  2004/07/02 14:18:11  sinnwell
#Initial revision
#
#Revision 1.5  2003/10/06 15:45:17  sinnwell
#dyn.load haplo.stats, not haplo.score
#
#Revision 1.4  2003/08/26 16:41:23  sinnwell
#change License statement
#
#Revision 1.3  2003/03/20 22:26:56  sinnwell
#remove source() command
#
#Revision 1.2  2003/03/17 16:45:41  sinnwell
#autoload data in hla.demo.q in .First.Lib()
#
#Revision 1.1  2003/03/06 23:25:55  sinnwell
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
# 
.First.lib <- function(lib, pkg) {
   library.dynam("haplo.stats", pkg, lib)
}
