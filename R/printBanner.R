#$Author: sinnwell $
#
#$Date: 2004/02/26 21:35:37 $
#
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/printBanner.R,v 1.1 2004/02/26 21:35:37 sinnwell Exp $
#
#$Id: printBanner.R,v 1.1 2004/02/26 21:35:37 sinnwell Exp $
#
#$Locker:  $
#
#$Log: printBanner.R,v $
#Revision 1.1  2004/02/26 21:35:37  sinnwell
#Initial revision
#
#Revision 1.2  2003/08/26 16:41:23  sinnwell
#change License statement
#
#Revision 1.1  2003/03/06 21:03:39  sinnwell
#Initial revision
#
#Revision 1.1  2002/09/09 19:53:18  sinnwell
#Initial revision
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
printBanner <- function(str, banner.width=80, char.perline=60, border = "=")
## created by Gregory Warnes, Groton Pfizer 2002
## used in haplo.score package by DJ Schaid, JP Sinnwell
## Mayo Biostatistics 3/2003
{

  
  ## Print Border
  cat( rep(border, banner.width), sep="" )
  cat( "\n\n" )

  # Split up words into lines
  lines <- strwrap( str, width=char.perline)

  # Roughly equivalent to, but better:
  # words <- strsplit(str, ' +', extended=T)
  # cumlen <- cumsum( sapply( words, function(x) nchar(x)+1 ) )
  # lineno <- cut( cumlen, c(seq(0, max(cumlen)-1, by=60), max(cumlen) ) )
  # lines <- sapply( split( unlist(words), lineno ), function(x) paste(x, collapse=" ") )
  #

  strcenter <- function( str )
  {
    spaces <- floor( ( banner.width - nchar(str) )/2 )
    if(spaces > 0)
      retval <- paste( rep(" ", spaces), collapse="" )
    else
      retval <- ""
    retval <- paste(retval, str, sep="" )
    return( retval )
  }

  # Print each line centerd
  cat( sapply(lines, strcenter), sep="\n")

  cat("\n")

  # Print Border
  cat( rep(border, banner.width), sep="" )
  cat( "\n\n" )
}
