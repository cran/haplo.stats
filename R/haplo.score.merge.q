#$Author: sinnwell $
#$Date: 2004/03/02 18:50:41 $
#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/haplo.score.merge.q,v 1.7 2004/03/02 18:50:41 sinnwell Exp $
#$Locker:  $
#$Log: haplo.score.merge.q,v $
#Revision 1.7  2004/03/02 18:50:41  sinnwell
#fix T to TRUE
#
#Revision 1.6  2003/08/28 15:02:48  sinnwell
#fix n.sim to simulate
#
#Revision 1.5  2003/08/26 16:39:04  sinnwell
#change license statement
#
#Revision 1.4  2003/06/19 20:53:23  sinnwell
#change "Hap-Prob" label to "Hap-Freq"
#
#Revision 1.3  2003/04/22 20:30:53  sinnwell
#use is.R() function
#
#Revision 1.2  2003/03/06 20:34:59  sinnwell
#revise to handle haplo.group list obj
#
#Revision 1.1  2003/01/17 16:29:47  sinnwell
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

haplo.score.merge <- function(score, group)
  # Created by JP Sinnwell/ DJ Schaid
  # Divisions of Biostatistics, Mayo Rochester
  # 10/2002
{

  ##  Combine the returned objects of haplo.group and haplo.score
  ##  Results will be sorted the same as haplo.score objects--by score.haplo. 
  ##  All haplotypes will be included, therefore some will not have
  ##  haplo.score info b/c ones w/ low freqs were left out (skip.haplo).
  if(!inherits(score, "haplo.score")) stop(
			"Not an object of class haplo.score!")
  if(!inherits(group, "haplo.group")) stop(
			"Not an object of class haplo.group!")

  n.loci <- group$n.loci
  df.group <- group$group.df

  # Build data frames including desired data from haplo.group and haplo.score
  group.labels <- c("Hap-Freq", names(df.group)[(n.loci+2):ncol(df.group)])

  df.score <- data.frame(score$haplotype, score$score.haplo,
                         score$score.haplo.p)

  score.labels <- c(score$locus.label, "Hap-Score", "p-val")

  #  attach on p.sim information if included in score object.
  if(score$simulate) {
     df.score <- data.frame(df.score, score$score.haplo.p.sim)
     score.labels <- c(score.labels,"sim p-val")
  }

  # merge score object onto group object by haplotypes.  Geno includes all
  # haplotypes, score will map onto those and have NA's where score left those out.
  df.merge <- merge(df.score,df.group,by=1:n.loci,all.x=TRUE,all.y=TRUE)

  # sort by score statistics
  ord <- order(df.merge[[n.loci+1]])
  df.merge <- df.merge[ord,]
  dimnames(df.merge) <- list( 1:nrow(df.merge), c(score.labels, group.labels))

  # Make a class for haplo.score.merge for print method.
  # Make it act just like a data.frame. Has one added attribute: n.loci

  if(exists("is.R") && is.function(is.R) && is.R()) {
    class(df.merge) <- c("haplo.score.merge","data.frame")
  } else {
    setOldClass(c("haplo.score.merge","data.frame"))
    oldClass(df.merge) <- "haplo.score.merge"
  }

  return(df.merge)

}
