#$Author: sinnwell $
#$Date: 2008/05/08 17:38:26 $
#$Header: /people/biostat3/sinnwell/Haplo/Make/RCS/plot.seqhap.q,v 1.3 2008/05/08 17:38:26 sinnwell Exp $
#$Locker:  $
#$Log: plot.seqhap.q,v $
#Revision 1.3  2008/05/08 17:38:26  sinnwell
#allow ylim, set minimum p as 1e-10, issue warning if pval not valid
#
#Revision 1.2  2007/05/25 18:17:15  sinnwell
#rm par options cex.axis and las/srt, reserve for users
#
#Revision 1.1  2007/05/23 17:41:36  sinnwell
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

plot.seqhap <- function(x, pval="hap", single=TRUE, ...) {
  # plot method of seqhap object
  #  1. x-axis: give locus labels in order, in "pos" chrom positions
  #  2. y-axis: plot the -log10 of the p-val chosen
  #  3. the pval parameter selects which ("hap","hap.sim","sum","sum.sim") to plot

  if(!inherits(x, "seqhap"))
	    stop("Not a seqhap object")

  ## Set p.small as the smallest "allowable" p-value.
  ## Any p-value smaller will be set to log10(p.small)
  small.p <- 1e-10
  
  ## if ylim is given in ..., pull it off and use it
  dots <- as.list(substitute(list(...)))[-1]
  ylim.indx <- match('ylim', names(dots))
  ylim.user <- if(!is.na(ylim.indx)) eval(dots[ylim.indx]$ylim) else NULL

  single.p <- 1-pchisq(x$chi.stat,1)

  if(is.na(match(pval, c("hap", "sum", "hap.sim", "sum.sim")))) {
    warning("Invalid value for pval, set to default. \n")
    pval <- "hap"
  }
  
  switch(pval,
         "hap"={seqp <- 1-pchisq(x$hap.stat,x$hap.df)
                ylabel <- '-log10(hap.pval)' 
                 },
         "sum"={seqp <- 1-pchisq(x$sum.stat,x$sum.df)
                ylabel <- '-log10(sum.pval)'
                 },
         "hap.sim"={seqp <- x$hap.p.point
                    ylabel <- '-log10(hap.sim.pval)'
                 },
         "sum.sim"={seqp <- x$sum.p.point
                    ylabel <- '-log10(sum.sim.pval)'
                 })

  if(any(single.p < small.p)) {
    cat("One or more single-marker p-values too small to plot, set to ", small.p, "\n")
    single.p <- ifelse(single.p < small.p, small.p, single.p)
  }
  if(any(seqp < small.p)) {
    cat("One or more multi-marker p-values too small to plot, set to ", small.p, "\n")
    seqp <- ifelse(seqp < small.p, small.p, seqp)
  }

  log.single.p <- if(single) -log10(single.p) else NULL
  log.seqp <- -log10(seqp)
  
  ylim <- if(is.null(ylim.user)) c(0, max(3, log.single.p, log.seqp)) else ylim.user
  
  plot(x$pos, log.single.p,ylim=ylim, type="n",xaxt="n",
       xlab='', ylab=ylabel) #, ...)

  for(i in 1:length(x$pos))
    {
      # plot a filled triangle for locus i at heigh of -log10p[i]
      points(x$pos[x$scanned.loci[i,1]], log.seqp[i], pch=17)
      if(sum(x$scanned.loci[i,]>0)>1)
        {
          scanned.loci.i <- sort(x$scanned.loci[i,])
          scanned.loci.i <- scanned.loci.i[scanned.loci.i>0]
          
          #plot a line for each set of combined loci at height of -log10p[i]
          lines(c(x$pos[min(scanned.loci.i)],x$pos[max(scanned.loci.i)]),rep(log.seqp[i],2))
          
          #plot an asterisk/circle on each locus combined at height of -log10p[i]
          points(x$pos[setdiff(scanned.loci.i, x$scanned.loci[i,1])],rep(log.seqp[i],length(scanned.loci.i)-1), ...)
         }
    }

  if(single) lines(x$pos,log.single.p,lty=2)
  if(is.R()) {  #mgp for title change
    axis(1,at=x$pos, labels=x$locus.label, ...) # adj=0, las=2, cex.axis=.8
  } else {
    axis(1,at=x$pos, labels=x$locus.label, ...)  #cex.axis=.8, srt=90, adj=0
  }

  invisible()

}

