#$Author: sinnwell $
#$Date: 2003/12/08 19:42:18 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.score.q,v 1.13 2003/12/08 19:42:18 sinnwell Exp $
#$Locker:  $
#$Log: haplo.score.q,v $
#Revision 1.13  2003/12/08 19:42:18  sinnwell
# changed F,T to FALSE,TRUE
#
#Revision 1.12  2003/09/17 22:56:25  schaid
#added check for var(u.score) < 0 when computing simulated haplo.score.sim
#
#Revision 1.11  2003/09/17 22:35:06  schaid
#modified how max-stat and score.haplo simulated p-values account for NA results (either NA for observed stat, or NA for simulated stats)
#
#Revision 1.10  2003/09/15 14:53:17  sinnwell
#add in na.rm in score.max.sim calculation.
#may be part of a bigger problem
#
#Revision 1.9  2003/09/11 21:27:54  sinnwell
#change simulations to target both global and max for precision
#
#Revision 1.8  2003/08/27 21:19:06  sinnwell
#*** empty log message ***
#
#Revision 1.7  2003/08/27 20:59:18  sinnwell
#add rows.rem back into return list, it was needed in haplo.group
#
#Revision 1.6  2003/08/22 18:05:24  sinnwell
#update for release of haplo.score 1.2.0 with 'PIN'
#
#Revision 1.4  2003/04/22 20:30:20  sinnwell
#include use of is.R()
#
#Revision 1.3  2003/03/06 21:48:45  sinnwell
#include license statement
#
#Revision 1.2  2003/01/17 16:57:06  sinnwell
#revision for haplo.score version 1.2
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

haplo.score <- function(y, geno, trait.type="gaussian",
                         offset = NA, x.adj = NA, skip.haplo=.005,
                         locus.label=NA, miss.val=c(0,NA),
                         simulate=FALSE, sim.control=score.sim.control(),
                         em.control = haplo.em.control())
{

# Choose trait type:
   trait.int <- charmatch(trait.type, c("gaussian", "binomial", "poisson",
                          "ordinal"))

   if(is.na(trait.int)) stop("Invalid trait type")
   if(trait.int == 0)   stop("Ambiguous trait type")

# Check dims of y and geno
  if(length(y)!=nrow(geno)) stop("Dims of y and geno are not compatible")
  n.loci <- ncol(geno)/2
  if(n.loci != (floor(ncol(geno)/2)) )stop("Odd number of cols of geno")

# Check if Adjusted by Reg
  adjusted <- TRUE
  if( all(is.na(x.adj)) ) adjusted <- FALSE
  if(adjusted){
    x.adj <- as.matrix(x.adj)
    if(nrow(x.adj)!=length(y)) stop("Dims of y and x.adj are not compatible")
  }
     

# General checks for missing data
   miss <- is.na(y) 
   if(adjusted) miss <- miss | apply(is.na(x.adj),1,any)

# If Poisson, check for errors and missing offset
   if(trait.int==3) {
       if(all(is.na(offset))) stop("Missing offset")
       miss <- miss | is.na(offset)
   }

# Subset to non-missing values:
       y <- as.numeric(y[!miss])
       geno <- geno[!miss,]
       if(adjusted) x.adj <- x.adj[!miss,,drop=FALSE]
       if(trait.int==3) offset <- offset[!miss]

# Create a haplo object (using EMhaps)
   haplo <- haplo.em(geno, locus.label, miss.val=miss.val, 
                     control = em.control)

# Check convergence of EM
   if(!haplo$converge) stop("EM for haplo failed to converge")

# Need to check if any rows of geno were removed.
# this could be due to excessive number of haplotype pairs and/or 
# because rm.geno.na was chosen as true.
# If any removed, need to remove corresonding values
# in y and rows of x.adj, and offset for poisson
   rows.rem <- haplo$rows.rem
   if(length(rows.rem) > 0){
     keep <- !apply(outer((1:length(y)),rows.rem,"=="),1,any)
     y <- y[keep]
     if(adjusted) x.adj <- x.adj[keep,,drop=FALSE]
     if(trait.int==3) offset <- offset[keep]
   }

# If binomial, check for errors
   if(trait.int==2) {
      if(!all(y==1|y==0)) stop("Invalid y values")
      if(all(y==1) | all(y==0)) stop("No variation in y values")
    }

# If Proportional Odds, recode and check y values
   if(trait.int==4){
      y <- factor(y)
      y.lev <- levels(y)
      y <- as.numeric(y)
      if(max(y) < 3) stop("Less than 3 levels for y values")
    }

   n.subj <- length(y)

# Compute haplotype score statistics

   hap1 <- haplo$hap1code
   hap2 <- haplo$hap2code
   indx <- haplo$indx.subj
   post <- haplo$post
   nreps <- as.vector(haplo$nreps)
   uhap  <- sort(unique(c(hap1,hap2)))

# Don't score haplotypes that have probabilities <  skip.haplo

   which.haplo <- haplo$hap.prob >= skip.haplo
   uhap <- uhap[which.haplo]

   x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")

   n.x <- ncol(x)
   x.post <- matrix(rep(NA, n.subj * n.x), ncol=n.x)

   for(j in 1:n.x){
      x.post[,j] <- tapply(x[,j]*post, indx, sum)
   }

# Scores for GLM's

   if(trait.int <= 3){

     # If not adjusted, use summaries of y

     if(!adjusted){
        mu <- switch(trait.int, mean(y), mean(y), sum(y)/sum(offset) )
        a  <- switch(trait.int, var(y), 1, 1)

        # if not adjusted, still need to set up x.adj, for intercept,
        # to be used by haplo.score.glm to compute score stat adjusted
        # for intercept

        x.adj <- matrix(rep(1,n.subj),ncol=1)
     }

    # If adjusted, use GLM to get fitted and residuals
   
     if(adjusted){

         reg.out <- glm(y ~ x.adj, family=trait.type)

         # bind col of 1's for intercept:
         x.adj <- cbind(rep(1,n.subj),x.adj)

         mu <- reg.out$fitted.values
         a  <- switch(trait.int, 
                 sum(reg.out$residuals^2)/reg.out$df.residual,
                 1, 1)
      }

     v <- switch(trait.int, 1/a, mu*(1-mu), mu )
  
     # Now compute score statistics
 
     tmp <- haplo.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
     u.score <- tmp$u.score
     v.score <- tmp$v.score
   }


# Scores for Proportional Odds
   if(trait.int ==4) {

      if(adjusted){
         library("Design")
         library("Hmisc")
         reg.out <- lrm(y ~ x.adj)
         K <- max(y)
         n.xadj <- ncol(x.adj)
         alpha <- reg.out$coef[1:(K-1)]
         beta <- reg.out$coeff[K:(K-1 + n.xadj)]

         tmp <- haplo.score.podds(y, alpha, beta, x.adj, nreps, x.post,
                              post, x)
       }

      if(!adjusted){
         tbl <- table(y)
         s <- 1- (cumsum(tbl)-tbl)/n.subj
         alpha <-  - log((1-s[-1])/s[-1])
         tmp <- haplo.score.podds(y, alpha, beta=NA, x.adj=NA, nreps, x.post,
                              post, x)
       }

      u.score <- tmp$u.score
      v.score <- tmp$v.score
    }
   
# Compute Score Statistics:

   tmp <- Ginv(v.score)
   df <- tmp$rank
   g.inv <- tmp$Ginv
   score.global <- u.score%*% g.inv %*%u.score
   score.haplo <- u.score / sqrt(diag(v.score))
   score.max <-  max(score.haplo^2, na.rm=TRUE)

# Compute empirical p-values if simulate=TRUE

   if(!simulate) {
      score.global.p.sim <- NA
      score.haplo.p.sim <- rep(NA,length(score.haplo))
      score.max.p.sim <- NA
      n.val.global <- NA
      n.val.haplo <- NA
   } else {
         

      # initialize rejection counts
      score.global.rej <- 0
      score.haplo.rej  <- rep(0,length(score.haplo))
      score.max.rej    <- 0
      
      # initialize valid simulation counts. A simulation for the
      # global test can be invalid if the global statistic results 
      # in an NA value. A simulation for the max stat, or the individual
      # haplotype scores, is declared invalid if any of the individual
      # haplotype score stats results in an NA value. However, there is
      # an exception. If the observed score statistic for an individual
      # haplotype is NA, this haplotype score is ignored in the simulations,
      # allowing the other haplotypes to have simulated p-values, and
      # the observed max stat ignores the score relating to this
      # haplotype with an NA value.
 
      n.val.global <- 0
      n.val.haplo <- 0

      score.haplo.sqr <- score.haplo^2
      score.haplo.ok <- !is.na(score.haplo.sqr)
           
      if(trait.int<=3){

         # Initialize mu.rand and v.rand, in case not adjusted (if adjusted, 
         # then these will be over-ridden in simulation loop

         mu.rand <- mu
         v.rand <- v
       }

  ## Compute sequential Monte Carlo p-values as given by
  ## Besag and Clifford, Biometrika Jun 1991
  ##   Choose global.sim or max.sim as target pval to sample
  ##   within desired pval error range.
  ##   Employ an extra rule: sample a min (default=1000) for enough
  ##   samples to get good estimate of haplo.p.sim pvals

      done <- FALSE
      while(!done) {

         # random order
         rand.ord <- order(runif(n.subj))

         if(trait.int <=3){

           if(adjusted){
              mu.rand <- mu[rand.ord]
              v.rand <- switch(trait.int, v, v[rand.ord], v[rand.ord])
            }

           tmp <- haplo.score.glm(y[rand.ord], mu.rand, a, v.rand, 
                                x.adj[rand.ord,], nreps, x.post, post, x)
         }

         if(trait.int ==4){

            if(adjusted) {
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta,
                               x.adj[rand.ord,,drop=FALSE],nreps, x.post, post, x)
            }

            if(!adjusted) {
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta=NA,
                               x.adj=NA,nreps, x.post, post, x)
             }
          }

         u.score <- tmp$u.score
         v.score <- tmp$v.score
         
         # Now compute score statistics
     
         tmp <- Ginv(v.score)
         g.inv <- tmp$Ginv  
         score.global.sim <- u.score %*% g.inv %*% u.score

         # note that score.haplo.sim is squared, where score.haplo is not
         # because we want to keep the sign of score.haplo for returned 
         # values

         score.haplo.sim  <- ifelse(diag(v.score) < 0, NA,  u.score^2 / diag(v.score) )
                  
         if(!is.na(score.global.sim)) {
            n.val.global <- n.val.global +1
            if(score.global.sim >= score.global) score.global.rej <- score.global.rej +1
          }

        # for max stat, and individual haplotype scores, we require the same number of valid
        # simulations, so that the haplotype-specific simuated p-values are all based on the
        # same number of valid simulations. Without this, it would be impossible to know the 
        # confidence in the individual p-values.

         if(!any(is.na(score.haplo.sim))){

            n.val.haplo <- n.val.haplo + 1

            score.haplo.rej[score.haplo.ok] <- score.haplo.rej[score.haplo.ok] +
                               ifelse(score.haplo.sim[score.haplo.ok] >= score.haplo.sqr[score.haplo.ok], 1, 0)

             score.max.sim <- max(score.haplo.sim[score.haplo.ok])

            if(score.max.sim >= score.max) score.max.rej <- score.max.rej +1

          }
 
         # h is target count of rejection obs to stop sampling
         #   formula from Besag and Clifford
         h.global <- 1/(sim.control$p.threshold^2 + 1/max(1,n.val.global))
         h.max <- 1/(sim.control$p.threshold^2 + 1/max(1,n.val.haplo))
         
         # print the info out to screen
         if(sim.control$verbose) {
           cat("h.global:", round(h.global,5), ", global: ", score.global.rej, ", count: ", n.val.global, "\n")
           cat("h.max:", round(h.max,5), ", max-score: ", score.max.rej, ", count: ", n.val.haplo, "\n")
         }
         # done if reach max or enough rej values.  If max is reached, 
         # keep results and calculate p-vals at that point.
         if( (max(n.val.global,n.val.haplo) == sim.control$max.sim) |
             (h.global<=score.global.rej & h.max<=score.max.rej)    )
           done <- TRUE

         # not done if min.sim not met
         # this control is for haplo.sim values
         if(min(n.val.global, n.val.haplo) < sim.control$min.sim) done <- FALSE
         
       }

      score.global.p.sim <- score.global.rej /  n.val.global

      score.haplo.p.sim <- rep(NA, length(score.haplo.rej))
     
      score.haplo.p.sim[score.haplo.ok] <- score.haplo.rej[score.haplo.ok] / n.val.haplo

      score.max.p.sim <- score.max.rej / n.val.haplo
      
    }


   score.global.p <- 1 - pchisq(score.global,df)
   score.haplo.p <- 1-pchisq(score.haplo^2,1)

# Create locus label if missing:
   if(all(is.na(locus.label))) {
      locus.label<- paste("loc-",1:n.loci,sep="")
    }

   obj <- (list(score.global=score.global, df=df,score.global.p=score.global.p,
       score.global.p.sim=score.global.p.sim,
       score.haplo=score.haplo,score.haplo.p=score.haplo.p,
       score.haplo.p.sim=score.haplo.p.sim,
       score.max.p.sim=score.max.p.sim,
       haplotype=haplo$haplotype[which.haplo,],
       hap.prob=haplo$hap.prob[which.haplo],
       locus.label=locus.label, simulate=simulate,
       sim.control=sim.control, n.val.global=n.val.global,
       n.val.haplo=n.val.haplo, rows.rem=rows.rem))

   if(exists("is.R") && is.function(is.R) && is.R()) {
     class(obj) <- "haplo.score"
   } else {
     oldClass(obj) <- "haplo.score"
   }

   return(obj)
   
 }


