#$Author: sinnwell $
#$Date: 2011/11/10 15:29:40 $
#$Header: /projects/genetics/cvs/cvsroot/haplo.stats/R/haplo.cc.q,v 1.14 2011/11/10 15:29:40 sinnwell Exp $
#$Locker:  $
#$Log: haplo.cc.q,v $
#Revision 1.14  2011/11/10 15:29:40  sinnwell
#major update to hapglm, minor changes to Rd files, prepare for version 1.5.0 release
#
#Revision 1.13  2008/04/10 20:58:36  sinnwell
#add eps.svd, only allow haplo.min.count in control()
#
#Revision 1.12  2007/03/30 20:02:39  sinnwell
#control random seed values before each call
#
#Revision 1.11  2007/03/23 19:08:37  sinnwell
#take out weights for haplo.glm, they don't work for haplo.glm in Splus
#
#Revision 1.10  2007/03/22 20:08:46  sinnwell
#insert line to check weights for NULL before passing to haplo.glm
#
#Revision 1.9  2007/03/13 18:14:40  sinnwell
#add weights parameter
#
#Revision 1.8  2005/06/03 14:32:58  sinnwell
#adjust haplo.rare if l.t. haplo.min.info
#
#Revision 1.7  2005/03/10 17:36:16  sinnwell
#make haplo.min.count persist in haplo.glm
#
#Revision 1.6  2005/01/04 20:22:29  sinnwell
#min.count to haplo.min.count
#
#Revision 1.5  2004/12/02 15:56:12  sinnwell
#logical vectors y==0 + y==1 to y==0 | y==1. R doesn't add logicals.
#
#Revision 1.4  2004/11/10 21:16:17  sinnwell
#problems with setupGeno before haplo.glm
#
#Revision 1.3  2004/07/13 22:02:27  sinnwell
#T to TRUE
#
#Revision 1.2  2004/06/09 18:57:55  sinnwell
#fix for when glm.object$haplo.rare is zero length
#
#Revision 1.1  2004/04/23 21:24:56  sinnwell
#Initial revision
#


## Sinnwell JP and Schaid DJ
## Mayo Clinic, Rochester MN
## 1/2004

haplo.cc <- function(y, geno, locus.label=NA, ci.prob=0.95,
                     miss.val=c(0,NA), weights=NULL, eps.svd=1e-5, simulate=FALSE,
                     sim.control=score.sim.control(),
                     control=haplo.glm.control())

{
  ## wrap the execution of haplo.score, haplo.group and haplo.glm
  ## into this function, specific for case-control (cc) studies.  

  ## 1) calculate Odds Ratios for haplotype effects.
  ## 2) keep in mind it is a ratio relative to the baseline
  ##    haplotype, which is often highly associated with the trait,
  ##    shown by freq's and the score.
  ## 3) You may choose different, 'neutral' baseline by setting
  ##    haplo.base to the index of haplo.unique from a previous run with same 
  ##    haps kept by haplo.freq.min, or haplo.min.count cutoffs
  
  # check for missing values in y and geno
  n.subj <- length(y)
  if(n.subj != nrow(geno))
    stop("Different number of rows for y and geno!")
 
  exclude <- which(is.na(y))
  tmp <- which(apply(is.na(geno),1,all))
  exclude <- unique(c(exclude,tmp))
  if(length(exclude)) {
    warning(paste("Subjects ", exclude,
                  " removed for missing y or all alleles missing."))
    geno <- geno[-exclude,]
    y <- y[-exclude]
  }

  if(is.null(control$em.control$iseed)) {
    runif(1) # use a random number, so the next .Random.seed value is used
    control$em.control$iseed <- .Random.seed
  }
 
  # check y to be 0's and 1's
  if(!all(y==0 | y==1))
    stop("y must be binary with 0=controls and 1=cases")
  
  n.subj <- nrow(geno)
  
  # set skip.haplo and hap.freq.min for haplo.glm to be the same, meeting
  # min. number of expected counts per haplotypes
  haplo.min.count <- control$haplo.min.count
  haplo.freq.min <- control$haplo.freq.min
  if(is.na(haplo.min.count)) haplo.min.count <- haplo.freq.min*2*nrow(geno)
  
  if( !is.null(haplo.min.count) && haplo.min.count < 1)
    warning("haplo.min.count may be too small for reliable results\n")

  set.seed(control$em.control$iseed)
  score.lst <- haplo.score(y=y, geno=geno, trait.type="binomial",
                        min.count=haplo.min.count, 
                        locus.label=locus.label,
                        miss.val=miss.val, haplo.effect=control$haplo.effect,
                        eps.svd=eps.svd, simulate=simulate,
                        sim.control=sim.control, em.control=control$em.c)
                           
  # get haplotype frequency estimates within the two groups
  set.seed(control$em.control$iseed)
  group.lst <- haplo.group(group=y, geno=geno, locus.label=locus.label,
                           miss.val=miss.val, weight=weights,
                           control=control$em.c)
  
  group.count <- group.lst$group.count
  names(group.count) <- c("control", "case")
  n.loci <- group.lst$n.loci
  
  # merge group and score data
  merge.df <- haplo.score.merge(score.lst, group.lst)

  names(merge.df)[(n.loci+score.lst$simulate+3):(n.loci+score.lst$simulate+5)] <-
                 c("pool.hf", "control.hf", "case.hf")
 
  # prepare data for haplo.glm with binomial family logit link.
  # extract for Odds Ratios, exp(betas)
  geno.glm <- setupGeno(geno, miss.val=miss.val)
  
  glm.data <- data.frame(geno.glm, y=y)

# weights are not supported for haplo.glm in S-PLUS--take out for now
#  if (is.null(weights)) weights=rep(1,n.subj)
   
  # will only get OR's for haplotypes with freq > haplo.freq.min, based on min.count
  # should be set same as skip.haplo so 1:1 merge of haps from merge.lst
  set.seed(control$em.control$iseed)
  fit.lst <- haplo.glm(y~geno.glm, family=binomial, data=glm.data, 
                       locus.label=locus.label,
                       miss.val=miss.val, control=control)

  ncoef <- length(fit.lst$coef)
  
  # get OR's and variance estimates based on var(b's)
  # make base haplo coef zero b/c or will be 1
  hap.coef <- (fit.lst$coefficients)[2:ncoef]
  se.bet <- sqrt(fit.lst$var.mat[cbind(2:ncoef, 2:ncoef)])

  # find CI's for OR's transormed from CI's on Betas 
  if(0 > ci.prob || ci.prob > 1) {
    warning("invalid value for ci.prob, setting to default 0.95")
    ci.prob <- 0.95
  }

  norm.constant <- -1*qnorm((1-ci.prob)/2)
  lower.ci <- hap.coef - norm.constant*se.bet
  lower.ci <- c(NA,exp(lower.ci))

  upper.ci <- hap.coef + norm.constant*se.bet
  upper.ci <- c(NA,exp(upper.ci))
  hap.OR <- c(1,exp(hap.coef))
  OR.df <- cbind(lower.ci, hap.OR, upper.ci)
  
  ## when a rare haplotype is below haplo.min.info, & not part of model
  ## need to remove from haplo.rare    --JPS 6/2005
  fit.lst$haplo.rare <- fit.lst$haplo.rare[fit.lst$haplo.freq[fit.lst$haplo.rare] > control$haplo.min.info]

  ## combine haps, OR's and CI's in a data.frame 
  fit.df <- as.data.frame(fit.lst$haplo.unique)[c(fit.lst$haplo.base,
                                          fit.lst$haplo.common, fit.lst$haplo.rare),]

  # if rare term replicate rare coeff for every rare hap
  if(length(fit.lst$haplo.rare) > 0) {
    rare.df <- matrix(rep(OR.df[ncoef,],length(fit.lst$haplo.rare)),ncol=3,byrow=TRUE)
    if(control$keep.rare.haplo) {
      OR.df <- rbind(OR.df[-ncoef,],rare.df)
    } else OR.df <- rbind(OR.df,rare.df)
  }
  
  class.vec <- c("Base", rep("Eff",length(fit.lst$haplo.common)), rep("R", length(fit.lst$haplo.rare)))
  fit.df <- cbind(fit.df,class.vec,OR.df)
  names(fit.df)[(n.loci+1):(n.loci+4)] <- c("glm.eff","OR.lower","OR", "OR.upper")

  # format the hap columns so the df's will merge on the haps
  fit.df[,1:n.loci] <- format(fit.df[,1:n.loci])
  merge.df[,1:n.loci] <- format(merge.df[,1:n.loci])

  cc.df <- merge(merge.df, fit.df,by=1:n.loci,all.x=TRUE,all.y=TRUE)
  cc.lst <- list(cc.df=cc.df, group.count=group.count, score.lst=score.lst,
                 fit.lst=fit.lst, ci.prob=ci.prob, exclude.subj=exclude)
  
  ## returned score and glm objects in the returned list for
  ## additional applications
  
  class(cc.lst) <- "haplo.cc"
  cc.lst
  
}
