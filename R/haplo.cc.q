#$Author: sinnwell $
#$Date: 2005/03/10 17:36:16 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.cc.q,v 1.7 2005/03/10 17:36:16 sinnwell Exp $
#$Locker:  $
#$Log: haplo.cc.q,v $
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

haplo.cc <- function(y, geno, haplo.min.count=5, locus.label=NA, ci.prob=0.95,
                     miss.val=c(0,NA), simulate=FALSE,
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
 
  exclude <- (1:length(y))[is.na(y)]
  tmp <- (1:n.subj)[apply(is.na(geno),1,all)]
  exclude <- unique(c(exclude,tmp))
  if(length(exclude)) {
    geno <- geno[-exclude,]
    y <- y[-exclude]
  }
  
  # check y to be 0's and 1's
  if(!all(y==0 | y==1))
    stop("y must be binary with 0=controls and 1=cases")
  
  n.subj <- nrow(geno)
  
  # set skip.haplo and hap.freq.min for haplo.glm to be the same, meeting
  # min. number of expected counts per haplotypes
  if(haplo.min.count < 1)
    warning("haplo.min.count may be too small for reliable results\n")
  
  skip.haplo <- control$haplo.freq.min <- haplo.min.count/(2*n.subj)
  control$haplo.min.count=haplo.min.count
  
  score.lst <- haplo.score(y=y, geno=geno, trait.type="binomial",
                           skip.haplo=skip.haplo, locus.label=locus.label,
                           miss.val=miss.val, simulate=simulate,
                           sim.control=sim.control, em.control=control$em.c)
                           
  # get haplotype frequency estimates within the two groups
  group.lst <- haplo.group(group=y, geno=geno, locus.label=locus.label,
                           miss.val=miss.val, control=control$em.c)
  
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
  save.alleles <- attributes(geno.glm)$unique.alleles
#  oldClass(geno) <- "model.matrix"
  glm.data <- data.frame(geno.glm, y=y)

# if (is.na(weights)) weights=rep(1,n.subj)
  
  # will only get OR's for haplotypes with freq > haplo.freq.min, based on min.count
  # should be set same as skip.haplo so 1:1 merge of haps from merge.lst
  fit.lst <- haplo.glm(y~geno.glm, family=binomial, data=glm.data, #weights=weights,
                       locus.label=locus.label, allele.lev=save.alleles,
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
                 fit.lst=fit.lst, ci.prob=ci.prob)
  
  ## returned score and glm objects in the returned list for
  ## additional applications
  
  oldClass(cc.lst) <- "haplo.cc"
  cc.lst
  
}