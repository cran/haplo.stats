#$Author: sinnwell $
#$Date: 2005/01/25 23:04:06 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.glm.q,v 1.12 2005/01/25 23:04:06 sinnwell Exp $
#$Locker:  $
#$Log: haplo.glm.q,v $
#Revision 1.12  2005/01/25 23:04:06  sinnwell
#use as.list(args(haplo.glm)) instead of formals when na.action not given
#formals is for R only
#
#Revision 1.11  2005/01/25 22:25:22  sinnwell
#take out print statements left from debugging
#
#Revision 1.10  2005/01/23 20:01:38  sinnwell
#use Thomas Lumley hint from R-help, use formals of haplo.glm to set
#default na.action to na.geno.keep.
#
#Revision 1.9  2005/01/19 16:14:39  sinnwell
#default na.action not set to na.geno.keep, force it if missing(na.action)
#m$na.action <- as.name("na.geno.keep")
#
#Revision 1.8  2005/01/04 17:37:22  sinnwell
#set haplo.min.freq as a function of haplo.min.count
#using a min expected count of haplotypes is better criteria
#
#Revision 1.7  2004/03/19 15:01:48  sinnwell
#consider .C(PACKAGE= as part of '...'
#
#Revision 1.6  2004/03/17 21:11:12  sinnwell
#separate calls of .C("groupsum" for R and Splus
#
#Revision 1.5  2004/03/03 22:13:16  schaid
#added allele.lev to pass allele labels, to allow haplo.glm to work in R for character alleles.
#
#Revision 1.4  2004/01/12 14:46:57  sinnwell
#get around warnings for binomial trait in R
#
#Revision 1.3  2003/12/08 19:37:09  sinnwell
#changed F,T to FALSE,TRUE
#
#Revision 1.2  2003/11/17 23:27:20  schaid
#made compatible with R
#
#Revision 1.1  2003/09/16 16:01:54  schaid
#Initial revision
#
haplo.glm     <- function(formula = formula(data),
                          family = gaussian, 
                          data = sys.parent(),
                          weights,
                          na.action="na.geno.keep",
                          start = eta,
                          miss.val=c(0,NA),
                          locus.label=NA,
                          allele.lev = NULL,
                          control = haplo.glm.control(),
                          method = "glm.fit",
                          model = FALSE,
                          x = FALSE,
                          y = TRUE,
                          contrasts = NULL,
                          ...){
  call <- match.call()
  
  m <- match.call(expand.dots = FALSE)

  if(is.null(m$na.action)) m$na.action=as.list(args(haplo.glm))$na.action

  m$family <- m$miss.val <-  m$locus.label <- m$control <- m$method <- m$model <- m$x <- m$y <- NULL
  m$contrasts <- m$print.iter <- m$allele.lev <- m$... <- NULL

  m$drop.unused.levels <- TRUE

  m[[1]] <- as.name("model.frame")
  
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")

  # set design matrix contrasts
  contrasts.tmp <- c(factor="contr.treatment", ordered = "contr.poly")
  contrasts.old <- options()$contrasts
  options(contrasts = contrasts.tmp)
  on.exit(options(contrasts=contrasts.old))

  # Extract  weights per subject

  if(exists("is.R") && is.function(is.R) && is.R()) {
    wt.subj <- model.weights(m)  
  } else
  {
    wt.subj <- model.extract(m, weights)
  }

  if(!length(wt.subj))
    wt.subj <- rep(1, nrow(m))
  else if(any(wt.subj < 0))
    stop("negative weights not allowed")

  # for control, base minimum haplotype frequencies on haplo.min.count
  # as precedent over haplo.freq.min, the min count is 5
  if(is.na(control$haplo.freq.min))
    control$haplo.freq.min <- control$haplo.min.count/(2*nrow(data))
  
  ####### Modify model.frame ####################################################
  # Translate from unphased genotype matrix to haplotype design matrix
  # Use haplo.model.frame to modify the model.frame. Note that haplo.model.frame
  # calls haplo.em to estimate haplotype frequencies.

  haplo.mf <- haplo.model.frame(m, locus.label=locus.label, allele.lev = allele.lev, 
                                miss.val=miss.val, control=control)


  # Now extract the model.frame created by haplo.model.frame
  m <- haplo.mf$m.frame

  if(method == "model.frame")
    return(m)

  # g.dat is a dataframe with indx.subj, hap1, hap2, where hap1 and hap2 are
  # numeric codes for haplotypes

  g.dat           <- haplo.mf$g.dat
  haplo.unique    <- haplo.mf$haplo.unique
  haplo.base      <- haplo.mf$haplo.base
  haplo.freq  <- haplo.freq.init <- haplo.mf$haplo.freq
  haplo.common    <- haplo.mf$haplo.common
  haplo.rare.term <- haplo.mf$haplo.rare.term
  haplo.rare      <- haplo.mf$haplo.rare

 
  # set up integer grouping variables (indices) for group sums

  haplo.group <- c(g.dat$hap1,g.dat$hap2)
  len.haplo.group <- length(haplo.group)
  n.haplo.group <- length(haplo.freq)
  subj.indx <- as.numeric(factor(g.dat$indx.subj))
  len.subj.indx <- length(subj.indx)
  n.subj <- length(unique(subj.indx))
 
  # Setup weights after expanding each subject into pairs of haplotypes
  # (i.e., possibly replicate weights according to replicates required
  # for enumerating pairs of haplotypes)

  wt.expanded <- wt.subj[subj.indx]

  # compute priors and posteriors of pairs of haplotypes, given observed data.

  prior.coef <- ifelse(g.dat$hap1!=g.dat$hap2, 2, 1)
  prior <- prior.coef * haplo.freq[g.dat$hap1] * haplo.freq[g.dat$hap2]

  pr.pheno.init <- tapply(prior,g.dat$indx.subj,sum)
  nreps <- tapply(g.dat$indx.subj,g.dat$indx.subj,length)
  den <- rep(pr.pheno.init,nreps)
  post <- prior/den

  lnlike.haplo <- sum(wt.subj*log(pr.pheno.init))
  

  # Keep initial post to return

  post.initial <- post
      
  # Keep count of true number of haplotypes for later use as denom 
  # of haplotype frequenices 

  n.hap <- 2 * sum(wt.subj)

  # First E-step: create new weights based on post. 

  # Note that these weights are based on post, which was computed assuming
  # that there are no regression terms (i.e., all regression
  # beta's = 0).

  # wt.expanded is the original input weight vector, wt.subj, 
  # (or set to 1's if not input),
  # after enumeration of all possible pairs of haplos.
  # Need to distinguish between wt.expanded and
  # working weights, which are wght.expanded * post

  w <- wt.expanded * post

  xvars <- as.character(attr(Terms, "variables"))

  if(exists("is.R") && is.function(is.R) && is.R()) {
     xvars <- xvars[-1]
     if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
   } 


  if(length(xvars) > 0) {
    xlevels <- lapply(m[xvars], levels)
    xlevels <- xlevels[!sapply(xlevels, is.null)]
    if(length(xlevels) == 0)
      xlevels <- NULL
  }
  else xlevels <- NULL

  a <- attributes(m)


  if(exists("is.R") && is.function(is.R) && is.R()) {
    Y <- model.response(m, "numeric")
   } else{
    Y <- model.extract(m, response)
  }


  X <- model.matrix(Terms, m, contrasts)

  start <- model.extract(m, start)

  if(exists("is.R") && is.function(is.R) && is.R()) {
      offset <- model.offset(m)
   } else{
      offset <- model.extract(m, offset)
  }


 if(exists("is.R") && is.function(is.R) && is.R()) {
   if (is.character(family)) 
        family <- get(family)
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
   } else{ 
    family <- as.family(family)
  }


  if(missing(method))
    method <- attr(family, "method")
  if(!is.null(method)) {
    if(!existsFunction(method))
      stop(paste("unimplemented method:", method))
  }
  else method <- "glm.fit"

  glm.fitter <- get(method)
   
  # Unlike the generic glm, we compute the Null model (Intercept-only) 
  # here,  instead of within the EM loop (to avoid recomputing the null
  # within each loop iteration). 

 

  if(exists("is.R") && is.function(is.R) && is.R()) {
    
    # switch fitter for binomial because it gives unneeded warnings
    if(family$family=="binomial") glm.fitter=glm.fit.nowarn
    
    fit.null <- glm.fitter(x = X[, "(Intercept)", drop = FALSE],
                            y = Y, 
                            weights = w, 
                            etastart = start, 
                            offset = offset, 
                            family = family, 
                            control=glm.control(maxit = control$maxit, epsilon = control$epsilon))
   } else{
     fit.null <- glm.fitter(x = X[, "(Intercept)", drop = FALSE],
                            y = Y, 
                            w = w,
                            start = start,
                            offset = offset, 
                            family = family, 
                            maxit = control$maxit, 
                            epsilon = control$epsilon)   
  }



  dfit <- dglm.fit(fit.null)

  prior <- prior*dfit
  pr.pheno <- tapply(prior, subj.indx, sum)
  lnlike.old <- lnlike.null <- sum(wt.subj*log(pr.pheno))


  # EM loop

  # Set up arrays for group sums within EM loop

  prior.tot <- rep(0,n.subj)
  post.tot <- rep(0,n.haplo.group)

  converge.em <- FALSE

  # need a separate control.epsilon for EM, because need more
  # stringent check for change in lnLike for EM than for 
  # glm.fit

  control.epislon.em <- min(c(control$epsilon,  0.000001))
  iter <- 0


  while(iter < control$maxit){

   iter <- iter + 1

   # M-step for regression beta's, using weights that depend on
   # earlier haplotype freqs and beta's



  if(exists("is.R") && is.function(is.R) && is.R()) {
     fit <- glm.fitter(x = X,
                       y = Y, 
                       weights = w, 
                       etastart = start,
                       offset = offset, 
                       family = family, 
                       control=glm.control(maxit = control$maxit, epsilon = control$epsilon))
   } else{
     fit <- glm.fitter(x = X,
                       y = Y,
                       w = w,
                       start = start,
                       offset = offset,
                       family = family,
                       maxit = control$maxit,
                       epsilon = control$epsilon,
                       trace = control$trace,
                       null.dev = NULL, qr=FALSE)
     
   }


    # Using new beta's, but earlier haplotype frequencies, update post, and 
    # compute lnLike

    dfit <- dglm.fit(fit)
    prior <- prior.coef * haplo.freq[g.dat$hap1] * haplo.freq[g.dat$hap2] * dfit

    # For the following, the use of tapply was originally used, 
    # as pr.pheno <- tapply(prior,subj.indx, sum), but this took too much
    # time within this EM loop, so a C function 'groupsum' is used instead.

    tmp.sum <- .C("groupsum",
                  x=as.double(prior),
                  indx=as.integer(subj.indx),
                  n=as.integer(len.subj.indx),
                  grouptot= as.double(prior.tot),
                  ngroup=as.integer(n.subj),
                  PACKAGE="haplo.stats")

    pr.pheno <- tmp.sum$grouptot
    den <- rep(pr.pheno,nreps)
    post <- prior/den
    lnlike.new <- sum(wt.subj*log(pr.pheno))


    # Create new weights based on post. These posteriors depend
    # on the ith iteration of beta's, but the (i-1)th iteration of
    # haplo.freq

    w <- wt.expanded * post

    # M-step for  haplotype frequencies, based on expected counts

   tmp.sum <- .C("groupsum",
                 x=as.double(c(post*wt.expanded, post*wt.expanded)),
                 indx=as.integer(haplo.group),
                 n=as.integer(len.haplo.group),
                 grouptot= as.double(post.tot),
                 ngroup=as.integer(n.haplo.group),
                 PACKAGE="haplo.stats")

   e.hap.count <- tmp.sum$grouptot
   haplo.freq <- e.hap.count/n.hap

    # convergence checks
    if(abs(lnlike.new - lnlike.old) < control.epislon.em){
      converge.em <- TRUE
      break
    }

    lnlike.old <- lnlike.new

  }


  if(!converge.em) warning("Failed to converge during EM-glm loop")

  lnlike.final <- lnlike.new

  fit$iter <- iter
  fit$lnlike  <- lnlike.final
  fit$lnlike.null <-  lnlike.null
  fit$lrt <- list(lrt = 2*(lnlike.final - lnlike.null), df = fit$rank-1 )

  fit$weights.expanded <- wt.expanded

  # Because prior.weights and weights depend only on the last iter
  # of the updated regression coefficients, they are not returned.

  fit$prior.weights <- NULL
  fit$weights <- NULL

  if(!is.null(xlevels))
    attr(fit, "xlevels") <- xlevels
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$call <- call
  if(model)
    fit$model <- m
  if(!y)
    fit$y <- NULL
  fit$control <- control
  if(!is.null(attr(m, "na.action")))
    fit$na.action <- attr(m, "na.action")

  fit$haplo.unique <- haplo.unique
  fit$haplo.base <- haplo.base
  fit$haplo.freq <- haplo.freq
  fit$haplo.freq.init <- haplo.freq.init
  fit$converge.em <- converge.em
  fit$haplo.common <- haplo.common
  fit$haplo.rare <- haplo.rare
  fit$haplo.rare.term <- haplo.rare.term
  fit$haplo.names <- haplo.mf$haplo.names

  # data.frame for info on posteriors
  haplo.post.info <- cbind(g.dat, post.initial, post)
  names(haplo.post.info) <- c("indx", "hap1", "hap2", "post.init", "post")
  fit$haplo.post.info <- haplo.post.info

  # estimation of the information and variance matrix.
  # requires X matrix, so attach to fit, and deattach
  # later if not desired to be returned

  fit$x <- X

  tmp <- louis.info(fit)
  fit$info <- tmp$info
  fit$var.mat <- tmp$var.mat
  fit$haplo.elim <- tmp$haplo.elim
  fit$rank.info <- tmp$rank
  
  if(!x) fit$x <- NULL

  if(exists("is.R") && is.function(is.R) && is.R()) {
     class(fit) <- "haplo.glm"
   } else {
     oldClass(fit) <- "haplo.glm"
   }
  
  return(fit)

}


