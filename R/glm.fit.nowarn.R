#$Author: sinnwell $
#$Date: 2004/01/15 22:09:38 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/glm.fit.nowarn.R,v 1.2 2004/01/15 22:09:38 sinnwell Exp $
#$Locker:  $
#$Log: glm.fit.nowarn.R,v $
#Revision 1.2  2004/01/15 22:09:38  sinnwell
#comment changes from glm.fit.R
#
#Revision 1.1  2004/01/12 14:47:45  sinnwell
#Initial revision
#

glm.fit.nowarn <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = glm.control(), intercept = TRUE) 

{
  # complete copy of glm.fit for haplo.glm
  # use this function when family is binom to bypass expected warnings
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    ynames <- names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- NCOL(x)
    if (nvars == 0) {
        cc <- match.call()
        cc[[1]] <- as.name("glm.fit.null")
        return(eval(cc, parent.frame()))
    }
    if (is.null(weights)) 
        weights <- rep(1, nobs)
    if (is.null(offset)) 
        offset <- rep(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE

    ## ----- JPS 1/15/2004----------------------------------
    ## warning expected when family$initialize is evaluated
    ## for non-integer counts, turn off warnings while initializing
    old.warn <- options("warn")
    options(warn=-1)

    if (is.null(mustart)) {
        eval(family$initialize)
        
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    
    ## restore warn option to what it was before.
    options(old.warn)
    ## --------  only changes from glm.fit.R -----------------    

    if (NCOL(y) > 1) 
        stop("y must be univariate unless binomial")
    coefold <- NULL
    eta <- if (!is.null(etastart)) 
        etastart
    else if (!is.null(start)) 
        if (length(start) != nvars) 
            stop("Length of start should equal ", nvars, " and correspond to initial coefs for ", 
                deparse(xnames))
        else {
            coefold <- start
            offset + as.vector(if (NCOL(x) == 1) 
                x * start
            else x %*% start)
        }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
        stop("Can't find valid starting values: please specify some")
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- FALSE
    for (iter in 1:control$maxit) {
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (any(is.na(varmu))) 
            stop("NAs in V(mu)")
        if (any(varmu == 0)) 
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) 
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            warning("No observations informative at iteration ", 
                iter)
            break
        }
        z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        ngoodobs <- as.integer(nobs - sum(!good))
        fit <- .Fortran("dqrls", qr = x[good, ] * w, n = ngoodobs, 
            p = nvars, y = w * z, ny = as.integer(1), tol = min(1e-07, 
                control$epsilon/1000), coefficients = double(nvars), 
            residuals = double(ngoodobs), effects = double(ngoodobs), 
            rank = integer(1), pivot = 1:nvars, qraux = double(nvars), 
            work = double(2 * nvars), PACKAGE = "base")
        if (any(!is.finite(fit$coefficients))) {
            conv <- FALSE
            warning("Non-finite coefficients at iteration ", 
                iter)
            break
        }
        if (nobs < fit$rank) 
            stop("X matrix has rank ", fit$rank, " but only ", 
                nobs, " observations")
        start[fit$pivot] <- fit$coefficients
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold)) 
                stop("no valid set of coefficients has been found:please supply starting values", 
                  call. = FALSE)
            warning("Step size truncated due to divergence", 
                call. = FALSE)
            ii <- 1
            while (!is.finite(dev)) {
                if (ii > control$maxit) 
                  stop("inner loop 1; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta <- drop(x %*% start)
                mu <- linkinv(eta <- eta + offset)
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
            if (control$trace) 
                cat("Step halved: new deviance =", dev, "\n")
        }
        if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit) 
                  stop("inner loop 2; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta <- drop(x %*% start)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Step halved: new deviance =", dev, "\n")
        }
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        else {
            devold <- dev
            coef <- coefold <- start
        }
    }
    if (!conv) 
        warning("Algorithm did not converge")
    if (boundary) 
        warning("Algorithm stopped at boundary value")
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) 
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
        if (any(mu < eps)) 
            warning("fitted rates numerically 0 occurred")
    }
    if (fit$rank < nvars) 
        coef[fit$pivot][seq(fit$rank + 1, nvars)] <- NA
    xxnames <- xnames[fit$pivot]
    residuals <- rep(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
        Rmat <- diag(nvars)
        Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
    }
    else Rmat <- fit$qr[1:nvars, 1:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    names(fit$effects) <- c(xxnames[seq(fit$rank)], rep("", sum(good) - 
        fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    resdf <- n.ok - fit$rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * fit$rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = fit$effects, R = Rmat, rank = fit$rank, qr = structure(fit[c("qr", 
            "rank", "qraux", "pivot", "tol")], class = "qr"), 
        family = family, linear.predictors = eta, deviance = dev, 
        aic = aic.model, null.deviance = nulldev, iter = iter, 
        weights = wt, prior.weights = weights, df.residual = resdf, 
        df.null = nulldf, y = y, converged = conv, boundary = boundary)
}
