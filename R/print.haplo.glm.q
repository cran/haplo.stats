#$Author: sinnwell $
#$Date: 2005/03/29 14:18:47 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/print.haplo.glm.q,v 1.9 2005/03/29 14:18:47 sinnwell Exp $
#$Locker:  $
#$Log: print.haplo.glm.q,v $
#Revision 1.9  2005/03/29 14:18:47  sinnwell
#fix call to print within widths, different for R/Splus
#®
#
#Revision 1.8  2004/10/22 22:08:27  sinnwell
#do not drop matrix to vector when subsetting to haplo.unique
#when only 1 haplotype
#
#Revision 1.7  2004/03/18 23:30:32  sinnwell
#keep matrix from converting to data.frame, and char vecs to factors
#
#Revision 1.6  2004/02/26 23:05:23  sinnwell
#print.banner to printBanner
#
#Revision 1.5  2004/02/06 16:34:17  sinnwell
#fix 1-sided pval to 2-sided
#
#Revision 1.4  2003/12/08 20:16:49  sinnwell
# changed T,F to TRUE,FALSE
#
#Revision 1.3  2003/11/17 23:28:19  schaid
#made compatible with R
#
#Revision 1.2  2003/10/15 21:13:30  schaid
#got rid of bug caused by use of 'fit' (should have been x)
#
#Revision 1.1  2003/09/16 16:03:15  schaid
#Initial revision
#
print.haplo.glm <- function(x, print.all.haplo=FALSE, digits = max(options()$digits - 4, 3), ...){

  if(exists("is.R") && is.function(is.R) && is.R()) {
    x$call <- deparse(x$call, width.cutoff=40)
    cat("\n  Call: ", x$call, sep="\n")
  }
  else {
    cat("\n  Call: \n")
    dput(x$call)
  }
  
  haplo.df<- function(x){
    z <- x$haplo.common
    df <- as.matrix(x$haplo.unique[z,,drop=FALSE])
    y <- x$haplo.freq[z]

    if(x$haplo.rare.term){
      df <- rbind(df, rep("*",ncol(df)))
      y <- c(y, sum(x$haplo.freq[x$haplo.rare]))
    }

    # use dimnames to change row names do not convert from matrix to df
    dimnames(df)[[1]] <- x$haplo.names
    df <- rbind(df,x$haplo.unique[x$haplo.base,])
    dimnames(df)[[1]][nrow(df)] <- "haplo.base"
    y <- c(y,x$haplo.freq[x$haplo.base])
    data.frame(df,hap.freq=y)
  }


  ncoef <- length(x$coef)
  coef <- x$coef
  se <- sqrt(x$var.mat[cbind(1:ncoef, 1:ncoef) ])

  wt <- x$weights.expanded * x$haplo.post.info$post
  df.residual <- sum(wt) - length(x$coef) 

  t.stat <- coef/se
  pval <- 2*(1-pt(abs(t.stat),  df.residual))

#  printBanner("Regression Coefficients")
  cat("\nCoefficients:\n")
  print(cbind(coef=coef, se=se, t.stat=t.stat, pval=pval), digits=digits)

#  cat("\n")
#  printBanner("Hapoltypes and their Frequencies")

  cat("\nHaplotypes:\n")
  print(haplo.df(x), digits=digits)
  

  if(print.all.haplo){
    haplo.type <- rep(NA,length(x$haplo.freq))
    haplo.type[x$haplo.common] <- "C"
    haplo.type[x$haplo.rare] <- "*"
    haplo.type[x$haplo.base] <- "B"
    df <- data.frame(x$haplo.unique, hap.freq = round(x$haplo.freq, digits), hap.type=haplo.type)
    cat("\n")
    printBanner("All Haplotypes")
    cat("B = base   haplotype\n")
    cat("C = common haplotype\n")
    cat("* = rare   haplotype\n\n")

    print(df)
  }
}

