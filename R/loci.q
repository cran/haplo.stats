#$Author: sinnwell $
#
#$Date: 2003/12/08 19:49:36 $
#
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/loci.q,v 1.4 2003/12/08 19:49:36 sinnwell Exp $
#
#$Locker:  $
#
#$Log: loci.q,v $
#Revision 1.4  2003/12/08 19:49:36  sinnwell
#changed T,F to TRUE,FALSE
#
#Revision 1.3  2003/09/19 15:01:56  sinnwell
#add class() for R dual compatibility
#
#Revision 1.2  2003/02/27 16:06:06  det01
#Removed setOldClass(c("loci", "model.matrix")) line of code
#and modified oldClass(loci) <- "loci" to oldClass(loci) <- "model.matrix"
#
#Revision 1.1  2002/12/13 17:52:43  det01
#Initial revision
#
#
loci <- function(geno,
                  locus.names,
                  chrom.label=NULL,
                  x.linked=FALSE, 
                  sex=NULL,
                  male.code="M",
                  female.code="F",
                  miss.val=NA, 
                  map=NA) {

# Title: Create an object of class model.matrix

   loci <- NULL
   unique.alleles <- list(NULL)

   if (missing(geno))
      stop("Error: Required argument geno is missing")

   if(!is.data.frame(geno)&!is.matrix(geno))
      stop("Error: Argument geno is not a data.frame or matrix")

   if (missing(locus.names))
      stop("Error: Required argument locus.names is missing")

   if (ncol(geno)%%2)
      stop("Error: Number of columns in geno object is not an even number")
  
    exp.num.loc <- ncol(geno)/2
    num.loc.names <- length(locus.names)
    if (exp.num.loc!=num.loc.names)
       stop("Error: The locus.names vector length is not equal to the expected number of loci")

    for(i in 1:exp.num.loc){
       tmp <- locus(geno[,(i*2-1)],geno[,i*2],chrom.label=chrom.label,locus.names[i],x.linked=x.linked,sex=sex,male.code=male.code,female.code=female.code,miss.val=miss.val)
       ualleles <- attr(tmp,"allele.labels")
       tmp <- as.matrix(tmp)
       dimnames(tmp)[[2]][1] <- paste(locus.names[i],dimnames(tmp)[[2]][1],sep=".")
       dimnames(tmp)[[2]][2] <- paste(locus.names[i],dimnames(tmp)[[2]][2],sep=".")
       loci <- cbind(loci,tmp)
       unique.alleles[[i]] <- ualleles
    }

    if(exists("is.R") && is.function(is.R) && is.R()) {
      class(loci) <- "model.matrix"
    } else {
      oldClass(loci) <- "model.matrix"
    }

    attr(loci,"locus.names") <- locus.names
    attr(loci,"map") <- map
    attr(loci,"x.linked") <- x.linked
    attr(loci,"unique.alleles") <- unique.alleles
    attr(loci,"male.code") <- male.code
    attr(loci,"female.code") <- female.code
    attr(loci,"chrom.label") <- chrom.label

  
    return(loci)
}



