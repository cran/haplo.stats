#$Author: sinnwell $
#$Date: 2004/03/10 21:21:16 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/setupGeno.q,v 1.3 2004/03/10 21:21:16 sinnwell Exp $
#$Locker:  $
#$Log: setupGeno.q,v $
#Revision 1.3  2004/03/10 21:21:16  sinnwell
#remove excess parms and attributes
#
#Revision 1.2  2004/03/10 15:45:56  sinnwell
#make the same as loci function instead
#

# original version had some similar contents to this function.
# This is copying much code from loci(), but it is better named
# for users as setupGeno

setupGeno <- function(geno, miss.val = c(0,NA))
{
	# Title: Create an object of class model.matrix
	loci <- NULL
	unique.alleles <- list(NULL)
	if(missing(geno))
		stop("Error: Required argument geno is missing")
	if(!is.data.frame(geno) & !is.matrix(geno))
		stop("Error: Argument geno is not a data.frame or matrix")
	if(ncol(geno) %% 2)
		stop("Error: Number of columns in geno object is not an even number")
	exp.num.loc <- ncol(geno)/2

	for(i in 1:exp.num.loc) {
          locus.names <- paste("loc",1:exp.num.loc, sep="-")
          tmp <- locus(geno[, (i * 2 - 1)], geno[, i * 2], chrom.label = NULL,
                       locus.names[i], x.linked = FALSE, sex = NULL,
                       male.code = "M", female.code = "F", miss.val = miss.val)
          ualleles <- attr(tmp, "allele.labels")
          tmp <- as.matrix(tmp)
          dimnames(tmp)[[2]][1] <- paste(locus.names[i], dimnames(tmp)[[2]][1], sep = ".")
          dimnames(tmp)[[2]][2] <- paste(locus.names[i], dimnames(tmp)[[2]][2], sep = ".")
          loci <- cbind(loci, tmp)
          unique.alleles[[i]] <- ualleles
	}
	if(exists("is.R") && is.function(is.R) && is.R()) {
          class(loci) <- "model.matrix"
	}
	else {
          oldClass(loci) <- "model.matrix"
	}

	attr(loci, "unique.alleles") <- unique.alleles

        return(loci)
}
