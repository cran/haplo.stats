#$Author: sinnwell $
#$Date: 2004/03/01 20:52:37 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo.em.q,v 1.11 2004/03/01 20:52:37 sinnwell Exp $
#$Locker:  $
#$Log: haplo.em.q,v $
#Revision 1.11  2004/03/01 20:52:37  sinnwell
#change T to TRUE for matrix()
#
#Revision 1.10  2004/02/02 23:00:01  sinnwell
#insert ghost runif(1) line to init .Random.seed for R bug
#
#Revision 1.9  2003/09/19 21:38:45  schaid
#fixed returned class to be R compatible
#
#Revision 1.8  2003/08/26 22:08:31  sinnwell
#add GPL License
#
#Revision 1.7  2003/08/26 21:04:49  schaid
#Major revision of haplo.em, by adding new functions and new C code to use progressive insertion of loci.
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
haplo.em  <- function(geno, locus.label=NA, miss.val=c(0,NA), weight=NULL, 
                          control = haplo.em.control()  ){



n.loci <- ncol(geno)/2
n.subject <- nrow(geno)
subj.id <- 1:n.subject

# set up weight
if(any(is.null(weight))){
  weight <- rep(1,n.subject)
}

if(any(weight<0)){
  stop("negative weights not allowed")
}

if(length(weight)!=n.subject){
  stop("Length of weight != number of subjects (nrow of geno)")
}

# Create locus label if not included
if(all(is.na(locus.label))) locus.label<- paste("loc-",1:n.loci,sep="")

if(length(locus.label)!=n.loci){
  stop("length of locus.label != n.loci")
}


# recode geno to integer values, accounting for missing values
temp.geno <- loci(geno,locus.names=locus.label,miss.val=miss.val)


# Compute the max number of pairs of haplotypes over all subjects
max.pairs <- geno.count.pairs(temp.geno)
max.haps <- 2*sum(max.pairs)

if(max.haps > control$max.haps.limit ) max.haps <- control$max.haps.limit


# check whether to delete some rows - now defunct, but use 
# dummy to not break code that uses this in returned list
rows.rem <- numeric(0) 
rows.keep <- subj.id

geno.vec <- as.vector(temp.geno)
geno.vec <- ifelse(is.na(geno.vec),0,geno.vec)

allele.labels <- attr(temp.geno, "unique.alleles")


if(length(allele.labels)!=n.loci)
  stop("Number of loci in alleles list != n.loci")

n.alleles <- numeric(n.loci)
a.freq <- vector("list",n.loci)

for(i in 1:n.loci){
  n.alleles[i] <- length(allele.labels[[i]])
  j <- (i-1)*2 + 1
  p <- table(temp.geno[,c(j, (j+1))], exclude=NA)
  p <- p/sum(p)
  a.freq[[i]] <- list(p=p)
}

if(is.null(control$loci.insert.order)) {
  control$loci.insert.order <- 1:n.loci
}


# need zero-offset for loci-insert-order when
# pass to C

loci.insert.order <- (control$loci.insert.order - 1)


if(length(loci.insert.order) != n.loci){
  stop("length of loci.insert.order != n.loci")
}

if(sum( abs(sort(loci.insert.order) - (0:(n.loci-1)))) > 0){
  stop("All loci are not accounted for in  loci.insert.order")
}

if(control$insert.batch.size > n.loci){
  control$insert.batch.size <- n.loci
}


if(!is.null(control$iseed)) {
  set.seed(control$iseed) } else
{  runif(1)
   control$iseed <- .Random.seed
 }


# The seeds for the ranAS183 random number generator used in the C function
# hwe_sim must be between 1 and 30000, but bigger is better (we think), so we
# add 10000

seed.array <- runif(3)

iseed1 = 10000 + 20000*seed.array[1]
iseed2 = 10000 + 20000*seed.array[2]
iseed3 = 10000 + 20000*seed.array[3]


fit <- haplo.em.fitter(
		       n.loci,
		       n.subject,
                       weight,
     		       geno.vec,
                       n.alleles,
                       max.haps,
		       max.iter=control$max.iter,
		       loci.insert.order,		      
		       min.posterior=control$min.posterior,
		       tol=control$tol,
		       insert.batch.size=control$insert.batch.size,
                       random.start=control$random.start,            
                       iseed1=iseed1,          
                       iseed2=iseed2,
                       iseed3=iseed3,
                       verbose=control$verbose)


  # if n.try > 1 try, remaining tries are random starts for posteriors
  if(control$n.try > 1){
    for(i in 2:control$n.try){
  
     seed.array <- runif(3)

     iseed1 = 10000 + 20000*seed.array[1]
     iseed2 = 10000 + 20000*seed.array[2]
     iseed3 = 10000 + 20000*seed.array[3]

     fit.new <- haplo.em.fitter(
		       n.loci,
		       n.subject,
                       weight,
     		       geno.vec,
                       n.alleles,
                       max.haps,
		       max.iter=control$max.iter,
		       loci.insert.order,		      
		       min.posterior=control$min.posterior,
		       tol=control$tol,
		       insert.batch.size=control$insert.batch.size,
                       random.start=1,            
                       iseed1=iseed1,          
                       iseed2=iseed2,
                       iseed3=iseed3,
                       verbose=control$verbose)

      if(fit.new$tmp1$lnlike > fit$tmp1$lnlike)
         { fit <- fit.new
         }
    } 
  }



tmp1 <- fit$tmp1
tmp2 <- fit$tmp2

u.hap <- matrix(tmp2$u.hap,nrow=tmp2$n.u.hap,byrow=TRUE)

# code alleles for haplotpes with original labels

haplotype <- NULL
for(j in 1:n.loci){
  haplotype <- cbind(haplotype, allele.labels[[j]][u.hap[,j]])
}

haplotype <- data.frame(haplotype)
names(haplotype) <- locus.label


# convert from 0-offset in C to 1-offset in S, and recode hap codes
# to 1,2,..., n_uhap

hap1code  <- tmp2$hap1code  + 1
hap2code  <- tmp2$hap2code  + 1
uhapcode  <- tmp2$u.hap.code + 1

n1 <- length(uhapcode)
n2 <- length(hap1code)

tmp <- as.numeric(factor(c(uhapcode, hap1code, hap2code)))
uhapcode <- tmp[1:n1]
hap1code <- tmp[(n1+1):(n1+n2)]
hap2code <- tmp[(n1+n2+1):(n1+2*n2)]

uhap.df <- data.frame(uhapcode, tmp2$hap.prob, u.hap)

names(uhap.df) <- c("hap.code","hap.prob",locus.label)

indx.subj = tmp2$indx.subj + 1
subj.used.id <-  subj.id[rows.keep][indx.subj]


# compute lnlike if no LD. This is a rough approximation which will
# be accurate only if all haplotypes are considered in the list
# of enumerated haplotypes. If there is no trimming
# (min.posterior = 0), and it is possible to enumerate all
# possible pairs of haplotypes, then hap.prob.noLD will sum
# to 1. But, with trimming, this may not occur, so we 
# rescale to force them to sum to 1. This may not lead to
# an accurate test for no LD by the likelihood ratio statistic.

hap.prob.noLD <- a.freq[[1]]$p[u.hap[,1]]
df.noLD <- length(a.freq[[1]]$p) - 1

for(j in 2:n.loci){
    hap.prob.noLD <- hap.prob.noLD *  a.freq[[j]]$p[u.hap[,j]]
    df.noLD <- df.noLD + length(a.freq[[j]]$p) - 1
  }
hap.prob.noLD <-  hap.prob.noLD/sum(hap.prob.noLD)

prior.noLD <- hap.prob.noLD[hap1code]*hap.prob.noLD[hap2code]
prior.noLD <- ifelse(hap1code!=hap2code, 2*prior.noLD, prior.noLD)
ppheno.noLD <- tapply(prior.noLD, indx.subj, sum)
lnlike.noLD <- sum(log(ppheno.noLD))

lr = 2*(tmp1$lnlike - lnlike.noLD)
df.LD <- sum(tmp2$hap.prob > 0.0000001) - 1
df.lr <-  df.LD - df.noLD


obj <- list(
  lnlike=tmp1$lnlike,
  lr = lr,
  df.lr = df.lr,
  hap.prob = tmp2$hap.prob,
  hap.prob.noLD = hap.prob.noLD,
  converge = tmp1$converge,
  locus.label = locus.label,
  indx.subj = indx.subj,    
  subj.id = subj.used.id, 
  post = tmp2$post,
  hap1code = hap1code,
  hap2code = hap2code,
  haplotype = haplotype,
  nreps = table(indx.subj),
  rows.rem = rows.rem,
  max.pairs=max.pairs,
  control=control)


   if(exists("is.R") && is.function(is.R) && is.R()) {
     class(obj) <- "haplo.em"
   } else {
     oldClass(obj) <- "haplo.em"
   }

return(obj)


}
