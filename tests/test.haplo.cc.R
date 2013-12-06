#$Author: sinnwell $
#$Date: 2011/12/05 20:53:27 $
#$Header: /projects/genetics/cvs/cvsroot/haplo.stats/test/test.haplo.cc.R,v 1.1 2011/12/05 20:53:27 sinnwell Exp $
#$Locker:  $
#$Log: test.haplo.cc.R,v $
#Revision 1.1  2011/12/05 20:53:27  sinnwell
#changed from .q to .R, to work with R check
#
#
#library(genetics) or library(haplo.stats) for R

## package: haplo.stats
## test script: haplo.cc

## settings
verbose=TRUE
require(haplo.stats)
Sys.setlocale("LC_ALL", "C")
Sys.getlocale()


# Jason Sinnwell 3/2004
# Mayo Clinic, Biostatistics

  if(verbose) cat("setting up data...\n")
  
  label <- c("DQB","DRB","B")

  if(exists("is.R") && is.function(is.R) && is.R()) data(hla.demo)
  y.bin <- 1*(hla.demo$resp.cat=="low")

  geno <- as.matrix(hla.demo[,c(17,18,21:24)])

  # commented code was to check data that goes into haplo.cc
  # and gets pasted together in the huge data.frame.
  
 # sink(file="results.haplo.cc.txt")  

  seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
  set.seed(seed)

#  group.hla <- haplo.group(y.bin, geno, miss.val=0, locus.label=label,
#                         control=haplo.em.control(n.try=1))
#  set.seed(seed)
  if(verbose) cat("hla data \n")
  cc.hla <- haplo.cc(y.bin, geno, miss.val=0,locus.label=label, 
                     control=haplo.glm.control(haplo.min.count=8,
                       em.c=haplo.em.control()))
  
#merge.hla <- haplo.score.merge(cc.hla$score.lst, group.hla)
#options(width=120)
#cat("\nHaplo.score output\n")
#print(cc.hla$score.lst)
#cat("\nhaplo.group output\n")
#print(group.hla)
#cat("\nhaplo.score.merge output\n")
#print(merge.hla, order.by="haplotype", digits=4)
#cat("\nhaplo.glm output\n")
#print(cc.hla$fit.lst)
#cat("\nThe whole thing, haplo.cc output\n")
#print(cc.hla, order.by="haplotype", digits=4)
   
  "geno.test" <- 
  matrix(c(1., 3., 1., 2., 1., 1., 2., 2., 2., 1., 3., 2., 1., 2., 2., 1., 1., 2., 3.,
        2., 1., 2., 1., 1., 1., 1., 1., 1., 2., 2., 2., 2., 2., 1., 1., 2.,
        1., 2., 3., 2., 2., 1., 2., 1., 1., 2., 2., 2., 1., 1., 2., 2., 1.,
        2., 2., 1., 1., 1., 2., 1., 1., 1., 2., 2., 1., 2., 2., 2., 1., 1.,
        2., 1., 2., 2., 1., 1., 2., 2., 2., 1., 2., 3., 2., 2., 3., 3., 1.,
        4., 3., 3., 2., 3., 3., 3., 2., 3., 2., 4., 2., 2., 3., 4., 3., 4.,
        2., 2., 3., 4., 2., 2., 4., 3., 1., 4., 2., 3., 3., 3., 3., 3.)
  , nrow = 20, ncol = 6)

  "y.test" <-  c(0., 0., 1., 1., 1., 1.,
               1., 0., 0., 1., 1., 1., 0., 0., 1., 1., 1., 0., 1.,1.)

  locus.label <- c("A", "B", "C")

  set.seed(seed)

  if(verbose) cat("small numeric data... \n")
  cc.test <- haplo.cc(y.test, geno.test, locus.label=locus.label,
                      ci.prob=.95, control=haplo.glm.control(haplo.min.count=2))

  geno.char <- ifelse(geno.test==1, 'A',ifelse(geno.test==2, 'T',
                        ifelse(geno.test==3, 'G', 'C')))

  set.seed(seed)
  if(verbose) cat("small char data with simulations... \n")
  cc.char.sim <- haplo.cc(y.test, geno.char, locus.label=locus.label, 
                          ci.prob=.90, simulate=TRUE,
                          control = haplo.glm.control(haplo.min.count=2))
  

  print.haplo.cc(cc.hla, digits=2, nlines=40)
  print.haplo.cc(cc.test, order.by="score", digits=2)
  print(cc.test, order.by="freq", digits=2)
  print(cc.test$fit.lst, digits=2)
  print.haplo.cc(cc.char.sim, order.by='haplotype', digits=2)
  print(cc.char.sim$score.lst, digits=2)
    
