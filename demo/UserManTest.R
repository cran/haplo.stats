#$Author: sinnwell $
#$Date: 2004/03/10 21:29:16 $
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/UserManTest.q,v 1.24 2004/03/10 21:29:16 sinnwell Exp $
#$Locker:  $
#$Log: UserManTest.q,v $
#Revision 1.24  2004/03/10 21:29:16  sinnwell
#remove two parms of setupGeno
#
#Revision 1.23  2004/03/09 23:18:36  sinnwell
#use setupGeno, which is now loci, add comments for haplo.glm
#
#Revision 1.22  2004/03/04 22:56:37  sinnwell
#add use of setupGeno, but miss.val to be added to it
#
#Revision 1.21  2004/03/02 17:23:54  sinnwell
#use setupData and fix T/F's to full words
#
#Revision 1.20  2004/02/18 20:17:31  sinnwell
#remove extra steps for R use, demo & data now work
#
#Revision 1.19  2004/02/02 17:38:48  sinnwell
#insert comment for haplo.group appearing in different order than .pdf file
#
#Revision 1.18  2004/01/12 22:26:38  sinnwell
#put fit.inter after fit.gaus, as in UserMan.pdf
#
#Revision 1.17  2004/01/12 21:01:20  sinnwell
#add comments for assigning model.matrix by class or oldClass
#
#Revision 1.16  2003/11/14 19:47:42  sinnwell
#check is.R to assign class/oldClass of model.matrix to geno
#
#Revision 1.15  2003/10/27 15:27:42  schaid
#added haplo.glm examples
#
#Revision 1.14  2003/10/24 18:48:19  sinnwell
#update plot of score.slide.gaus call
#
#Revision 1.13  2003/10/20 19:52:44  schaid
#added more examples, including haplo.glm
#
#Revision 1.12  2003/10/02 20:03:10  sinnwell
#change haplo.score package name to haplo.stats
#
#Revision 1.11  2003/10/02 19:33:58  sinnwell
#update for hla.demo has 11 loci, and use haplo.score.slide, new simulations
#
#Revision 1.10  2003/03/26 16:23:42  sinnwell
#*** empty log message ***
#
#Revision 1.9  2003/03/20 22:43:01  sinnwell
#as in v1.6, make comments to source data when in R
#
#Revision 1.8  2003/03/19 20:12:13  sinnwell
#fix library() line
#
#Revision 1.7  2003/03/17 18:26:49  sinnwell
#take out R option to load hla.demo.  Now done in zzz.R
#
#Revision 1.6  2003/03/14 18:19:17  sinnwell
#make usable for both R and S, following instructions
#
#Revision 1.5  2003/03/10 20:31:32  sinnwell
#edit comments and install path
#
#Revision 1.4  2003/03/10 19:49:25  sinnwell
#update for changes made with R package
#
#Revision 1.3  2003/03/06 23:01:30  sinnwell
#add license text
#
#Revision 1.2  2003/01/30 22:23:16  sinnwell
#update for haplo.score v1.1
#
#Revision 1.1  2003/01/17 16:29:47  sinnwell
#Initial revision
#
# Test Splus code, and steps used to create User Manual

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

# If library is not installed, attach haplo.stats routines and hla.demo,
# by use of a local library path.  BE SURE TO CHANGE "INSTALL DIR PATH HERE" TO
# THE DIRECTORY PATH WHERE haplo.stats IS INSTALLED.

#library(haplo.stats, lib.loc="/INSTALL DIR PATH HERE/")

setupData(hla.demo)

names(hla.demo)
attach(hla.demo)
# Create Dataframe of columns of alleles for pairs of loci:
geno <- hla.demo[,c(17,18,21:24)]
label <-c("DQB","DRB","B")

geno.desc <- summaryGeno(geno, miss.val=c(0,NA))
print(geno.desc)

# Simulations are used in several of the functions (to determine random 
# starting values for haplo.em, and to compute permutation p-values in 
# haplo.score). In order to reproduce results in this user guide, 
# set the .Random.seed before any function which uses random numbers.
# In practice, the user would not ordinarily reset the seed.


seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

# Demo of haplo.em
set.seed(seed)
save.em <- haplo.em(geno=geno, locus.label=label, miss.val=c(0,NA))

# demo how to exclude subjects with any missing genotypes, when missing
# values are either NA ro 0:

keep <- !apply(is.na(geno) | geno==0, 1, any)
set.seed(seed)
save.em2 <- haplo.em(geno=geno[keep,], locus.label=label)


# print and summary of haplo.em

print(save.em)

summary(save.em)

summary(save.em, show.haplo=TRUE)

### haplo.group note ###
# At this point in the HaploStatsUserMan.pdf file we demo haplo.group
# The code for this is later with the demo of haplo.score.merge
######

# Haplotype analyses for a quantitative trait, resp with gaussian errors
set.seed(seed)
score.gaus <- haplo.score(resp, geno, trait.type="gaussian",
                          skip.haplo=.005, locus.label=label, simulate=FALSE)

print(score.gaus)

#plot(score.gaus)
#title("Figure 1. Haplotype Score Statistics\nQuantitative Response",cex=.6)
#loc.gaus <- locator.haplo(score.gaus)


# now create an eps file:

#postscript(file="Fig1.1.eps",onefile=TRUE, width=5, height=4, print.it=FALSE, horizo=FALSE)
#plot(score.gaus,cex=.6)
#title("Figure 1. Haplotype Score Statistics\nQuantitative Response",cex=.6)
#text(loc.gaus$x.coord,loc.gaus$y.coord,loc.gaus$hap.txt,cex=.5)
#dev.off()

# analysis with skip.haplo = .01, will keep only haplotypes with frequency
# above skip.haplo

set.seed(seed)
score.gaus.01 <- haplo.score(resp, geno, trait.type="gaussian",
                         offset = NA, x.adj = NA, skip.haplo=.01,
                         locus.label=label, miss.val=0, simulate=FALSE)
                      
print(score.gaus.01)


# Now try adjusted scores, adjusted for sex and age.

x.ma <- cbind(male, age) # set-up design matrix, with male=1/0 (male/female) and age in months
set.seed(seed)
score.gaus.adj <- haplo.score(resp, geno, trait.type="gaussian",
                        offset = NA, x.adj = x.ma, skip.haplo=.005,
                        locus.label=label, miss.val=0, simulate=FALSE)

print(score.gaus.adj)

# Haplotype analyses for ordinal trait, based on three categories for
# the quantitative response

y.ord <- as.numeric(resp.cat)
set.seed(seed)
score.ord <- haplo.score(y.ord, geno, trait.type="ordinal",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=label, miss.val=0, simulate=FALSE)

print(score.ord)

#motif()
#plot(score.ord)
#title("Haplotype Score Statistics\n Ordinal Response")
#loc.ord <- locator.haplo(score.ord)
#dev.off()

# now create an eps file:

#postscript(file="Fig1.2.eps",onefile=TRUE, width=5, height=4, print.it=FALSE, horizo=FALSE)
#plot(score.ord,cex=.6)
#title("Figure 1.2. Haplotype Score Statistics\nOrdinal Response",cex=.6)
#text(loc.ord$x.coord,loc.ord$y.coord,loc.ord$hap.txt,cex=.5)
#dev.off()


# Haplotype analyses for binary trait, y.bin
# 'high' and 'normal' code to 0, 1 is 'low'
y.bin <-ifelse(y.ord==1,1,0)

set.seed(seed)
score.bin <- haplo.score(y.bin, geno, trait.type="binomial",
                         offset = NA, x.adj = NA, skip.haplo=.005,
                         locus.label=label, miss.val=0, simulate=FALSE)
                      
print(score.bin)

#motif()
#plot(score.bin)
#title("Haplotype Score Statistics\nBinomial Response")
#loc.bin <- locator.haplo(score.bin)
#dev.off()

# create an eps file:

#postscript(file="Fig1.3.eps",onefile=TRUE, width=5, height=4, print.it=FALSE, horizo=FALSE)
#plot(score.bin,cex=.6)
#title("Figure 1.3. Haplotype Score Statistics\nBinary Response",cex=.6)
#text(loc.bin$x.coord,loc.bin$y.coord,loc.bin$hap.txt,cex=.5)
#dev.off()


# simulation p-values for binary trait
set.seed(seed)
score.bin.sim <- haplo.score(y.bin, geno, trait.type="binomial",
                        offset = NA, x.adj = NA, skip.haplo=.005,
                        locus.label=label, miss.val=0, simulate=TRUE,
                        sim.control=score.sim.control())
                       
print(score.bin.sim)


# Find estimated haplotype frequencies for both groups
# of the binary response.  Results stored in a haplo.group object
set.seed(seed)
group.bin <- haplo.group(y.bin, geno, locus.label=label, miss.val=0)

# group.bin is a large object, one row for all ~170 haplotypes.
# for a subset (recommended), use print(group.bin$group.df[1:50,]) or other row numbers
print(group.bin)

# merge haplo.group output with haplo.score output

merge.bin <- haplo.score.merge(score.bin, group.bin)

# similar to group.bin, merge.bin is very large.
# by default, the print method only prints those haplotypes with a score statistic
print(merge.bin)


# Slide haplo.score across all 11 loci to find group with most association with
# quantitative trait, resp

# set up the data frame with all 11 loci and the labels for them
geno.11 <- hla.demo[,-c(1:4)]
label.11 <- c("DPB","DPA","DMA","DMB","TAP1","TAP2","DQB","DQA","DRB","B","A")

set.seed(seed)
score.slide.gaus <- haplo.score.slide(resp, geno.11, trait.type="gaussian",
                          n.slide=3, skip.haplo=.005, locus.label=label.11)
                          
print(score.slide.gaus)


# plot the -log10 p-values against the 3 loci for which it was computed

#motif()
#plot(score.slide.gaus,cex=.6)
#title("Figure 2. Global p-values for\n sub-haplotypes; Gaussian Response",cex=.6)
#dev.off()

# create an eps file:
#postscript(file="Fig2.eps",onefile=TRUE, width=5, height=4, print.it=FALSE, horizo=FALSE)
#plot(score.slide.gaus,cex=.6)
#title("Figure 2. Global p-values for\n sub-haplotypes; Gaussian Response",cex=.6)
#dev.off()


# Demo of haplo.glm

label <-c("DQB","DRB","B")

y <- hla.demo$resp
y.bin <- 1*(hla.demo$resp.cat=="low")

# set up a genotype array as a model.matrix for inserting into data frame
# Note that hla.demo is a data.frame, and we need to subset to columns of interest.
# Also also need to convert to a matrix object, so that setupGeno can code alleles
# and convert geno to 'model.matrix' class.

geno <- as.matrix(hla.demo[,c(17,18,21:24)])

geno <- setupGeno(geno, miss.val=c(0,NA))

  # geno now has an attribute 'unique.alleles' which must be passed to
  # haplo.glm as allele.lev=attributes(geno)$unique.alleles, see below

my.data <- data.frame(geno=geno, age=hla.demo$age, male=hla.demo$male,
                      y=y, y.bin=y.bin)

seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

fit.gaus <- haplo.glm(y ~ male + geno, family = gaussian,  na.action="na.geno.keep",
                 allele.lev=attributes(geno)$unique.alleles, data=my.data,
                 locus.label=label,control = haplo.glm.control(haplo.freq.min=0.02))
fit.gaus

# Splus users may execute haplo.glm as above, but the allele.lev parameter
# is not required if geno and my.data are set up as:
# > oldClass(geno) <- "model.matrix"

# > my.data <- data.frame(geno=tmp$geno, age=hla.demo$age, male=hla.demo$male,
# >                       y=y, y.bin=y.bin)

# Since the first usage works in both Splus and R, we continue to assign a
# 'allele.lev' in the haplo.glm calls so the script runs in R & Splus

set.seed(seed)
fit.inter <- haplo.glm(y ~ male * geno, family = gaussian,  na.action="na.geno.keep",
                       allele.lev=attributes(geno)$unique.alleles, data=my.data, locus.label=label,
                       control = haplo.glm.control(haplo.freq.min = 0.02))
fit.inter

set.seed(seed)
fit.bin <- haplo.glm(y.bin ~ male + geno, family = binomial,  na.action="na.geno.keep",
                     allele.lev=attributes(geno)$unique.alleles, data=my.data,
                     locus.label=label, x=TRUE,
                     control = haplo.glm.control(haplo.freq.min = 0.02))
fit.bin

set.seed(seed)
fit.dom <- haplo.glm(y ~ male + geno, family = gaussian,  na.action="na.geno.keep",
                     allele.lev=attributes(geno)$unique.alleles, data=my.data, locus.label=label,
                     control = haplo.glm.control(haplo.freq.min = 0.02, 
                                                 haplo.effect="dom"))
fit.dom


