###################################################
### chunk number 1: 
###################################################
# load the library, load and preview at demo dataset
library(haplo.stats)


###################################################
### chunk number 2: 
###################################################
options(width=60)
cat("## if local library, use: \n\t ## library(haplo.stats, lib.loc=\"/local/install/path/\")\n")


###################################################
### chunk number 3:  eval=FALSE
###################################################
## help(haplo.em)


###################################################
### chunk number 4: 
###################################################
# load and preview demo dataset stored in ~/haplo.stats/data/hla.demo.tab
data(hla.demo)
names(hla.demo)
# attach hla.demo to make columns available in the session
attach(hla.demo)



###################################################
### chunk number 5: 
###################################################
geno <- hla.demo[,c(17,18,21:24)]
label <-c("DQB","DRB","B")



###################################################
### chunk number 6: 
###################################################
geno.desc <- summaryGeno(geno, miss.val=c(0,NA))
print(geno.desc[c(1:10,80:85,135:140),])



###################################################
### chunk number 7: 
###################################################

# find if there are any people missing all alleles
table(geno.desc[,3])



###################################################
### chunk number 8:  eval=FALSE
###################################################
## # create an index of people missing all alleles
## miss.all <- which(geno.desc[,3]==3)
## 
## # use index to subset hla.demo
## hla.demo.updated <- hla.demo[-miss.all,]
## 


###################################################
### chunk number 9: 
###################################################
# this is how to set the seed for reproducing results where haplo.em is 
# involved, and also if simulations are run. In practice, don't reset seed.
seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)


###################################################
### chunk number 10: seed
###################################################
set.seed(seed)


###################################################
### chunk number 11: 
###################################################
save.em <- haplo.em(geno=geno, locus.label=label, miss.val=c(0,NA))
print(save.em, nlines=10)


###################################################
### chunk number 12: 
###################################################
summary(save.em, nlines=7)


###################################################
### chunk number 13: 
###################################################
# show full haplotypes, instead of codes
summary(save.em, show.haplo=TRUE, nlines=7)


###################################################
### chunk number 14:  eval=FALSE
###################################################
## # demonstrate only the syntax of control parameters
## save.em <- haplo.em(geno=geno, locus.label=label, miss.val=c(0, NA),
##      control = haplo.em.control(n.try = 20, insert.batch.size=2))


###################################################
### chunk number 15: seed
###################################################
set.seed(seed)


###################################################
### chunk number 16: 
###################################################
## run haplo.em on sub-groups
## create ordinal and binary variables
y.bin <- 1*(resp.cat=="low")
group.bin <- haplo.group(y.bin, geno, locus.label=label, miss.val=0)
print(group.bin, nlines=15)


###################################################
### chunk number 17: seed
###################################################
set.seed(seed)


###################################################
### chunk number 18: 
###################################################
# score statistics w/ Gaussian trait
score.gaus.add <- haplo.score(resp, geno, trait.type="gaussian",
               haplo.effect="additive", min.count=5, 
               locus.label=label, simulate=FALSE)
print(score.gaus.add, nlines=10)


###################################################
### chunk number 19: seed
###################################################
set.seed(seed)


###################################################
### chunk number 20: 
###################################################
# scores w/ ordinal trait
y.ord <- as.numeric(resp.cat)
score.ord <- haplo.score(y.ord, geno, trait.type="ordinal",
                     x.adj = NA, min.count=5,
                     haplo.effect="additive", locus.label=label, 
                     miss.val=0, simulate=FALSE)
print(score.ord, nlines=7)


###################################################
### chunk number 21: seed
###################################################
set.seed(seed)


###################################################
### chunk number 22: 
###################################################
# scores, binary trait
y.bin <- 1*(resp.cat=="low")
score.bin <- haplo.score(y.bin, geno, trait.type="binomial",
                        x.adj = NA, min.count=5,
                        haplo.effect="additive", locus.label=label, 
                        miss.val=0, simulate=FALSE)
print(score.bin, nlines=10)


###################################################
### chunk number 23: 
###################################################
## plot score vs. frequency, gaussian response
plot(score.gaus.add)
 
## locate and label pts with their haplotypes
## works similar to locator() function
#> pts.haplo <- locator.haplo(score.gaus)

pts.haplo <- list(x.coord=c(0.05098, 0.03018, .100), 
                  y.coord=c(2.1582, 0.45725, -2.1566), 
                  hap.txt=c("62:2:7", "51:1:35", "21:3:8"))

text(x=pts.haplo$x.coord, y=pts.haplo$y.coord, labels=pts.haplo$hap.txt)


###################################################
### chunk number 24: seed
###################################################
set.seed(seed)


###################################################
### chunk number 25: 
###################################################
# increase skip.haplo, expected hap counts = 10  
score.gaus.min10 <- haplo.score(resp, geno, trait.type="gaussian",
                          x.adj = NA, min.count=10,
                          locus.label=label, miss.val=0, simulate=FALSE)
print(score.gaus.min10)


###################################################
### chunk number 26: seed
###################################################
set.seed(seed)


###################################################
### chunk number 27: 
###################################################
# score w/gaussian, adjusted by covariates
x.ma <- cbind(male, age)
score.gaus.adj <- haplo.score(resp, geno, trait.type="gaussian",
                        x.adj = x.ma, min.count=5,
                        locus.label=label, miss.val=0, simulate=FALSE)
print(score.gaus.adj, nlines=10)


###################################################
### chunk number 28: seed
###################################################
set.seed(seed)


###################################################
### chunk number 29: 
###################################################
# score w/gaussian, dominant effect

score.gaus.dom <- haplo.score(resp, geno, trait.type="gaussian",
                        x.adj=NA, min.count=5,
                        haplo.effect="dominant", locus.label=label, 
                        miss.val=0, simulate=FALSE)
print(score.gaus.dom, nlines=10)


###################################################
### chunk number 30: seed
###################################################
set.seed(seed)


###################################################
### chunk number 31: 
###################################################
# simulations when binary response
score.bin.sim <- haplo.score(y.bin, geno, trait.type="binomial",
              x.adj = NA, locus.label=label, min.count=5,
              haplo.effect="additive", miss.val=0, simulate=TRUE, 
              sim.control = score.sim.control() )
print(score.bin.sim)


###################################################
### chunk number 32: 
###################################################
# set up data for haplo.glm, include geno.glm, 
# covariates age and male, and responses resp and y.bin
geno <- hla.demo[,c(17,18,21:24)]
geno.glm <- setupGeno(geno, miss.val=c(0,NA), locus.label=label)
y.bin <- 1*(resp.cat=="low")
my.data <- data.frame(geno.glm, age=age, male=male, y=resp, y.bin=y.bin)


###################################################
### chunk number 33: seed
###################################################
set.seed(seed)


###################################################
### chunk number 34: 
###################################################
# glm fit with haplotypes, additive gender covariate on gaussian response
fit.gaus <- haplo.glm(y ~ male + geno.glm, family=gaussian, data=my.data, 
        na.action="na.geno.keep", locus.label = label, 
        allele.lev = attributes(geno.glm)$unique.alleles,
        control=haplo.glm.control(haplo.min.count=5))
print(fit.gaus)



###################################################
### chunk number 35: seed
###################################################
set.seed(seed)


###################################################
### chunk number 36: 
###################################################
# glm fit haplotypes with covariate interaction
fit.inter <- haplo.glm(formula = y ~ male * geno.glm, 
                   family = gaussian, data=my.data, 
                   na.action="na.geno.keep", 
                   locus.label = label, allele.lev = 
                   attributes(geno.glm)$unique.alleles, 
                   control = haplo.glm.control(haplo.min.count = 10))
print(fit.inter)


###################################################
### chunk number 37: seed
###################################################
set.seed(seed)


###################################################
### chunk number 38: 
###################################################
# gender and haplotypes fit on binary response, 
# return model matrix
fit.bin <- haplo.glm(y.bin ~ male + geno.glm, family = binomial, 
             data=my.data, na.action = "na.geno.keep",
             locus.label=label,
             allele.lev=attributes(geno.glm)$unique.alleles,
             control = haplo.glm.control(haplo.min.count=10))
print(fit.bin)



###################################################
### chunk number 39: seed
###################################################
set.seed(seed)


###################################################
### chunk number 40: 
###################################################
# control dominant effect of haplotypes (haplo.effect) 
# by using haplo.glm.control
fit.dom <- haplo.glm(y ~ male + geno.glm, family = gaussian,
       data = my.data, na.action = "na.geno.keep", 
       locus.label = label, 
       allele.lev = attributes(geno.glm)$unique.alleles, 
       control = haplo.glm.control(haplo.effect='dominant', 
       haplo.min.count=8))
print(fit.dom)


###################################################
### chunk number 41: seed
###################################################
set.seed(seed)


###################################################
### chunk number 42: 
###################################################
# control baseline selection, perform the same exact run as fit.bin, 
# but different baseline by using haplo.base chosen from haplo.common
fit.bin$haplo.common
fit.bin$haplo.freq.init[fit.bin$haplo.common]
fit.bin.base77 <- haplo.glm(y.bin ~ male + geno.glm, family = binomial,
       data = my.data, na.action = "na.geno.keep",
       locus.label = label, 
       allele.lev = attributes(geno.glm)$unique.alleles, 
       control = haplo.glm.control(haplo.base=77, 
       haplo.min.count=8))
print(fit.bin.base77)


###################################################
### chunk number 43: 
###################################################
# merge haplo.score and haplo.group results
merge.bin <- haplo.score.merge(score.bin, group.bin)
print(merge.bin, nlines=10)


###################################################
### chunk number 44: seed
###################################################
set.seed(seed)


###################################################
### chunk number 45: 
###################################################

# demo haplo.cc where haplo.min.count is specified
# use geno, and this function prepares it for haplo.glm
y.bin <- 1*(hla.demo$resp.cat=="low")
cc.hla <- haplo.cc(y=y.bin, geno=geno, haplo.min.count=8,
                   locus.label = label,
                   control=haplo.glm.control(em.c=haplo.em.control(iseed=10)))
print(cc.hla, nlines=25, digits=2)
names(cc.hla)


###################################################
### chunk number 46: seed
###################################################
set.seed(seed)


###################################################
### chunk number 47: 
###################################################
# haplo.score on 11 loci, slide on 3 consecutive loci at a time
geno.11 <- hla.demo[,-c(1:4)]
label.11 <- c("DPB","DPA","DMA","DMB","TAP1","TAP2","DQB","DQA","DRB","B","A")
score.slide.gaus <- haplo.score.slide(hla.demo$resp, geno.11, trait.type =
                "gaussian", n.slide=3, min.count=5, locus.label=label.11)
print(score.slide.gaus)


###################################################
### chunk number 48: 
###################################################
# plot global p-values for sub-haplotypes from haplo.score.slide
plot(score.slide.gaus, las=2)


###################################################
### chunk number 49: seed
###################################################
set.seed(seed)


###################################################
### chunk number 50: 
###################################################
geno.11 <- hla.demo[,-c(1:4)]
y.bin <- 1*(hla.demo$resp.cat=="low")
hla.summary <- summaryGeno(geno.11, miss.val=c(0,NA))

# track those subjects with too many possible haplotype pairs ( > 50,000)
many.haps <- (1:length(y.bin))[hla.summary[,4]>50000]

# For speed, or even just so it will finish, make y.bin and geno.scan 
# for genotypes that don't have too many ambigous haplotypes
geno.scan <- geno.11[-many.haps,]
y.scan <- y.bin[-many.haps]

# scan haplotypes for regions within width of 3 for each locus.
# test statistic measures difference in haplotype counts in cases and controls
# p-values are simulated for each locus and the maximum statistic, 
# we do 100 simuations here, should use default settings for analysis

scan.hla <- haplo.scan(y.scan, geno.scan, width=3,
        sim.control=score.sim.control(min.sim=100, max.sim=100),
        em.control=haplo.em.control())

print(scan.hla)


###################################################
### chunk number 51: seed
###################################################
set.seed(seed)


###################################################
### chunk number 52: 
###################################################
# define binary response and genotype matrix
data(seqhap.dat)
data(seqhap.pos)
y <- seqhap.dat$disease
geno <- seqhap.dat[,-1]
# get vector with chrom position
pos <- seqhap.pos$pos
seqhap.out <- seqhap(y=y, geno=geno, pos=pos,
                     n.sim=1000, r2.threshold=.95, 
                     mh.threshold=3.84, miss.val=c(0,NA))

print(seqhap.out)



###################################################
### chunk number 53: 
###################################################
# plot global p-values for sub-haplotypes from haplo.score.slide
plot(seqhap.out, pval="hap", single=TRUE, las=2)


###################################################
### chunk number 54: 
###################################################
# create a matrix of haplotype effect columns from haplo.em result
hap.effect.frame <- haplo.design(save.em)

names(hap.effect.frame)

hap.effect.frame[1:10,1:8]



###################################################
### chunk number 55: 
###################################################
# create haplotype effect cols for haps 4 and 138
hap4.hap138.frame <- haplo.design(save.em, hapcodes=c(4,138), 
                                 haplo.effect="dominant")

hap4.hap138.frame[1:10,]

dat.glm <- data.frame(resp, male, age, 
                      hap.4=hap4.hap138.frame$hap.4, 
                      hap.138=hap4.hap138.frame$hap.138)

glm.hap4.hap138 <- glm(resp ~ male + age + hap.4 + hap.138, 
                       family="gaussian", data=dat.glm)
summary(glm.hap4.hap138)


