###################################################
### chunk number 1: 
###################################################
options(width=60)
cat("## if local library, use: \n\t ## library(haplo.stats, lib.loc=\"/local/install/path/\")\n")


###################################################
### chunk number 2: 
###################################################
# load the library, load and preview at demo dataset
library(haplo.stats)
setupData(hla.demo)
attach(hla.demo)
names(hla.demo)


###################################################
### chunk number 3: 
###################################################
geno <- hla.demo[,c(17,18,21:24)]
label <-c("DQB","DRB","B")


###################################################
### chunk number 4: 
###################################################
# this is how to set the seed for reproducing results where haplo.em is 
# involved, and also if simulations are run. In practice, don't reset seed.
seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)


###################################################
### chunk number 5: 
###################################################
geno.desc <- summaryGeno(geno, miss.val=c(0,NA))
print(geno.desc[c(1:10,80:85,135:140),])


###################################################
### chunk number 6: seed
###################################################
set.seed(seed)


###################################################
### chunk number 7: 
###################################################
save.em <- haplo.em(geno=geno, locus.label=label, miss.val=c(0,NA))
print(save.em, nlines=10)


###################################################
### chunk number 8: 
###################################################
# if you discard the geno rows with missings, haplo.em runs in a matter 
# of seconds as opposed to ~5 minutes
keep <- !apply(is.na(geno) | geno==0, 1, any)
set.seed(seed)
save.em2 <- haplo.em(geno=geno[keep,], locus.label=label)


###################################################
### chunk number 9: 
###################################################
summary(save.em, nlines=7)


###################################################
### chunk number 10: 
###################################################
# show full haplotypes, instead of codes
summary(save.em, show.haplo=TRUE, nlines=7)


###################################################
### chunk number 11: 
###################################################
# demonstrate only the syntax of control parameters
save.em <- haplo.em(geno=geno, locus.label=label, miss.val=c(0, NA),
     control = haplo.em.control(n.try = 20, insert.batch.size=2))


###################################################
### chunk number 12: seed
###################################################
set.seed(seed)


###################################################
### chunk number 13: 
###################################################
## run haplo.em on sub-groups
## create ordinal and binary variables
y.bin <- 1*(hla.demo$resp.cat=="low")
group.bin <- haplo.group(y.bin, geno, locus.label=label, miss.val=0)
print(group.bin, nlines=15)


###################################################
### chunk number 14: seed
###################################################
set.seed(seed)


###################################################
### chunk number 15: 
###################################################
# score statistics w/ Gaussian trait
score.gaus <- haplo.score(resp, geno, trait.type="gaussian",
               skip.haplo=5/(2*nrow(geno)), 
               locus.label=label, simulate=FALSE)
print(score.gaus, nlines=10)


###################################################
### chunk number 16: seed
###################################################
set.seed(seed)


###################################################
### chunk number 17: 
###################################################
# scores w/ ordinal trait
y.ord <- as.numeric(resp.cat)
score.ord <- haplo.score(y.ord, geno, trait.type="ordinal",
                     offset = NA, x.adj = NA, skip.haplo=5/(2*nrow(geno)),
                     locus.label=label, miss.val=0, simulate=FALSE)
print(score.ord, nlines=7)


###################################################
### chunk number 18: seed
###################################################
set.seed(seed)


###################################################
### chunk number 19: 
###################################################
# scores, binary trait
y.bin <- 1*(hla.demo$resp.cat=="low")
score.bin <- haplo.score(y.bin, geno, trait.type="binomial",
                        offset = NA, x.adj = NA, skip.haplo=5/(2*nrow(geno)),
                        locus.label=label, miss.val=0, simulate=FALSE)
print(score.bin, nlines=10)


###################################################
### chunk number 20: 
###################################################
## plot score vs. frequency, gaussian response
plot(score.gaus)

## locate and label pts with their haplotypes
## works similar to locator() function
# pts.haplo <- locator.haplo(score.gaus)

## and the results from locating three distinct pts can be 
## re-created by these two steps
cat("These next two steps substitute for doing: \n", 
 "\t > locator.haplo(score.gaus)\n")

pts.haplo <- list(x.coord=c(0.05098, 0.03018, .100), 
                  y.coord=c(2.1582, 0.45725, -2.1566), 
                  hap.txt=c("62:2:7", "51:1:35", "21:3:8"))

text(x=pts.haplo$x.coord, y=pts.haplo$y.coord, labels=pts.haplo$hap.txt)


###################################################
### chunk number 21: seed
###################################################
set.seed(seed)


###################################################
### chunk number 22: 
###################################################
# increase skip.haplo, expected hap counts = 2*220*(.02)=8.8
score.gaus.02 <- haplo.score(resp, geno, trait.type="gaussian",
                          offset = NA, x.adj = NA, skip.haplo=.02,
                          locus.label=label, miss.val=0, simulate=FALSE)
print(score.gaus.02)


###################################################
### chunk number 23: seed
###################################################
set.seed(seed)


###################################################
### chunk number 24: 
###################################################
# score w/gaussian, adjusted by covariates
x.ma <- cbind(male, age)
score.gaus.adj <- haplo.score(resp, geno, trait.type="gaussian",
                        offset = NA, x.adj = x.ma, skip.haplo=5/(2*nrow(geno)),
                        locus.label=label, miss.val=0, simulate=FALSE)
print(score.gaus.adj, nlines=10)


###################################################
### chunk number 25: seed
###################################################
set.seed(seed)


###################################################
### chunk number 26: 
###################################################
# simulations when binary response
score.bin.sim <- haplo.score(y.bin, geno, trait.type="binomial",
              offset = NA, x.adj = NA, locus.label=label, 
              miss.val=0, simulate=TRUE, 
              sim.control = score.sim.control() )
print(score.bin.sim)


###################################################
### chunk number 27: 
###################################################
# set up data for haplo.glm, include geno.glm, 
# covariates age and male, and responses resp and y.bin
geno <- as.matrix(hla.demo[,c(17,18,21:24)])
geno.glm <- setupGeno(geno, miss.val=c(0,NA))
y.bin <- 1*(hla.demo$resp.cat=="low")
my.data <- data.frame(geno.glm, age=age, male=male, y=resp, y.bin=y.bin)


###################################################
### chunk number 28: seed
###################################################
set.seed(seed)


###################################################
### chunk number 29: 
###################################################
# glm fit with haplotypes, additive gender covariate on gaussian response
fit.gaus <- haplo.glm(y ~ male + geno.glm, family=gaussian, data=my.data, 
        na.action="na.geno.keep", locus.label = label, 
        allele.lev = attributes(geno.glm)$unique.alleles,
        control=haplo.glm.control(haplo.min.count=5))
print(fit.gaus)



###################################################
### chunk number 30: seed
###################################################
set.seed(seed)


###################################################
### chunk number 31: 
###################################################
# glm fit haplotypes with covariate interaction
fit.inter <- haplo.glm(y ~ male * geno.glm, family = gaussian,
           data=my.data, na.action = "na.geno.keep", 
           locus.label = label, allele.lev = 
           attributes(geno.glm)$unique.alleles, 
           control = haplo.glm.control(haplo.min.count = 5))
print(fit.inter)


###################################################
### chunk number 32: seed
###################################################
set.seed(seed)


###################################################
### chunk number 33: 
###################################################
# gender and haplotypes fit on binary response, 
# return model matrix, x=TRUE
fit.bin <- haplo.glm(y.bin ~ male + geno.glm, family = binomial, 
             data=my.data, na.action = "na.geno.keep",
             locus.label= label,   #  x=TRUE,
             allele.lev=attributes(geno.glm)$unique.alleles,
             control = haplo.glm.control(haplo.min.count=5))
print(fit.bin)


###################################################
### chunk number 34: 
###################################################
# determine how many subject have homozygous copies of haplotypes
# hardly any, so a recessive model doesn't work
is.homzyg <- save.em$hap1code == save.em$hap2code
homzyg.counts <- table(save.em$hap1code[is.homzyg])
homzyg.counts


###################################################
### chunk number 35: seed
###################################################
set.seed(seed)


###################################################
### chunk number 36: 
###################################################
# control dominant effect of haplotypes (haplo.effect) 
# by using haplo.glm.control
fit.dom <- haplo.glm(y ~ male + geno.glm, family = gaussian,
       data = my.data, na.action = "na.geno.keep", 
       locus.label = label, 
       allele.lev = attributes(geno.glm)$unique.alleles, 
       control = haplo.glm.control(haplo.effect='dominant', 
       haplo.min.count=5))
print(fit.dom)


###################################################
### chunk number 37: seed
###################################################
set.seed(seed)


###################################################
### chunk number 38: 
###################################################
# control baseline selection, perform the same exact run as fit.bin, 
# but different baseline by using haplo.base chose from haplo.common
fit.bin$haplo.common
fit.bin$haplo.freq.init[fit.bin$haplo.common]
fit.bin.base77 <- haplo.glm(y.bin ~ male + geno.glm, family = binomial,
       data = my.data, na.action = "na.geno.keep",
       locus.label = label, 
       allele.lev = attributes(geno.glm)$unique.alleles, 
       control = haplo.glm.control(haplo.base=77, 
       haplo.min.count=5))
print(fit.bin.base77)


###################################################
### chunk number 39: 
###################################################
# merge haplo.score and haplo.group results
merge.bin <- haplo.score.merge(score.bin, group.bin)
print(merge.bin, nlines=10)


###################################################
### chunk number 40: seed
###################################################
set.seed(seed)


###################################################
### chunk number 41: 
###################################################

# demo haplo.cc where haplo.min.count is specified
# use geno, and this function prepares it for haplo.glm
cc.hla <- haplo.cc(y=y.bin, geno=geno, haplo.min.count=5,
       locus.label = label)
print(cc.hla, nlines=25, digits=2)
names(cc.hla)


###################################################
### chunk number 42: seed
###################################################
set.seed(seed)


###################################################
### chunk number 43: 
###################################################
# haplo.score on 11 loci, slide on 3 consecutive loci at a time
geno.11 <- hla.demo[,-c(1:4)]
label.11 <- c("DPB","DPA","DMA","DMB","TAP1","TAP2","DQB","DQA","DRB","B","A")
score.slide.gaus <- haplo.score.slide(resp, geno.11, trait.type =
                "gaussian", n.slide=3, skip.haplo=5/(2*nrow(geno.11)),
                locus.label=label.11)
print(score.slide.gaus)


###################################################
### chunk number 44: 
###################################################
# plot global p-values for sub-haplotypes from haplo.score.slide
plot(score.slide.gaus)


###################################################
### chunk number 45: seed
###################################################
set.seed(seed)


###################################################
### chunk number 46: 
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


