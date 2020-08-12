# haplo.stats
Routines for the analysis of indirectly measured haplotypes. The statistical 
methods assume that all subjects are unrelated and that haplotypes are ambiguous 
(due to unknown linkage phase of the genetic markers).

The main functions are described below.


## haplo.em

Estimation of haplotype frequencies, and posterior probabilities of haplotype 
pairs for a subject, conditional on the observed marker data.

## haplo.glm

GLM regression models for the regression of a trait on haplotypes, 
possibly including covariates and interactions. S3 methods for anova and summary 
have been implemented.

## haplo.score

Score statistics to test associations between haplotypes and a wide 
variety of traits, including binary, ordinal, quantitative, and Poisson.

## seqhap

Sequential haplotype scan methods to perform association analyses
for case-control data.  When evaluating each locus, loci that 
contribute additional information to haplotype associations with 
disease status will be added sequentially.

## haplo.design

Uses as input the result from haplo.em(), and makes a design matrix for haplotype 
dosage, such that modeling haplotypes is similar to how it would be done within haplo.glm(), but without the iteratetively re-weighted least squares steps.

## haplo.cc

Runs simple haplo.score and haplo.glm without covariates with combined results for case-control (binomial family) response.
