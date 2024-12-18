# Selective genotyping and phenotyping strategies in a complex trait context S Sen, F Johannes, and KW Broman

Selective genotyping and phenotyping strategies can reduce the cost of
QTL (quantitative trait loci) experiments. We analyze selective
genotyping and phenotyping strategies in the context of multi-locus
models, and non-normal phenotypes. Our approach is based on calculations
of the expected information of the experiment under different
strategies. Our central conclusions are the following. (1) Selective
genotyping is effective for detecting linked and epistatic QTL as long
as no locus has a large effect. When one or more loci have large
effects, the effectiveness of selective genotyping is unpredictable \--
it may be heightened or diminished relative to the small effects case.
(2) Selective phenotyping efficiency decreases as the number of unlinked
loci used for selection increases, and approaches random selection in
the limit. However, when phenotyping is expensive, and a small fraction
can be phenotyped, the efficiency of selective phenotyping is high
compared to random sampling, even when over 10 loci are used for
selection. (3) For time-to-event phenotypes such as lifetimes, which
have a long right tail, right-tail selective genotyping is more
effective than two-tail selective genotyping. For heavy-tailed phenotype
distributions, such as the Cauchy distribution, the most extreme
phenotypic individuals are not the most informative. (4) When the
phenotype distribution is exponential, and a right-tail selective
genotyping strategy is used, the optimal selection fraction (proportion
genotyped) is less than 20% or 100% depending on genotyping cost. (5)
For time-to-event phenotypes where followup cost increases with the
lifetime of the individual, we derive the optimal followup time that
maximizes the information content of the experiment relative to its
cost. For example, when the cost of following up an individual for the
average lifetime in the population is approximately equal to the fixed
costs of genotyping and breeding, the optimal strategy is to follow up
approximately 70% of the population.

## **Disclaimer:** 

The source code and software distributed in this web page has no
implied warranty. Use at your own risk.

## **Symbolic computation code:** 

Some of the results in the paper were derived using symbolic
calculations in [Maxima](http://maxima.sourceforge.net). Start Maxima
in your system, and then you can cut and paste the contents of the
files below into the command window. The files are commented, so you
should be able to follow the steps.

-   Information calculations for linked QTL:
    [2qtl-linked.max](2qtl-linked.max)
-   Information calculations for epistatic QTL:
    [2qtl-epistatic.max](2qtl-epistatic.max)
-   Information gain functions, and expected information for non-normal
    distributions
    [info-calc.max](info-calc.max)

## **Figures:**

-   Figures 1 and 2:
    [2qtl-epistatic.R](2qtl-epistatic.R)
-   Figures 3 and 4:
    [info.pheno.R](info.pheno.R)
-   Figures 5 and 6:
    [non-normal.R](non-normal.R)
-   Figures 7 and 8:
    [info-calc.R](info-calc.R)
    [exponential-censor0.R](exponential-censor0.R)
-   Figure 9:
    [survival.R](survival.R)
-   Figure 10:
    [info-gain-censor.R](info-gain-censor.R)
-   Figure 11:
    [exponential-censor1.R](exponential-censor1.R)

