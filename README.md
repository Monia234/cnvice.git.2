# cnvice

CNVice (Inbreeding Coefficients Estimation for CNV data) is a freely available R script for population genetics applications.
CNVice performs the following analyses:
1.	Estimates allele frequencies for a given population assuming Hardy-Weinberg equilibrium (HWE), using the expectation maximization algorithm as implemented in CoNVEM (Gaunt et al. 2010), conditioning on the observed diplotype frequencies. 
2.	Estimates the population structure parameter  CNV and allele frequencies, using the profiled likelihood function and the Expectation Maximization (EM) algorithm (Dempster et al. 1977)
3.	Estimates the population genotype frequency, conditioning on the observed diplotype distribution and the estimated fCNV and allele frequencies.
4.	Uses trio information to improve the inference of an offspring’s genotype, by considering the parents’ diplotypes and the population genotype frequency.

Citing us: the CNVice paper has been submitted for publication and will be soon available here. 

#Requirements
CNVice is implemented in R. Packages aylmer and hwriter are installed when CNVice is executed.
 
#Input
The main input for CNVice is the distribution of observed diplotype frequencies of a given loci in a population. The distribution must always start with diplotype 0 (see the Examples section).  

#Running CNVice
The main function to be executed is: 
executeCnvice(Nj,fpar,rept,document)

Please read the user guide for more details.
