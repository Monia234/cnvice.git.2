source("cnvice.r")

Nj = c(0,8,10,10,1) 
# Option 1
executeCnvice(Nj,F,document="output-report-HW")
# Option 2
executeCnvice(Nj,T,document="output-report-F")	
# Option 3
executeCnvice(Nj,T,rept=10,document="output-report-repeat10")

# Batch use
CNVice<-executeCnvice(Nj,T)
CNVice$allele_frequencies
CNVice$genotype_frequencies2
CNVice$genotype_frequencies

# Option 3: For genotype estimation using trio information
Nj= c(78, 829, 2510, 756, 83, 9, 1)
CNVice<-executeCnvice(Nj,T)
trio(4,3,3,CNVice$genotype_frequencies)