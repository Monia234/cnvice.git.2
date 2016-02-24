source("cnvice.r")
Nj= c(78, 829, 2510, 756, 83, 9, 1)
# Option 1
executeCnvice(Nj,T,rept=10,document="output-report-withF.doc")
# Option 2
executeCnvice(Nj,T,document="output-report-withF.doc")

# Option 3: For genotype estimation using trio information
egmatrix<-executeCnvice(Nj,T,rept=10,document="output-report-withF.doc")
trio(4,3,3,egmatrix)