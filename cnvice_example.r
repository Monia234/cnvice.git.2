source("cnvice.r")
Nj = c(0,8,10,10,1)
# Option 1
executeCnvice(Nj,T,rept=10,document="output-report-withF")
# Option 2
executeCnvice(Nj,T,document="output-report-withF")
# Option 3
executeCnvice(Nj,F,document="output-report-withoutF")

# Option 4: For genotype estimation using trio information
Nj= c(78, 829, 2510, 756, 83, 9, 1)
egmatrix<-executeCnvice(Nj,T,rept=10,document="output-report-trio-withF")
trio(4,3,3,egmatrix)
