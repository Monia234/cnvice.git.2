######################################################################################################################
######## CNVice
######################################################################################################################

########  Install packages "hwriter" and "aylmer"
if( !require(aylmer) )	install.packages("aylmer")
if( !require(hwriter) )	install.packages("hwriter")

rm(list = ls()) #will remove ALL objects 

######################################################################################################################
######## This function calculates genotypic frequencies
########
######## j : copy number (integer)
######## nj: number of individuals with j copies (integer)
######## n : total number of individuals (integer)
######## jmaximo: maximum number of copies (integer)
######## pzero: initial values for allele frequencies (Real vector)
######## return ?: squared matrix of real numbers (jmaximo+1)X(jmaximo+1)
######## return : P(Hk,Hl)g according to paper (for only one j)
######################################################################################################################

func.ptil<-function(pzero,f){#calculates P~ based on F and pzero 
  jmaximo=length(pzero)-1
  Ptil <- matrix(0,length(pzero),length(pzero))
  for (k in 0:jmaximo){
    for (l in 0:jmaximo){
      if ((l>=k)&&(l+k<=jmaximo)){
        if ((k == l)){
          Ptil[k+1,l+1] <- f*pzero[k+1]+(1-f)*pzero[k+1]^2
         }
        if ((k != l)){
          Ptil[k+1,l+1] <- (1-f)*2*pzero[k+1]*pzero[l+1]
         }
        Ptil[l+1,k+1]=Ptil[k+1,l+1]
      }
    }
  }
  return(Ptil)
}

func.ptil_j=function(ptil0,j){#zeros Ptil cells with line+column != J
  Ptil=ptil0
  for(l in 0:(nrow(ptil0))-1)
  for(c in 0:(ncol(ptil0))-1){
    if(l+c!=j){
      Ptil[l+1,c+1]=0
    }
  }
  return(Ptil)
}

func.pj=function(ptil1){#calculates Pj based on Ptil
  jmaximo=nrow(ptil1)-1
  Pj <- 0
  pj = rep(0,jmaximo+1)

  for(j in 0:jmaximo){
    Ptilj=func.ptil_j(ptil1,j)
    Pj <- 0
    if(j>0){## Starts calculating Pj using Ptil previously calculated and updates Ptil to Pkl using the calculated Pj
      if (j%%2==1){
        Pj = sum(Ptilj[,1:((j/2)+1)])
      }
      if (j%%2==0){
        Pj <- sum(Ptilj[,1:((j+2)/2)])
      }
    }else{
      Pj <- sum(Ptilj[,])
    }
    pj[j+1]=Pj
  }
  return(pj)
}

func.pkl=function(nj,ptil3){#calculates Pkl based on Ptil

  jmaximo=length(nj)-1
  Pklj  <- matrix(0,jmaximo+1,jmaximo+1)
  Pkl=Pklj
  Pj <- 0
  n=sum(nj)

  for(j in 0:jmaximo){
    Ptilj=func.ptil_j(ptil3,j)
    Pj <- 0
    if(j>0){
      if (j%%2==1){
          Pj = sum(Ptilj[,1:((j/2)+1)])
          if(Pj==0){
            Pklj = 0
          }else{
            Pklj = (nj[j+1]/n)*(Ptilj[,]/Pj)
          }
      }
      if (j%%2==0){
           Pj <- sum(Ptilj[,1:((j+2)/2)])
           if(Pj==0){
              Pklj <- 0
           }else{
              Pklj <- (nj[j+1]/n)*(Ptilj[,]/Pj)
           }
      }
    }else{
      Pj <- sum(Ptilj[,])
      if(Pj==0){
        Pklj <- 0
      }else{
        Pklj <- (nj[j+1]/n)*(Ptilj[,]/Pj)
      }
    }
    Pkl=Pkl+Pklj
  }
  return(Pkl)
}

func.logPj=function(f,nj,p0){#calculates the product operator of Pj^Nj
  PTIL=func.ptil(p0,f)
  PJ=func.pj(PTIL)
  logPj = log(PJ)
  for (i in 1:length(PJ)){
    if (logPj[i] == -Inf){
    logPj[i] = NaN
    }
    }
  logPj = sum(nj*logPj,na.rm=T) 
  return(-logPj)
}

func.logPjTRV=function(f,nj,p0){#calculates the product operator of Pj^Nj for the TRV calculation
  PTIL=func.ptil(p0,f)
  PJ=func.pj(PTIL)
  logPj = log(PJ)
  for (i in 1:length(PJ)){
    if (logPj[i] == -Inf){
    logPj[i] = NaN
    }
    }
  logPj = sum(nj*logPj,na.rm=T) 
  return(logPj)
}

######################################################################################################################
######## This function calculates the final allele frequency
########
######## Nj : vector with number of individuals per copy number (starts with 0)
######## p0: initial values for allele frequencies
######## return : p as in paper (before iterations)
######################################################################################################################

func.p.alelo<-function(Nj,pzero1,f,valor){
  Pt2.j <- matrix(0,length(Nj),length(Nj))
  Palelo.j  <- rep(0,length(Nj))

  if(f==TRUE){
    F_OTIMO=optimise(func.logPj,nj=Nj,p0=pzero1,interval=c(0,1))$minimum
  }
  else{
    if(valor!=0){
      F_OTIMO=valor
    }
    else{
      F_OTIMO=0 #HWE
    }
  }
  PTIL=func.ptil(pzero1,F_OTIMO)
  PKL=func.pkl(Nj,PTIL)

    for (j in 1:length(Nj)){
      Palelo.j[j] = (1/2) * (sum(PKL[,j]))
      Palelo.j[j] =  Palelo.j[j] + (1/2)*(PKL[j,j])
   }
   Palelo=list(f_otimo=F_OTIMO,freq.alelica=Palelo.j)
   print(F_OTIMO)
  return(Palelo)
}


######################################################################################################################
######## Convergence function - until stable frequency
########
######## Nj : vector with number of individuals per copy number (starts with 0)
######## p0: initial values for allele frequencies
######## return: veter with p0 after convergence 
######################################################################################################################

CNVice <- function(Nj, pzero, maxit,f) {
  contf=0
  f_aux=f
  valor=0
  eps=1e-12
  out = matrix(0, nrow=maxit+1, ncol=length(Nj))
  dif = matrix(0, nrow=maxit+1, ncol=length(Nj))   
  final = length(Nj)
  out[1,1:final] =  func.p.alelo(Nj,pzero,f,valor)$freq.alelica
  i=1
  continue=T
  count = 0
  old=pzero
  while (continue) {
    novo      <- out[i,1:final]
    Pt.alelo  <- func.p.alelo(Nj,novo,f,valor)
    out[i+1,1:final]= Pt.alelo$freq.alelica
    dif[i,] <- round(abs(out[i+1,1:final]-out[i,1:final]),9)

    for (m in 1:length(Nj)){
      if (dif[i,m]<eps) {
         count = count + 1
      }
    }

    if(valor!=Pt.alelo$f_otimo){
      valor=Pt.alelo$f_otimo
      contf=0
    }
    else{
      contf=contf+1
    }
    if(contf>5){
      f=FALSE
    }

    i=i+1
    continue=(i <= maxit & count < length(Nj))
    old=novo
    count = 0
  }
  resp=list(cnvice=novo, f=func.p.alelo(Nj,novo,f,valor)$f_otimo)
  print(novo)
  return(resp)
}


######################################################################################################################
######## End of main functions
######################################################################################################################

monta.p0=function(Nj,repeticoes){## generates the matrix with random p0
  p0.rep=matrix(0,repeticoes,length(Nj))
  for(i in 1:repeticoes){
    p0.rep[i,]=rep(runif(length(Nj), min=0, max=1))
    p0.rep[i,]=p0.rep[i,]/sum(p0.rep[i,])
    #p0.rep[i,]=Nj/sum(Nj)
  }
  return(p0.rep)
}

######################################################################################################################
######## Repetitions
######################################################################################################################

# runs the CNVice function N times and returns a matrix with one result per line (it saves each iteration)
# filename: file RData for backup, or backup files with data already generated
# apagarRData: TRUE->deletes backup file after the process finishes
# if the file already exists and the number of repetitions are smaller than the original file, the algorithm does not run
CNVice.rep=function(Nj,repeticoes,f,filename='CNVice.RData',apagarRData=TRUE){
  options(warn=-1)
  print('begin: CNVice.rep')
  print(date())
  p0.rep=monta.p0(Nj,repeticoes)
  if(file.exists(filename)){#carrega o arquivo
	load(filename)
	print('CNVice file loaded')
	file.copy(filename,'bkpCNViceRData')
  }else{#inicializa as variáveis
	print('CNVice file generated')
	print('CNVice repetition: 1')
	aux=CNVice(Nj, p0.rep[1,], 1000,f)
	Fs = aux$f
	p.rep = round(aux$cnvice,7)
	save(p.rep,Fs,file=filename)
	if(file.exists('tmpCNVice.RData'))
		file.remove('tmpCNVice.RData')
  }
  if(repeticoes>length(Fs))
  for (m in (length(Fs)+1):repeticoes){
    print(paste('CNVice repetition: ',m))
    aux=CNVice(Nj, p0.rep[m,], 1000,f)
	p.rep=rbind(p.rep,round(aux$cnvice,7))
	Fs=append(Fs,aux$f)
	file.copy(filename,'tmpCNVice.RData')
	save(p.rep,Fs,file=filename)
	if(file.exists('tmpCNVice.RData'))
		file.remove('tmpCNVice.RData')
  }
  if(file.exists('bkpCNVice.RData'))
	file.remove('bkpCNVice.RData')
  print('end: CNVice.rep')
  print(date())
  lista=list(cnvice=p.rep,f=Fs)
  if(apagarRData)
	file.remove(filename)
  return(lista)
}

CNVice.rep.TRV=function(Nj,repeticoes,f){
  options(warn=-1)
  p0.rep=monta.p0(Nj,repeticoes)
  aux=CNVice(Nj, p0.rep[1,], 1000,f)
  Fs = aux$f
  p.rep = round(aux$cnvice,7)
  
  if(repeticoes>length(Fs))
    for (m in (length(Fs)+1):repeticoes){
      print(paste('CNVice repetition: ',m))
      aux=CNVice(Nj, p0.rep[m,], 1000,f)
      p.rep=rbind(p.rep,round(aux$cnvice,7))
      Fs=append(Fs,aux$f)
    }  
  lista=list(cnvice=p.rep,f=Fs)
  
  return(lista)
}

######################################################################################################################
######## Calculates allele frequencies based on allele probabilities
######################################################################################################################
calc.esp=function(p,Nj){ #calculates "Expected bin counts" (random count)
  total.esp = rep(0,length(p))
  CNV.esp = rep(0,length(p))
  for (j in 0:length(p)-1){
    if(j==0){
      total.esp[j+1] = total.esp[j+1] + p[j+1]^2
    }else{
      if (j%%2==0){
        for (k in 0:(j/2)){
          l = j-k
          if (l!=k){
            total.esp[j+1] = total.esp[j+1] + 2.0*p[k+1]*p[l+1]
          }else{
            total.esp[j+1] = total.esp[j+1] + p[k+1]^2
          }
        }
      }
      if (j%%2==1){
        for (k in 0:(((j+1)/2)-1)){
          l= j-k
          total.esp[j+1] = total.esp[j+1] + 2.0*p[k+1]*p[l+1]
        }
      }
    }
    CNV.esp[j+1] = round((total.esp[j+1] * sum(Nj)), digits = 0)
  }
  return(CNV.esp)
}


#####################################################################
### Stores individual genotype frequency
Pgenot.ind <- function(j,matriz){
  linhas=NULL
  valores=NULL
  if(j==0){
    linhas=append(linhas,c("00"))
    valores=append(valores,c(matriz[1,1]))
  }
  else
  {
    for(l in 1:nrow(matriz))
      for(c in 1:ncol(matriz)){
        if(l+c-2==j){
          linhas=append(linhas,paste(l-1,c-1,sep=""))
          valores=append(valores,matriz[l,c])
        }
      }
  }
  if(sum(valores)!=0)
    valores=valores/sum(valores)
  ind=data.frame(valores,row.names=linhas)
  names(ind)="Prob"
  return(ind)
}# to look at result, resultado["02",] gets line 02

################################################################################
## Improve estimates of child given parents

pai.mae<-function(x,y){ #child with father x copies and mother y copies return matrix where col1=alelle1 e col2=alelle2)
  cont=0
  paimae=NULL
  if(x==0&&y==0){
    paimae=rbind(paimae,c(x,y))
    cont=cont+1
  }
  else{
    for(i in 0:x)
    for(j in 0:y){
      cont=cont+1
      paimae=rbind(paimae,c(i,j))
      if(i==x/2){
        cont=cont+1
        paimae=rbind(paimae,c(i,j))
      }
      if(j==y/2){
        cont=cont+1
        paimae=rbind(paimae,c(i,j))
      }
      if(i==x/2&&j==y/2&&i==j){
        cont=cont+1
        paimae=rbind(paimae,c(i,j))
      }
    }
  }
  return(paimae)
}

trio <- function(pai,mae,filho,matriz){ # Child probability given parents and genotype matrix - return(data.frame)
  filhos=NULL
  paimae=pai.mae(pai,mae)
  for(i in 1:nrow(paimae)){
    if(sum(paimae[i,])==filho){
      filhos=append(filhos,paste(paimae[i,],sep="",collapse=""))
    }
  }
  filhos=as.data.frame(table(filhos))

  for(i in 1:nrow(filhos)){#identifies the pair of mother and father
    pri=as.integer(substr(filhos$filhos[i],1,1))
    seg=abs(pri-pai)
    if(pri>seg){
      pri=seg
      seg=abs(pri-pai)
    }
    filhos$pai[i]=paste(pri,seg,sep="")##father's pair

    pri=as.integer(substr(filhos$filhos[i],2,2))
    seg=abs(pri-mae)
    if(pri>seg){
      pri=seg
      seg=abs(pri-mae)
    }
    filhos$mae[i]=paste(pri,seg,sep="")##mother's pair
  }

  filhos$paimae=paste(filhos$pai,filhos$mae,sep="")##probability calculation
  paimae=as.data.frame(table(filhos$paimae))
  names(paimae)=c("paimae","Freq2")
  filhos=merge(filhos,paimae)
  filhos$paimae=NULL
  filhos$prob.filho.dado.paimae=filhos$Freq/(filhos$Freq*filhos$Freq2)
  filhos$Freq2=NULL

  res.pai=Pgenot.ind(pai,matriz)##probability of father, given number of copies
  res.pai$pai=rownames(res.pai)
  names(res.pai)=c("prob.pai","pai")
  filhos=merge(filhos,res.pai)

  res.mae=Pgenot.ind(mae,matriz)##probability of mother, given number of copies
  res.mae$mae=rownames(res.mae)
  names(res.mae)=c("prob.mae","mae")
  filhos=merge(filhos,res.mae)

  for(i in 1:nrow(filhos)){##corrects probabilities for mother and father
    if(substr(filhos$mae[i],1,1)!=substr(filhos$mae[i],2,2)){
      filhos[i,]$prob.mae=filhos[i,]$prob.mae*2
    }

    if(substr(filhos$pai[i],1,1)!=substr(filhos$pai[i],2,2)){
      filhos[i,]$prob.pai=filhos[i,]$prob.pai*2
    }
  }

  cat("Offspring?s genotype probabilities:\n")
  filhos$prob=filhos$prob.filho.dado.paimae*filhos$prob.pai*filhos$prob.mae
  filhos$prob.final=round(filhos$prob/sum(filhos$prob),7)##final probability
  filhos$prob.mae=round(filhos$prob.mae,7)
  filhos$prob.pai=round(filhos$prob.pai,7)
  return(filhos[,c(2,1,3,6,7,9)])
}

# rows of A that are in B 
rowmatch <- function(A,B) {
    f <- function(...) paste(..., sep=":") 
   if(!is.matrix(B)) B <- matrix(B, 1, length(B)) 
    a <- do.call("f", as.data.frame(A)) 
    b <- do.call("f", as.data.frame(B)) 
    match(b, a) 
} 

#resumes the result of CNVice.rep to the one more probable
resumeCnvice=function(resultadoCnvice,decimaisProb=3,decimaisF=3){
	cnvice=round(resultadoCnvice$cnvice,decimaisProb)
	linhas=apply(cnvice,1,paste,collapse=',') #convert rows to strings for counting
	moda=subset(table(linhas),table(linhas)==max(table(linhas)))
	pzeroE=eval(parse(text=paste("c(",names(moda),")")))
	#f=round(resultadoCnvice$f,decimaisF)
	#moda=subset(table(f),table(f)==max(table(f)))
	#F=eval(parse(text=paste("c(",names(moda),")")))
	F=round(resultadoCnvice$f[rowmatch(cnvice,pzeroE)],decimaisF)
	return(list("cnvice"=pzeroE,"f"=F))
}

resumeCnvice_TRV=function(resultadoCnvice,decimaisProb=3,decimaisF=3){
  cnvice_res=round(resultadoCnvice$cnvice,decimaisProb)
  linhas=apply(cnvice_res,1,paste,collapse=',') #convert rows to strings for counting
  
  Fs_est <- round(resultadoCnvice$f,decimaisProb)
  F_max=max(Fs_est)   #subset(table(f),table(f)==max(table(f)))
  moda_F=subset(linhas,Fs_est==F_max)
  
  moda=subset(table(linhas),table(linhas)==max(table(moda_F)))
  pzeroE=eval(parse(text=paste("c(",names(moda),")")))
  
  return(list("cnvice_res"=pzeroE,"f"=F_max))
}

##################################################################
# TRV
# d1 = -2*func.logPjTRV(0,NjObservado,pzeroE.0) -- with F=0
#      +2*func.logPjTRV(F,NjObservado,pzeroE.F) -- with F<>0
##################################################################
## Calculates TRV between F=0 and F<>0

TRV=function(Nj,repeticoes){
	p.repF0=CNVice.rep.TRV(Nj,repeticoes,FALSE) #calls "Estimated allele frequencies" function (allele probabilities)
	resultadoF0=resumeCnvice(p.repF0,7,7)

	p.repFs=CNVice.rep.TRV(Nj,repeticoes,TRUE) #calls "Estimated allele frequencies" function (allele probabilities)
	resultadoFs=resumeCnvice_TRV(p.repFs,7,7)
# 	
# 	if(apagarRData){
# 		file.remove(filename_FS)
# 		file.remove(filename_F0)
# 	}
	
	print('')
	print('Generated values')
	print('')
	
	print('Prof with F=0')
	print(resultadoF0$cnvice)
	print('')
	
	print('Prob with F<>0')
	print(resultadoFs$cnvice_res)	
	print('')
	
	print('F values')	
	print(resultadoFs$f)
	print('')
	
	print('D1:')
	d1 = -2*(func.logPjTRV(0,Nj,resultadoF0$cnvice)-func.logPjTRV(resultadoFs$f,Nj,resultadoFs$cnvice_res))
	print(d1)
	
	#print('pValor:')
	#print(0.5+pchisq(d1,1)/2)
	#maira: changed 2.42075 to 1.9205 -- 5% confidence interval 
	return(list("trv"=d1,"rejeita"=(d1>1.9205),"resultadoF0"=resultadoF0,"resultadoFs"=resultadoFs))
}

######################################################################################################################
######## Execution
######################################################################################################################
executeCnvice <- function(Nj,estimar_F,rept=100,document="CNViceReport",decimaisProb=3,decimaisF=3,result=1,filename_backup='backup.RData',apagarRData=TRUE){
#estimar_F: estimate F statistic or consider it zero
#document: document name where the report will be written
#result=1 : show best result
#filename_backup: back up file name
#apagarRData: delete backup file? (TRUE or FALSE)

  #rm(list = ls()) #will remove ALL objects 
  
  starttime = Sys.time()
  if(file.exists('bkpCNVice.RData'))
    file.remove('bkpCNVice.RData')
  if(file.exists(filename_backup))
    file.remove(filename_backup)
  
  resultado=result
  ## calculates allele probabilities and stores
  #estimar_F=TRUE
  p.rep=CNVice.rep(Nj,rept,estimar_F,filename_backup,FALSE) #calls function that calculates the estimated allele frequencies (alleleic probabilities)
  p.rep$cnvice=round(p.rep$cnvice,decimaisProb)
  p.rep$f=round(p.rep$f,decimaisF)
  linhas=apply(p.rep$cnvice,1,paste,collapse=',') #change to string
  #table(linhas) #counts
  #moda=subset(table(linhas),table(linhas)==max(table(linhas)))
  #moda=eval(parse(text=paste("c(",names(moda),")")))
  respostas=names(table(linhas)) #generates matrix with allelic probabilities 
  contagem=as.vector(table(linhas))
  ps=matrix(0,nrow=length(respostas),ncol=ncol(p.rep$cnvice)) 
  for(i in 1:length(respostas)){ #generates matrix with allelic probabilities (as number)
    ps[i,]=eval(parse(text=paste("c(",respostas[i],")")))
  }
  fs=p.rep$f[rowmatch(p.rep$cnvice,ps)]
  ######################################################################################################################

  ######################################################################################################################
  ## Prints statistics and graphs with allelic probabilities

 #table(linhas) 
 # table(round(fs,7)) #F counts
  esperado=ps
  matrizes = array(0, dim=c(length(Nj),length(Nj),length(ps[,1])))
  pvalores=rep(0,length(respostas))
  #par(mfrow=c(1,length(esperado[,1])+1))
  #par(mfrow=c(length(esperado[,1])+1,1))
 #barplot(rev(Nj),horiz=TRUE,main="Observado",names.arg=(length(Nj)-1):0)
  #hist(rep(0:(length(Nj)-1),Nj),labels=T,breaks=c(-1:(length(Nj)-1)),main="Observado",ylab=NULL,xlab=NULL,ylim=c(0,sum(Nj)))
  for(i in 1:length(ps[,1])){
    esp=calc.esp(ps[i,],Nj)
    pvalores[i]=chisq.test(Nj,esp)$p.value
    #pvalores[i]=ks.test(Nj,esp)$p.value
    #barplot(rev(esp),horiz=TRUE,names.arg=(length(esp)-1):0,main=paste("Esperado",i))
    #hist(rep(0:(length(esp)-1),esp),labels=T,breaks=c(-1:(length(esp)-1)),main=paste("Esperado",i),ylab=NULL,xlab=NULL,ylim=c(0,sum(Nj)))
    #text(0,max(esp)+20,paste("p-valor:",round(pvalores[i],3)))
   # print(paste("############### Resultado ",i," ###############"))
    #print(ps[i,])
    cat('\nF: ',fs[i],'\n')
    cat('observed: ',Nj)
    cat('\nestimated: ',esp)
    #try(print(chisq.test(Nj,p=esp,rescale.p=T)))   #Teste qui-quadrado
    #try(print(fisher.test(matrix(c(Nj,esp),nrow=length(ps[,1]),byrow=T))))  #Teste Exato de Fisher
    #try(print(ks.test(Nj,esp)))
    esperado[i,]=esp
    ptil_temp=func.ptil(ps[i,],fs[i])##fazer calcular o ptil com o F final estimado
    matrizes[,,i]=func.pkl(Nj,ptil_temp)##matrizes[,,i]=func.pkl(Nj,ps[i,])##fazer calcular o pkl com o ptil estimado a cima
  }
  matrizes=round(matrizes,decimaisProb)
  
  if(estimar_F){ #TRV=function(Nj,repeticoes,filename_F0='TRV_F0.RData',filename_FS='TRV_Fs.RData',apagarRData=TRUE)
	  trv=TRV(Nj,rept)
  }else{
	print("TRV não precisa ser feito") #trv=TRV(Nj,rept,filename_F0=filename_backup,apagarRData=apagarRData)
  }

  #dev.new()
  #par(mfrow=c(1,length(esperado[,1])+1))
  #barplot(rev(Nj),horiz=TRUE,main="Observado",names.arg=(length(Nj)-1):0)
  #for(i in 1:length(ps[,1])){
  # esp=calc.esp(ps[i,],Nj)
  #  barplot(rev(esp),horiz=TRUE,names.arg=(length(esp)-1):0,main=paste("Esperado",i))
  # text(max(esp)/2,1,paste("p-valor:",round(pvalores[i],3)))#text(0,max(esp),paste("p-valor:",round(pvalores[i],3)))
  #  esperado[i,]=esp
  #}

  Nj #observed number of individuals with j copies (integer)
  ps #estimated (relative) allelic frequencies
  fs #F's corresponding to ps
  esperado #estimated expected values (based on ps)
  contagem #ps count (how many appearences of the result)
  matrizes #Estimated genotype frequencies
  


  ######################################################################################
  ## choosing the best result
  #resultado=1 
  expected.bin=esperado[resultado,]
  estimated.allele=ps[resultado,]
  estimated.genotype.mult2=matrizes[,,resultado]


  #######################################################################################
  ## correcting matrix of "Estimated genotype frequencies"
  estimated.genotype=estimated.genotype.mult2/2
  for( i in 1:nrow(estimated.genotype))
    estimated.genotype[i,i]=estimated.genotype[i,i]*2


  ##data formating
  matriz=estimated.genotype.mult2
  for(l in 1:nrow(estimated.genotype)){
    for(c in 1:ncol(estimated.genotype)){
      if(l>c){
        matriz[l,c]=NA
      }
    }
  }
  #maira
  nomes=NULL
  for(l in 1:nrow(estimated.genotype)){
	nomes[l]=paste(l-1)
  }
  matriz = as.data.frame(matriz)
  names(matriz)=nomes
  row.names(matriz)=nomes
  #end-maira
  
  tabela=data.frame(0:(length(Nj)-1),Nj,t(esperado),t(ps))
  nomes=NULL
  nomes[1]='j'
  nomes[2]='Nj'
  for(i in 1:nrow(esperado)){
    nomes[i+2]=paste('expected.bin',i,sep='_')
  }
  j=1
  i=i+2
  while(i<ncol(tabela)){
    i=i+1
    nomes[i]=paste('freq.alel.est',j,sep='_')
    j=j+1
  }
  names(tabela)=nomes

  if(nrow(esperado)==2)
  tabela=tabela[,c(1,2,3,5,4,6)]
  if(nrow(esperado)==3)
  tabela=tabela[,c(1,2,3,6,4,7,5,8)]

  endtime = Sys.time()

  if(document != ""){
    library(hwriter)
    ##printing data
    a=openPage(paste(document,".doc",sep=""))
    hwrite('Data Summary',a,a)
    hwrite(tabela,a,a,br=TRUE,borders=0) #data summary
  
    hwrite('Alternative estimations count: ',a,a)
    hwrite(matrix(c(1:length(contagem),contagem),ncol=length(contagem),byrow=T),a,a,br=TRUE) #ocurrence count
  
    hwrite('F values found for each estimation: ',a,a)
    hwrite(matrix(c(1:length(fs),fs),ncol=length(fs),byrow=T),a,a,br=TRUE) #F values
  
    hwrite('Likelihood ratio (LR) test: ',a,a,br=TRUE)
    hwrite('Value of statistic found by LR: ',a,a)
    hwrite(round(trv$trv,5),a,a,br=TRUE)
    hwrite('Reject H0, i.e., F different than ZERO?',a,a,br=TRUE)
    if(trv$rejeita){
  	hwrite('Yes',a,a,br=TRUE)
    }else{
  	hwrite('No',a,a,br=TRUE)
    }
  
    hwrite('<br>P-values for chi-squared tests:',a,a)
    hwrite(matrix(c(1:length(pvalores),round(pvalores,decimaisProb)),ncol=length(pvalores),byrow=T),a,a,br=TRUE)
  
    hwrite('<br>Population genotype probability:',a,a)
    hwrite(matriz,a,a,br=TRUE) #estimated.genotype
    hwrite("Individual genotype probability (given the individual's copy number)",a,a)
    for(i in 0:(length(Nj)-1)){
      hwrite((t(round(Pgenot.ind(i,estimated.genotype),decimaisProb))),a,a)
    }
    hwrite('<br>',a,a)
    #hwrite('Probabilidade Genotipica dado que o filho tem x copias, dado que o Pai tem y copias e a Mae z copias.',a,a)
    #hwrite(trio(4,3,3,estimated.genotype),a,a,br=TRUE)
    hwrite('Runtime (start,end):<br>',a,a)
    hwrite(toString(starttime),a,a)
    hwrite('<br>',a,a)
    hwrite(toString(endtime),a,a,br=TRUE)
    print("Total time: ")
    print(starttime)
    print(endtime)
    closePage(a,splash=F)
  }
  
  file.remove(filename_backup)
  print("Backup files removed.")
  
  # OUTPUT:
  # Estimated allele frequencies (tabela)
  # Estimated population genotpye frequencies (estimated.genotype)
  # Estimated individual genotype frequencies (?)
  
  returnList <- list("allele_frequencies" = tabela, "genotype_frequencies" = estimated.genotype, "genotype_frequencies2" = estimated.genotype.mult2)
  return(returnList)
}