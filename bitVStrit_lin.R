#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
library(gplots)
source("simBigMemory_lin.R")
source("simulation3_lin.R")
library(doParallel)

#request number of CPUs accordingly
os=system("cat ../os.txt",intern = T)
if(os=="mac"){
  registerDoParallel(cores=8)
}else if(os=="linux"){
  registerDoParallel(cores=6)
}

#SET PARAMETERS
barcodes = c(20,40,80)
#barcodes= c(2,10,50,100)
generations=c(4)
#generations=c(8)
mus = c(0.8,0.5,0.3)

generations=c(1,2,4)
#mus = c(0.99,0.6,0.5,0.4,0.1,0.01)
#mus=c(0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001)
#barcodes = c(6,7)
#generations = c(3,4,5)
nRepeats=100

types=c('trit')
#types=c('trit')

muVariation=list()
for(m in 1:length(mus)){
  simType=list()
  mu=mus[m]
  for(st in 1:length(types)){
    simulationType=types[st]
    barcodeData=list()
    for(bc in 1:length(barcodes)){
       barcodeLength=barcodes[bc]
       genData=list()
       for(ng in 1:length(generations)){
          nGen=generations[ng]
          #test :
          nIntegrases = nGen
            nGen =4 #always keep the same number of generations:
          xxx=compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats,nIntegrases = nIntegrases)
          genData[[ng]]= xxx[[1]] #the first element in this list is the results.matrix object that contains all the distance measures calcualted inside the funciton 
            
          print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength)," N_ints=",toString(nIntegrases),sep=""))
       }
       barcodeData[[bc]]=genData
    }
    simType[[simulationType]]=barcodeData
  }
  muVariation[[m]]=simType
  save(muVariation,file=paste("muVar_mu_",toString(mu),"_.rda",sep=""))
}
