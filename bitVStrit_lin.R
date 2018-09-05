# Testing variations in cascade recording


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
barcodes = c(20,40,60,80,100)
#barcodes= c(2,10,50,100)
integrases = c(1,2,4)
#generations=c(8)
mus = c(0.5,0.4,0.3,0.15)
mus=c(0.4)
generations=c(4,6,8,10,12) #,10,12)
#mus = c(0.99,0.6,0.5,0.4,0.1,0.01)
#mus=c(0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001)
#barcodes = c(6,7)
#generations = c(3,4,5)
nRepeats=6

types=c('trit')
#types=c('trit')

#Main variation of all parameters including:
#Edit rate, barcodeLength , number of integrases , number of generations
#
cascadevar=list() #this is the main object
for(ca in 1:length(integrases)){
    muVariation=list()
    nIntegrases = integrases[ca]
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
              # this call is parallelized:
              xxx=compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats,nIntegrases = nIntegrases)
              genData[[ng]]= xxx[[1]] #the first element in this list is the results.matrix object that contains all the distance measures calcualted inside the funciton

              print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength)," N_ints=",toString(nIntegrases),sep=""))
           }
           barcodeData[[bc]]=genData
        }
        simType[[simulationType]]=barcodeData
      }
      muVariation[[m]]=simType
    }
    cascadevar[[ca]] = muVariation
    save(muVariation,file=paste("muVar_mu_",toString(mu),"_.rda",sep=""))
}
