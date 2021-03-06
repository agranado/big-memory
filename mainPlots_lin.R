
#Oct 9th,2018
#First round of cascade simulations is Done!
#Here we will generate basic plots for comparing 1 vs 2 integrases as independent channels using a new reconstruction method
#the goal:show as a proof of principle that by applying smarter memory usage we can optimize lineage recording


#this script should load the object muVariation or it should be executed afeter running bitVStrit.R
os=system("cat ../os.txt",intern = T) #Local Mac repository (laptop)
if( os=="linux"){
  data.path = "/home/alejandrog/MEGA/Caltech/lineage/simulation_data/" #Location of simulation data from AWS
}else if(os=="mac"){

  data.path = "./"
}
#parameters from AWS: large object cascadevar 10-Sep-2018
#SET PARAMETERS
barcodes = c(20,30,40,50,60,80,100,200)
integrases = c(1,2)
mus = c(0.5,0.4,0.3,0.2,0.1,0.05)
generations=c(4,5,6,7,10)#,11,12)
nRepeats=72
types=c('trit')
#large file with all the simulations: distances and trees.
# 1 muVariation object per integrase number, for cascade recording
file.name.cascade1="intVar_int_1_.rda"
file.name.cascade2="intVar_int_2_.rda"
load(paste(data.path,file.name.cascade1,sep=""))
cascadevar1=cascadevar; rm(cascadevar)
load(paste(data.path,file.name.cascade2,sep=""))
cascadevar2=cascadevar ;rm(cascadevar)

cascadevar=list()
cascadevar[[1]] = cascadevar1[[1]]
cascadevar[[2]] = cascadevar2[[1]]
rm("cascadevar1","cascadevar2")
#name of the object is cascadevar
# struct[[int_number]][[mu]][[typeSim]][[bc]][[ng]]
# DIMS : cascadevar[[1]][[5]][[1]][[6]][[6]]

# RANDOM distance: for normalizaiton

ran.generations = c(4,5,6,7,10)#,11,12)
a=list()
  #RUN
  #results= foreach(i=1:length(ran.generations)) %dopar% simMemoirRandomGuess(ran.generations[i],mu,alpha,barcodes[1],nRepeats)
  #pre-runned
  a.long=c(0,0,9.2,   24.0,   56.0,  120.0,  248.4,  504.4, 1016.8, 2042.0)
  a = a.long[generations]
#plotting ----------------------------------------------------------------
manual.idx=3
upgam.idx=1
alpha=2/3

distBit=array(0,dim=c(length(barcodes),length(generations)))
distTrit=array(0,dim=c(length(barcodes),length(generations)))

# # # # # # #
pathPlots="/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay1st/"
#plot for each generation, RF dist vs BC, overlap blue and black line for trit vs bit
#estimated edit rate 0.4 edits per site per generation
a<-c(  26,   120, 506)

int_index = 1

x11()
#pdf(paste(pathPlots,"RFdist_vsBC_empiricMu_compareGenerations.pdf",sep=""))
par(mfrow=c(length(integrases),length(generations)))


cascade.matrices=list()
#method number 3 only applies for more than 1 integrase
method.idx=c(2,3)
for (ni in 1:length(integrases)){
  distTrit=array(0,dim=c(length(barcodes),length(generations)))

  muVariation = cascadevar[[ni]]
  mIdx=3 #only one edit rate
  simType=muVariation[[mIdx]]
    for(ng in 1:length(generations)){

      for(bc in 1:length(barcodes)){
        #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
        distTrit[bc,ng]=apply(simType[['trit']][[bc]][[ng]],2,mean)[method.idx[ni]] # 12 -> RF.dist using manualDist + UPGMA

      }
    }
    cascade.matrices[[ni]]=distTrit
 }


dev.off()

#normalize distances with the random guess and plot heatmaps
distBit = cascade.matrices[[1]]
distTrit = cascade.matrices[[2]]
distBitNorm=array(0,dim=dim(distBit));for(d in 1:dim(distBit)[2]){distBitNorm[,d]=(a[d]-distBit[,d])/a[d]}

distTritNorm=array(0,dim=dim(distTrit));for(d in 1:dim(distTrit)[2]){distTritNorm[,d]=(a[d]-distTrit[,d])/a[d]}

#rename columns

rownames(distBitNorm)=as.character(barcodes)
colnames(distBitNorm)=as.character(generations)

rownames(distTritNorm)=as.character(barcodes)
colnames(distTritNorm)=as.character(generations)

#plot heatmaps
#normalized distances (by random tree) using the empirical edit rate of 0.4
plotPath= "/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay1st/"
pdf(paste(plotPath,"heatmap_Distance_BitVSTrit_muNOTRACE",toString(mus[mIdx]),".pdf",sep=""))
heatmap.2(t(distBitNorm),dendrogram='none', Rowv=FALSE, Colv=FALSE,col=paste("gray",1:99,sep=""),main="Bit",trace="none",xlab="Barcodes",ylab="Generations")
heatmap.2(t(distTritNorm),dendrogram='none', Rowv=FALSE, Colv=FALSE,col=paste("gray",1:99,sep=""),main="Trit",trace="none",xlab="Barcodes",ylab="Generations")

vals=unique(scales::rescale(c(distBitNorm)))
o<-order(vals,decreasing=F)
cols<-scales::col_numeric("Blues",domain=NULL)(vals)
colz<-setNames(data.frame(vals[o],cols[o]),NULL)



p1=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=distBitNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))
p2=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=distTritNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p2,limits=c(0,1))
subplot(p1,p2)


dev.off()

#compare the different distance measures:
distAllCompare = array(0,dim=c(12,length(generations),length(mus)))
par(mfrow=c(2,3))
for(ng in 1:length(generations)){
  for (mIdx in 1:length(mus)){
    simType=muVariation[[mIdx]]

    for(d in 1:12){
      #distBitAll[bc,ng,mIdx]=apply(simType[['binary']][[bc]][[ng]],2,mean)[12] #12 is normal UPGMA + as.dist
      distAllCompare[d,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[d]
      #distBitAll[bc,ng, mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[12]

    }
  }
  #plot(mus,1-distTritAll[3,ng,]/a[ng],ylim=c(0,1),
  #     type="o",col="blue",ylab="norm dist",xlab="edit rate p/site p/gen",
  #     main=paste("g=",toString(generations[ng]),sep=""),
  #     cex.axis=1.5,cex.lab=1.6,cex.main=2);
  #lines(mus,1-distBitAll[3,ng,]/a[ng],type="o")

  #lines(mus,1-distTritAll[5,ng,]/a[ng],col="#3399FF",pch=17,lty=6)
  #lines(mus,1-distBitAll[5,ng,]/a[ng],col="gray",pch=17,lty=6)

}

# # # # # # 
 # # # # 



# Sep 10th
# convert the big list into a more friendly matrix
compare.cascade=list()
for (ca in 1:length(integrases)){
    distBitAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
    distTritAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
    par(mfrow=c(2,3))
    upgma.idx=1
    manualDist.idx=3
    BCn = 3
    muVariation = cascadevar[[ca]]
    for(ng in 1:length(generations)){
      for (mIdx in 1:length(mus)){
        simType=muVariation[[mIdx]]
        for(bc in 1:length(barcodes)){
          distBitAll[bc,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[upgma.idx] #12 is normal UPGMA + as.dist
          distTritAll[bc,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[manualDist.idx]
        }
      }



      plot(mus,1-distTritAll[BCn,ng,]/a[ng],ylim=c(0,1),
           type="o",col="blue",ylab="norm dist",xlab="edit rate p/site p/gen",
           main=paste("g=",toString(generations[ng]),sep=""),
           cex.axis=1.5,cex.lab=1.6,cex.main=2);
      lines(mus,1-distBitAll[BCn,ng,]/a[ng],type="o")
    }
    compare.cascade[[ca]]=distBitAll
}


#take minimun across all edit rates
#optimal rate
# rows are barcodes
# cols are generations
compare.cascade.optimal = list()
for(ca in 1:length(integrases)){

    distTritAll = compare.cascade[[ca]]
    optimMatBit=apply(distTritAll,c(1,2),min)
    compare.cascade.optimal[[ca]] = optimMatBit
}


#for a given number of generations: how does reconstruction improves with increasing number of barcodes

ng = 2
# R is for reconstructability
x11()
par(mfrow=c(2,3))
for( ng in 1:length(generations)){
  R_vs_barcode = matrix(0,length(integrases),dim(compare.cascade.optimal[[1]])[1] )
  for (ca in 1:length(integrases)){

    optimMatBit = compare.cascade.optimal[[ca]]
    R_vs_barcode[ca,] = optimMatBit[,ng]
  }
    R_vs_barcode=1-R_vs_barcode/a[ng]
    matplot(barcodes,t(R_vs_barcode),type="l",main =paste("g=",toString(generations[ng])),log="x")

}
lends <- c("1","2","4")
text(cbind(30, 50*c(1,3,5)-.4), lends, col= 1:3, cex = 1.5)



pdf(paste(pathPlots,"histograms_tenmer_3gen_muVar.pdf"))
par(mfrow=c(2,2))
ng=1;bc=2
for(m in 1:4){
  medDist = mean(1-muVariation[[m]][['trit']][[bc]][[1]][,13]/a[ng])
  medianDist = median(1-muVariation[[m]][['trit']][[bc]][[1]][,13]/a[ng])
  hist(1-muVariation[[m]][['trit']][[bc]][[ng]][,13]/a[ng],
       main=paste("u=",toString(mus[m]),", <d>=",toString(medDist),", me=",toString(medianDist),sep=""),
       cex.lab=1.5,cex.axis=1.5,cex.main=1.7,
       xlab="norm distance",xlim=c(0,1),ylim=c(0,45))
}
dev.off()



# across generations


pdf(paste(pathPlots,"RFdist_vsBC_optimalRate.pdf",sep=""))
par(mfrow =c(2,3))
barcodeAxis= log(barcodes)
for(ng in 1:length(generations)){
  plot(barcodeAxis,(a[ng]-optimMatTrit[,ng])/a[ng],
       cex.main=2,cex.lab=1.5,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N scratchpads",ylab="Dist norm to random guess",lwd=1.5);
  lines(barcodeAxis,(a[ng]-optimMatBit[,ng])/a[ng],col="black",type="o",lwd=1.5)
  abline(h=a[ng],col="red")
}
dev.off()

pdf(paste(pathPlots,"RFdist_vsBC_optimalRate_vsEmpiricalBIT.pdf",sep=""))
par(mfrow =c(2,3))
barcodeAxis= log(barcodes)
for(ng in 1:length(generations)){
  plot(barcodeAxis,(a[ng]-optimMatBit[,ng])/a[ng],
       cex.main=2,cex.lab=1.5,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N scratchpads",ylab="Dist norm to random guess",lwd=1.5);
  lines(barcodeAxis,distBitNorm[,ng],col="orange",type="o",lwd=1.5)
  abline(h=a[ng],col="red")
}
dev.off()




#heatmaps of RF.dist for BC vs nG
optimMatBitNorm=array(0,dim=dim(optimMatBit));
optimMatTritNorm=array(0,dim=dim(optimMatTrit));
for(o in 1:dim(optimMatBit)[2]){
  optimMatBitNorm[,o]=optimMatBit[,o]/a[o]
  optimMatTritNorm[,o]=optimMatTrit[,o]/a[o]
}

#subplots
#once the optimal reates were taken into account *by considering the minimum of the distance
#across all values that were simulated, we can plot the distance in the color axis and consider
#the number of scratchpads and the number of gnerations in the x and y axis.
#this plot is from c(0,1), 1->perfect reconstruction; 0->random guess

pdf(paste(pathPlots,"heatmaps_blue_BitvsTrit_atOptimalRates.pdf"))
vals=unique(scales::rescale(c(1-optimMatBitNorm)))
o<-order(vals,decreasing=F)
cols<-scales::col_numeric("Blues",domain=NULL)(vals)
colz<-setNames(data.frame(vals[o],cols[o]),NULL)

p1=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=1-optimMatBitNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))
p2=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=1-optimMatTritNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))

subplot(p1,p2)
dev.off()

x11()
plotList=list()
par(mfrow =c(2,3))
muVar = array(0,dim =c(length(barcodes),length(mus)))
for (ng in 1:length(generations)){
  for (muIdx in 1:length(mus)){
    for (bc in 1:length(barcodes)){
      muVar[bc,muIdx]= apply(muVariation[[muIdx]][['trit']][[bc]][[ng]],2,mean)[12]
    }
  }
  plotList[[ng]]=plot_ly(y=as.character(log(barcodes)),x=as.character(mus),z=muVar,type="heatmap")
  matplot(log(barcodes),muVar)
}
legend("right",as.character(mus),col=seq_len(length(mus)),cex=0.8,fill=seq_len(length(mus)),title="edit rate")

#matplot(mus,t(apply(muVarMatrix,c(1,2),mean)),type="o",pch=1:4)
#legend("right",distNames,col=seq_len(length(distNames)),cex=0.8,fill=seq_len(length(distNames)),title="Generations")





#this function calculates the area under the curve of the ecdf of the vector. In this case, the cumulative fraction
#of trees (as the plot in Memoir 1.0), that have a certain degree of accuracy (RF.dist)
auc.ecdf<-function(x){
  emp.cdf = ecdf(x)
  return(max(cumsum(emp.cdf(min(x):max(x)))))
}








#paper plot
# HHMI plots / retreat plots
#how many generations you can do with N barcodes
par(mfrow=c(2,2))
optimMatNorm= array(0, dim=dim(optimMatTrit) ); for(i in 1:dim(optimMatTrit)[2]){ optimMatNorm[,i] = (a[i]-optimMatTrit[,i])/a[i] }
optimMatBitNorm= array(0, dim=dim(optimMatTrit) ); for(i in 1:dim(optimMatBit)[2]){ optimMatBitNorm[,i] = (a[i]-optimMatBit[,i])/a[i] }
t=c(0.80,0.85,0.90,0.95)
#ylims=c(100,100,100,200
ylims=c(60,60,60,100)
legend.y.pos = c(90,90,90,180)
for(tt in 1:length(t)){


    threshold=t[tt]
    #bc.per.gen.bit = array();for(i in 1:dim(optimMatBitNorm)[2]){ bc.per.gen.bit[i] = min( which(optimMatBitNorm[,i]>threshold) ) }
    bc.per.gen.bit = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.bit[i] = min( which(optimMatNorm[,i]>t[3]) ) }
    bc.per.gen.bit2 = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.bit2[i] = min( which(optimMatNorm[,i]>t[1]) ) }
    bc.per.gen.trit = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.trit[i] = min( which(optimMatNorm[,i]>threshold) ) }

    plot(generations,barcodes[bc.per.gen.trit],type="l",lwd=2.5,
         xlim=c(2,10),ylim=c(2,ylims[tt]),cex.lab=1.5,cex.axis=1.2,ylab="Array length",xlab="Generations",
         main=paste(toString(t[tt])," % reconstructability",sep=""))
    grid (10,10, lty = 6, col = "cornsilk2")
    lines(generations,barcodes[bc.per.gen.bit],type="l",col="black",lwd=2)
    lines(generations,barcodes[bc.per.gen.trit],type="l",col="blue",lwd=2)
    lines(generations,barcodes[bc.per.gen.bit2],type="l",col="gray",lwd=2)
    #legend(3,legend.y.pos[tt],c("Trit","Bit"),fill=c("blue","black"))
}
