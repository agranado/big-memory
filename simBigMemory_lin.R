
#commands that execute this function in parallel

#200 repeats, using the indicated parameters,
# the list of methods contains the distance metrics that we want to test in the stringdist function
# >  results= foreach(i=1:200) %dopar% simMemoirStrdist(nGen,mu,alpha,barcodeLength,methods)

#convert list to matrix
# >  results.matrix=do.call(rbind,results)
# >  apply(results.matrix,2,mean)

#Arguments of the function:
#barcodeLength<-16
#nGen =5
#mu= (1/barcodeLength)*3;
#mu=1/5
#alpha= 2/3;

library(phangorn)
library(stringdist)
library(doParallel)
library(gplots)
source("../lineageSoftware/MLfunctions.R") #load from the lineageSoftware repository
source("simulation3_lin.R")
source("additionalFunctions.R")


library(gplots)


rand.dist<-c(10,  26,  58, 120, 250, 506)

compareDist <- function(simulationType='trit',nGen=4,mu=0.3,alpha_=2/3,barcodeLength=20,nRepeats=20,methods=c('osa','lv','dl','hamming','lcs','qgram','cosine','jaccard','jw','soundex'),recType="integrase"){


  results= foreach(i=1:nRepeats) %dopar% simMemoirStrdist(nGen=nGen,mu=mu,alpha=alpha_,barcodeLength=barcodeLength,methods=methods,simulationType=simulationType)

 #let's unlist the results from the parallel calculations:
  results_=list()
  tree.list = list()
  for(i in 1:length(results)){
    results_[[i]] =results[[i]][[1]]
    tree.list[[i]] = results[[i]][[2]]
  } #this takes the first element of the results list which is the array with distance calculations

  #save the simulated trees to the hard drive
  simulation.file = paste("recType_",recType,"_mu_",toString(mu),
     "_BC_",toString(barcodeLength),
     "_nG_",toString(nGen),"_Nrep_",toString(nRepeats),"_.rda",sep="")



  results.matrix=do.call(rbind,results_)
  #Optional when only interested in the mean
  #apply(results.matrix,2,mean)
  return(results.matrix)



}

#functions to use the proportion of perfect trees as the measure.
#May 8
eq.zero<-function(r,x){sum(r[,x]==0)}



#nGen=3;mu=0.4;alpha=1/2;barcodeLength=10;methods=c();simulationType='trit';
#April 8th
#Test stringdistance measures using the stringdist R library
#use the same format as before but testing different methods included in the stringdist function
simMemoirStrdist<-function(nGen=3,mu=0.4,alpha=1/2,barcodeLength=10,methods=c(),simulationType='trit',recType="integrase"){
  #load necessary libraries and functions
  #detection of OS

  #update (trying to create a single branch that works in AWS and in my laptop)
  os=system("cat ../os.txt",intern = T)
  if(os=="mac"){ #local Alejandro's laptop
    pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
    pathName2="/Users/alejandrog/MEGA/Caltech/trees/simulation"
  }else if(os=="linux"){ #AWS server or any other linux machine (path is for AWS)
    pathName = "/home/alejandrog/MEGA/Caltech/lineage/"
    pathName2= "/home/alejandrog/MEGA/Caltech/lineage"

  }else if(os=="linux_local"){
    #linux desktop
    file.dir = "/home/alejandrog/MEGA/Caltech/trees/GraceData/integrase-data/10mer/"
    setwd("/home/alejandrog/MEGA/Caltech/trees/macbookBranch/lineageSoftware")
    pathName2= "/home/alejandrog/MEGA/Caltech/lineage"
    pathName= "/home/alejandrog/MEGA/Caltech/lineage/"
    #At this point we are supposed to be in the local branch already, because we read the os.txt...
  }
  #clear the variable (since it behaves as global)
  if(exists("firstCell")){
    rm(firstCell)
  }

 firstCell<-Node$new("1"); firstCell$barcode <-paste(rep("u",barcodeLength),collapse="");
  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

# SIMULATION
# # # # # # # # #
 # # # # # # # #
# # # # # # # # #
  actIntegrase=1 # start recording using integrase 1
  nIntegrases =4 # how many integrases comprise the cascade
  act_times = c(1,1,2,2,3,3,3,4,4,4)
  for (g in 1:nGen){
    #this function simulates one generation of the tree
    actIntegrase = act_times[g]
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType,recType=recType,actIntegrase,nIntegrases)
  }
# # # # # # # # #
 # # # # # # # # #
# # # # # # # # #


  #prints only the barcodes for all leaves
  #  print(firstCell,"barcode")
  #  print("Tree simulation completed")
  #save to file as newick tree
  # #save the length of branches plus the ID (which so far is a number)
  # newickTree<-ToNewick(firstCell)
  # #Generate unique ID for writing file to disc (we'll erase it later)
  # fileID = toString(runif(1))
  #
  # #firstCellFile = paste(pathName,"trees/firstCell",fileID,".nwk",sep="")
  #
  # firstCellFile =tempfile("trees/firstCell",tmpdir = pathName2)
  # firstCellFile =paste(firstCellFile,fileID,".nwk",sep="")
  #
  #
  # write(newickTree,file=firstCellFile)
  # #load the tree from the file as a tree structure.
  # trueTree<-read.tree(file=firstCellFile)

  #file is now deleted
  #  print("True tree read")
  # plot(trueTree,main=paste("True tree ",sep=""))

  trueTree<-as.phylo.Node(firstCell)


  #get the sequences from the simulated tree + names
  barcodes<-firstCell$Get("barcode")

  #now we have the patters for ALL cells, but we need only the leaves, since that is what
  #we are going to use for reconstruction.
  #The way the tree was built, only the last 2^g cells are leaves; where g is the number of generations
  #take the number ID for the leaves.
  leavesID<-(length(barcodes)-2^nGen+1):length(barcodes)
  #grab those cells from the tree

  barcodeLeaves = array()
  namesLeaves=array()
  for (l in 1:length(leavesID)){
    barcodeLeaves[l] <-barcodes[names(barcodes)==leavesID[l]]
    namesLeaves[l] <-names(barcodes)[names(barcodes)==leavesID[l]]
  }

  #take the barcodes from internal nodes, not useful here but for saving the complete tree simulations
  #
  nodesID=1:(leavesID[1]-1)

  barcodeNodes=array()
  namesNodes=array()
  for(l in 1:length(nodesID)){
    barcodeNodes[l]<-barcodes[names(barcodes)==nodesID[l]]
    namesNodes[l]<-names(barcodes)[names(barcodes)==nodesID[l]]
  }
  names(barcodeNodes)<-namesNodes
  ### NOTE names for the new tree are not correct 
  names(barcodeNodes)<-namesNodes[as.numeric(trueTree$node.label)]
  #names(barcodeLeaves)<-namesLeaves
  #now barcodeLeaves has all the leaves of the tree, with their original ID from the data.tree structure.
  #create Fasta file using only the leaves of the tree (n= 2^g)
  fastaBarcodes<-convertSimToFasta(barcodeLeaves)
  #convert name of variable to string

  #3  varName<-deparse(substitute(firstCell))

  fasID =toString(runif(1))
  #fasIN <-paste(pathName,"fasta/firstCell_",fasID,".fas",sep="")
  fasIN =tempfile("fasta/firstCell",tmpdir = pathName2)
  fasIN = paste(fasIN,fasID,".fas",sep="")

  #randomize the barcodes such that order is not a factor in the lineage reconstruction
  #the real tree (because the way is constructed, has an order 1,2,3,...N), if the barcodes are not
  #randomized, the order will "help" to the reconstruction which is not good!
  rand.barcode.order<-sample(1:length(fastaBarcodes))

  #We re-order the barcodes using a fixed (but random) order
  fastaBarcodes<-fastaBarcodes[rand.barcode.order]
  barcodeLeaves<-barcodeLeaves[rand.barcode.order]
  write(fastaBarcodes,file=fasIN)



  sed.return<-convertMemoirToDNA(fasIN)
  #now we can use the phyDat
  if(sed.return){
    memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
    #for distance based trees
    #from the phangorn tutorial pdf:
  }
  #Apr 8th: this is where the distance comes into place

  #Apr 8th:
  #we calculate string distances only for the leaves ( which is the data we actually get)

  allDistances = array()
  # for(m in 1:length(methods)){
  #   stringTree= upgma(stringdistmatrix(barcodeLeaves,method=methods[m]))
  #   stringTree$tip.label<-trueTree$tip.label
  #   allDistances[m]= RF.dist(stringTree,trueTree)
  # }

  #control against default dist.ml function + UPGMA, which so far is the best method
  m=0
  if(sed.return==1){
    dm<-dist.ml(memoirfas)
    dm.ham=dist.hamming(memoirfas)
    if(simulationType=='binary'){
      treeUPGMA<-upgma(dm.ham)
    }else{
      treeUPGMA<-upgma(dm)
    }
    treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
    allDistances[m+1]= RF.dist(treeUPGMA_noSeq,trueTree)
  }else{
    allDistances[m+1]= NaN
  }


  matdist_ = manualDistML(barcodeLeaves,mu,alpha,nGen)
  manualTree_ =upgma(as.dist(t(matdist_)))
  manualTree_$tip.label= treeUPGMA$tip.label

  allDistances[m+2]= RF.dist(removeSeqLabel(manualTree_),trueTree)

  #try new distance using the built in dendrogram of heatmap2

  h=heatmap.2(matdist_+t(matdist_),trace="none",dendrogram = 'column')
  # h=heatmap.2(matdist_+t(matdist_),Colv="Rowv")
  # heatmap.tree=as.phylo(as.hclust(h$colDendrogram))
  # heatmap.tree$tip.label = treeUPGMA$tip.label
  # allDistances[m+3]= RF.dist(removeSeqLabel(heatmap.tree),trueTree)
  #

  #alternative w/o plotting the actual heatmap, only hclust method
  hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
  hclust.tree$tip.label = treeUPGMA$tip.label
  allDistances[m+3]= RF.dist(removeSeqLabel(hclust.tree),trueTree)

  #system(paste("rm ",firstCellFile,sep=""))
  if(sed.return){
    system(paste("rm ",paste(fasIN,".bak",sep=""),sep=""))
    system(paste("rm ",fasIN,sep=""))
  }
  #system(paste("rm ",fasIN,".bak",sep=""))


  return(list(allDistances,trueTree))
}
#END of simulation function
# # # # # # # # # #
 # # # # # # # # #
# # # # # # # # # #




runThisScript <- function (){
  #run memoirSim.R once and plot both trees.
  x11()
  par(mfrow=c(1,2))
  source("/Users/alejandrog/MEGA/Caltech/trees/simulation/memoirSim.R")
  ratioMat = manualDist(barcodeLeaves,mu,alpha,nGen );ratioTree = upgma(as.dist(t(ratioMat))); ratioTree$tip.label= treeUPGMA$tip.label;

  plot(ratioTree,main=paste("Manual tree",toString(RF.dist(trueTree,removeSeqLabel(ratioTree))) ))
  plot(treeUPGMA,main=paste("UPGMA dist tree",toString(RF.dist(trueTree,removeSeqLabel(treeUPGMA))) ))

  #execute multiple trees to get statistics
  # results= foreach(i=1:200) %dopar% simMemoirStrdist(3,0.3,alpha,barcodeLength,methods); results.matrix=do.call(rbind,results); meanDist =apply(results.matrix,2,mean)
}



#run function to calculate unique states

simMemoirBarcodes<-function(nGen=3,mu=0.4,alpha=1/2,barcodeLength=10,methods=c(),simulationType='trit'){
  #load necessary libraries and functions
  #detection of OS

  #update (trying to create a single branch that works in AWS and in my laptop)
  os=system("cat ../os.txt",intern = T)
  if(os=="mac"){ #local Alejandro's laptop
    pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
    pathName2="/Users/alejandrog/MEGA/Caltech/trees/simulation"
  }else if(os=="linux"){ #AWS server or any other linux machine (path is for AWS)
    pathName="/home/ubuntu/alejandrog/Caltech/lineage/"
    pathName2="/home/ubuntu/alejandrog/Caltech/lineage"

  }
  #clear the variable (since it behaves as global)
  if(exists("firstCell")){
    rm(firstCell)
  }
  #create intialize the tree using 1 as the ID
  #Node is a function from data.tree
  firstCell <- Node$new("1")

  #number of letters (this will change for simulation)
  #A goal is to understand how the lenght of the pattern affects the ability to reconstruc lineages

  #intialitze the barcode using "u->unchanged
  firstCell$barcode <-paste(rep("u",barcodeLength),collapse="")

  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

  for (g in 1:nGen){
    #this function simulates one generation of the tree
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType)
  }

  #prints only the barcodes for all leaves
  #  print(firstCell,"barcode")
  #  print("Tree simulation completed")
  #save to file as newick tree
  # #save the length of branches plus the ID (which so far is a number)
  # newickTree<-ToNewick(firstCell)
  # #Generate unique ID for writing file to disc (we'll erase it later)
  # fileID = toString(runif(1))
  #
  # #firstCellFile = paste(pathName,"trees/firstCell",fileID,".nwk",sep="")
  #
  # firstCellFile =tempfile("trees/firstCell",tmpdir = pathName2)
  # firstCellFile =paste(firstCellFile,fileID,".nwk",sep="")
  #
  #
  # write(newickTree,file=firstCellFile)
  # #load the tree from the file as a tree structure.
  # trueTree<-read.tree(file=firstCellFile)

  #file is now deleted
  #  print("True tree read")
  trueTree<-as.phylo.Node(firstCell)
  # plot(trueTree,main=paste("True tree ",sep=""))

  #get the sequences from the simulated tree + names
  barcodes<-firstCell$Get("barcode")

  #now we have the patters for ALL cells, but we need only the leaves, since that is what
  #we are going to use for reconstruction.
  #The way the tree was built, only the last 2^g cells are leaves; where g is the number of generations
  #take the number ID for the leaves.
  leavesID<-(length(barcodes)-2^nGen+1):length(barcodes)
  #grab those cells from the tree

  barcodeLeaves = array()
  namesLeaves=array()
  for (l in 1:length(leavesID)){
    barcodeLeaves[l] <-barcodes[names(barcodes)==leavesID[l]]
    namesLeaves[l] <-names(barcodes)[names(barcodes)==leavesID[l]]
  }
  names(barcodeLeaves)<-namesLeaves
  #now barcodeLeaves has all the leaves of the tree, with their original ID from the data.tree structure.
  #create Fasta file using only the leaves of the tree (n= 2^g)
  fastaBarcodes<-convertSimToFasta(barcodeLeaves)
  #convert name of variable to string

  #3  varName<-deparse(substitute(firstCell))

  fasID =toString(runif(1))
  #fasIN <-paste(pathName,"fasta/firstCell_",fasID,".fas",sep="")
  fasIN =tempfile("fasta/firstCell",tmpdir = pathName2)
  fasIN = paste(fasIN,fasID,".fas",sep="")

  #randomize the barcodes such that order is not a factor in the lineage reconstruction
  #the real tree (because the way is constructed, has an order 1,2,3,...N), if the barcodes are not
  #randomized, the order will "help" to the reconstruction which is not good!
  rand.barcode.order<-sample(1:length(fastaBarcodes))

  #We re-order the barcodes using a fixed (but random) order
  fastaBarcodes<-fastaBarcodes[rand.barcode.order]
  barcodeLeaves<-barcodeLeaves[rand.barcode.order]
  # write(fastaBarcodes,file=fasIN)
  #  print("writting fasta file, simulated tree")
  #this is the format:
  #>1_uuuuuu
  #uuuuuu
  #>2_uxuxuu
  #uxuxuu
  #>4_uxuxxc
  #uxuxxc
  #>8_uxuxxc
  #uxuxxc

  #we can also convert the data.tree format to igraph format
  #igraph has more algorithms and plotting function, so it's worth.
  #LIBRARY
  #library(igraph)
  #as.igraph.Node converts data.tree to igraph object
  #  firstCell_igraph<-as.igraph.Node(firstCell)
  #rename all vertices in the network using the barcodes.
  # V(firstCell_igraph)$name<-barcodes
  #the object has now all info.
  #REMEBER: this is the REAL tree, still we need to perform alignment and reconstruction


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  #### FROM here comes the alignment and reconstruction
  #previously saved fasta file
  #LIBRARY\


  #this a new type of object, from the package phangorn "phyDat"

  #we can apply lineage reconstruction here so later we can compare with the real tree.
  #before loading the fasta file we need to convert the characters to DNA, so that it is compatible
  #sed.return ==1 if replacement went through (Apr 25)
  # sed.return<-convertMemoirToDNA(fasIN)
  # #now we can use the phyDat
  # if(sed.return){
  #   memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
  #   #for distance based trees
  #   #from the phangorn tutorial pdf:
  # }
  # #Apr 8th: this is where the distance comes into place
  #
  # #Apr 8th:
  # #we calculate string distances only for the leaves ( which is the data we actually get)
  #
  # allDistances = array()
  # # for(m in 1:length(methods)){
  # #   stringTree= upgma(stringdistmatrix(barcodeLeaves,method=methods[m]))
  # #   stringTree$tip.label<-trueTree$tip.label
  # #   allDistances[m]= RF.dist(stringTree,trueTree)
  # # }
  #
  # #control against default dist.ml function + UPGMA, which so far is the best method
  # m=0
  # if(sed.return==1){
  #   dm<-dist.ml(memoirfas)
  #   dm.ham=dist.hamming(memoirfas)
  #   if(simulationType=='binary'){
  #     treeUPGMA<-upgma(dm.ham)
  #   }else{
  #     treeUPGMA<-upgma(dm)
  #   }
  #   treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
  #   allDistances[m+1]= RF.dist(treeUPGMA_noSeq,trueTree)
  # }else{
  #   allDistances[m+1]= NaN
  # }
  #
  # #Apr 9th
  # #Manual distance calculation (v beta1.0)
  # #  matdist = manualDist(barcodeLeaves,mu,alpha,nGen)
  # #  manualTree =upgma(as.dist(t(matdist)))
  # #  manualTree$tip.label= treeUPGMA$tip.label
  #
  # #  allDistances[m+2]= RF.dist(removeSeqLabel(manualTree),trueTree)
  #
  # # allDistances[m+2]=0
  #
  # #manualdist from MLfunctinos.R
  #
  # matdist_ = manualDistML(barcodeLeaves,mu,alpha,nGen)
  # manualTree_ =upgma(as.dist(t(matdist_)))
  # manualTree_$tip.label= treeUPGMA$tip.label
  #
  # allDistances[m+2]= RF.dist(removeSeqLabel(manualTree_),trueTree)
  #
  # #try new distance using the built in dendrogram of heatmap2
  #
  # h=heatmap.2(matdist_+t(matdist_),trace="none",dendrogram = 'column')
  # # h=heatmap.2(matdist_+t(matdist_),Colv="Rowv")
  # # heatmap.tree=as.phylo(as.hclust(h$colDendrogram))
  # # heatmap.tree$tip.label = treeUPGMA$tip.label
  # # allDistances[m+3]= RF.dist(removeSeqLabel(heatmap.tree),trueTree)
  # #
  #
  # #alternative w/o plotting the actual heatmap, only hclust method
  # hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
  # hclust.tree$tip.label = treeUPGMA$tip.label
  # allDistances[m+3]= RF.dist(removeSeqLabel(hclust.tree),trueTree)
  #
  # # print("All distances calcualted")
  #
  # # allDistances[m+5] = calcDstRF(as(removeSeqLabel(treeUPGMA),'TreeMan'),as(trueTree,'TreeMan'))
  # # allDistances[m+6] = calcDstRF(as(removeSeqLabel(manualTree_),'TreeMan'),as(trueTree,'TreeMan'))
  # # allDistances[m+7] = calcDstRF(as(removeSeqLabel(heatmap.tree),'TreeMan'),as(trueTree,'TreeMan'))
  # # allDistances[m+8] = calcDstRF(as(removeSeqLabel(hclust.tree),'TreeMan'),as(trueTree,'TreeMan'))
  # # allDistances
  #
  # #delete files

  system(paste("rm ",firstCellFile,sep=""))
  # if(sed.return){
  #   system(paste("rm ",paste(fasIN,".bak",sep=""),sep=""))
  #   system(paste("rm ",fasIN,sep=""))
  # }
  #system(paste("rm ",fasIN,".bak",sep=""))


  return(barcodeLeaves)
}



#Apr 23
#GIT update for analyzing performance of reconstruction method with simulated data.
#these updates happened before the comparison between binary and tri simulations
# code section for ATOM execution and examples
#registerDoParallel()
#results<-compareDist()
