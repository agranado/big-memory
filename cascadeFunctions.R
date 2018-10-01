

# For a given number of generations and a given number of integrases (or recording mechanisms, gRNA etc)
# we need to indicate at which point in time (generation units), the next element will become active.
cascadeActivation<-function(nGen,nIntegrases){
    c( nGen %% nIntegrases,ceiling(nGen/nIntegrases),  nGen%/%nIntegrases)
    act_time =array()
    a =nGen %% nIntegrases
    b =ceiling(nGen/nIntegrases)
    c = nGen%/%nIntegrases

    if(a>0){
      for (i in 1:a){
          act_time =c(act_time,rep(i,b))
      }
      act_time = act_time[-1]
      act_time2=array()
      for(i in (a+1):nIntegrases){
          act_time2 =c(act_time2, rep(i,c))
      }
      act_time2 = act_time2[-1]
      act_time = c(act_time,act_time2)
    }else{

      for(i in 1:nIntegrases){
        act_time = c(act_time,  rep(i,c))
      }
      act_time = act_time[-1]
    }

    return(act_time)
  }

 partitionBarcodes <- function (fullBarcode,totalInts, currentInts){
   barcodeLength = length(fullBarcode)
   editableIndx =(barcodeLength / totalInts *(currentInts-1) +1) : (barcodeLength /totalInts * currentInts);
   a = fullBarcode[editableIndx]
   return(list(a,editableIndx))
}
library(phytools) #for pasting trees together

#this function works well for nIntegrases=2
#for more integrases, it will require recursive reconstruction
#We might not need more than 4 integrases for the paper (or in reality)
cascadeReconstruction<-function(barcodeLeaves,totalInts,currentInts,nGen,mu,alpha){

    # 1 Partition barcode into the different integrase-specific elements
      ##  barcodeLeaves= this.tree$tip.label
      #let's separate the barcode string into arrays, so we can index using the integrase indexes.
      # te=apply(as.array(barcodeLeaves),1,strsplit,"")
      # all.barcodes =as.array(unlist(te))
      # dim(all.barcodes)<-c(nchar(barcodeLeaves[1]) ,length(barcodeLeaves)) #Colums are the cells,
      # # all.barcoes[,1] == barcodesLeaves[1]

      #indexes from the function are wrong, so the first 1:n-1 elements are used then n + 1 :end
      # FIX THIS

      # for(t in 1){
      #     sub.barcodes= matrix(0,dim(all.barcodes)[2],length(partitionBarcodes(all.barcodes[,1],totalInts,1)))
      #     for(c in 1:dim(all.barcodes)[2]) { #for each cell
      #         sub.barcodes[c,]= partitionBarcodes(all.barcodes[,c],totalInts,t)
      #
      #     }
      #
      #   # 2 Reconstruct each subtree
      # }
      act_time = cascadeActivation(nGen,totalInts)

      #calcute the index then use substr
      #this function returns both the sub.barcode and the indexes
      editableIndx=partitionBarcodes(strsplit(barcodeLeaves[1],"")[[1]],totalInts,currentInts)[[2]]
    #  editableIndx= editableIndx[-1]
      # NEXT LINE is WRONG, fix it after correcting the indexes in the next dataset
      sub.barcodes = substr(barcodeLeaves,range(editableIndx)[1],range(editableIndx)[2])

      #generate a small subtree based on unique profiles
      unique.sub.barcodes = unique(sub.barcodes)

      sub.nGen = sum(act_time ==currentInts)
      matdist_ = manualDistML(unique.sub.barcodes,mu,alpha,sub.nGen)
      colnames(matdist_)<-unique.sub.barcodes
      manualTree_1 =upgma(as.dist(t(matdist_))) #now the tree has names
      #this works
      big.tree = manualTree_1
      # for each unique.sub.barcode, get all the daughters coming from that clone
      for(i in 1:length(unique.sub.barcodes)){
        daughter.index= grep(unique.sub.barcodes[i],sub.barcodes) #get all the cells that have same barcode for first integrase (theu belong to the same clone)
        daughters.barcodes = barcodeLeaves[daughter.index] #get the barcodes

        matdist_2 = manualDistML(daughters.barcodes,mu,alpha,sub.nGen) # submatrix for distances
        colnames(matdist_2)<-daughters.barcodes # name the matrix with the actual barcodes, will be used later for distance between the real and reconstructed trees
        manualTree_2 =upgma(as.dist(t(matdist_2))) #now the tree has names

        # AND REPEAT (recursively)
        # 3 Join subtrees into a big tree
        which.leaf=which(big.tree$tip.label==unique.sub.barcodes[i])
        big.tree$tip.label[which.leaf]<-"NA"
        manualTree_2$root.edge<-0
        big.tree = paste.tree(big.tree,manualTree_2)
      }



      #using phylo tools to join elements into a branch/leaf



        return(big.tree)

  }
