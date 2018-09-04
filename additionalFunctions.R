

detectOS <-function(){
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
      return(list(pathName,pathName2))
}
