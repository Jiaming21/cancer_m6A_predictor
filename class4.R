electronIonInteraction <- function(Data,NI=3){
  for(i in 1:NI){
    if(i==1){
      NP <- c("A","G","C","T")
    }else{
      NP1 <- NULL
      for(j in c("A","G","C","T")){
        for(k in 1:length(NP)){
          NP1 <- c(NP1,paste0(NP[k],j))
        }
      }
      NP <- NP1
    }
  }
  for(i in 1:NI){
    if(i==1){
      NPE <- c(0.1260,0.0806,0.1340,0.1335)
    }else{
      NP1E <- NULL
      for(j in c(0.1260,0.0806,0.1340,0.1335)){
        for(k in 1:length(NPE)){
          NP1E <- c(NP1E,NPE[k]*j)
        }
      }
      NPE <- NP1E
    }
  }
  MA2 <- matrix(NA,ncol = length(NP),nrow = length(Data))
  colnames(MA2) <- NP
  for(i in 1:length(NP)){
    MA2[,i]<- vcountPattern(NP[i],DNAStringSet(Data))
  }
  M3 <- MA2
  for(i in 1:nrow(M3)){
    M3[i,] <-  MA2[i,]/ncol(MA2)
  }  
  M4 <- M3
  for(i in 1:ncol(M3)){
    M4[,i] <-  M3[,i]*NPE[i]
  }
  return(M4)
}