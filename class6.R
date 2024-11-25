PseKNC <- function(Data,NI=2,NTYPE="RNA"){
  MA <- NC1(Data) 
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  for(i in 1:NI){
    if(i==1){
      NP <- c("A","G","C",U)
    }else{
      NP1 <- NULL
      for(j in c("A","G","C",U)){
        for(k in 1:length(NP)){
          NP1 <- c(NP1,paste0(NP[k],j))
        }
      }
      NP <- NP1
    }
  }
  MA.A <- MA[,-((ncol(MA)+1)/2)]
  MA1 <- matrix(NA,ncol = (ncol(MA.A)/2),nrow = nrow(MA))
  for(i in 1:(ncol(MA.A)/2)){
    for(j in 1:NI){
      if(j==1){
        MA1Z <- MA.A[,(2*(i-1)+1)]
      }else{
        MA1Z <- paste0(MA1Z, MA.A[,(2*i)])
      }
    }
    MA1[,i] <- MA1Z
  }
  MA2 <- matrix(NA,ncol = length(NP),nrow = nrow(MA1))
  for(i in 1:length(NP)){
    MA2[,i] <- rowCumsums(matrix(as.numeric(MA1 == NP[i] ), ncol = ncol(MA1), byrow = F))[,ncol(MA1)]
  }
 colnames(MA2) <- NP
  return(MA2)
}