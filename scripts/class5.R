RNA_SCORE <- readRDS("./RNA_score.rds")


autoCovariance <- function(Data,NI=2,NTYPE="RNA",descriptor="Shift",LG=10,convert=2){
  MA <- NC1(Data)
  MAF <- MA
  if(convert==1){
    MAF[which(MA=="A")] <- "U"
    MAF[which(MA=="T")] <- "A"
    MAF[which(MA=="C")] <- "G"
    MAF[which(MA=="G")] <- "C"
  }
  MA <- MAF
  NTYPE="RNA"
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
  descriptor_score <- RNA_SCORE[,which(colnames(RNA_SCORE)==descriptor)]
  RNAname <- rownames(RNA_SCORE)
  # RNAname <- c("AA","TT","AT","TA","CA","GT","CT","GA","CG","GC","GG","CC","AG","AC","TG","TC")
  MA2 <- MA1
  for(i in 1:length(RNAname)){
    MA2[MA1==RNAname[i]] <- descriptor_score[i]
  }
  MA3 <- matrix(as.numeric(MA2), ncol = ncol(MA2),nrow = nrow(MA2),byrow = F) 
  MA3M <- rowMeans(MA3)
  XXXA <- matrix(NA,ncol =LG,nrow = nrow(MA3))
  for(lg in 1:LG){
    for(i in 1:(ncol(MA3)-lg)){
      if(i==1){
        XXX <- (MA3[,i]-MA3M)*(MA3[,(i+lg)]-MA3M)
      }else{
        XXX <- (MA3[,i]-MA3M)*(MA3[,(i+lg)]-MA3M) + XXX
      }
    }
    XXXA[,lg] <- XXX*(1/(ncol(MA3)-lg)) 
  }
  return(XXXA)
}


crossCovariance <- function(Data,NI=2,NTYPE="RNA",descriptor1="Shift",descriptor2="Slide",LG=10){
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
  MA1 <- matrix(NA,ncol = (ncol(MA)-NI+1),nrow = nrow(MA))
  for(i in 1:(ncol(MA)-NI+1)){
    for(j in 1:NI){
      if(j==1){
        MA1Z <- MA[,i]
      }else{
        MA1Z <- paste0(MA1Z, MA[,(i+j-1)])
      }
    }
    MA1[,i] <- MA1Z
  }
  descriptor_score <- RNA_SCORE[,which(colnames(RNA_SCORE)==descriptor1)]
  # RNAname <- rownames(RNA_SCORE)
  RNAname <- c("AA","TT","AT","TA","CA","GT","CT","GA","CG","GC","GG","CC","AG","AC","TG","TC")
  MA2 <- MA1
  for(i in 1:length(RNAname)){
    MA2[MA1==RNAname[i]] <- descriptor_score[i]
  }
  MA3 <- matrix(as.numeric(MA2), ncol = ncol(MA2),nrow = nrow(MA2),byrow = F) 
  MA3M <- rowMeans(MA3)
  
  descriptor_score <- RNA_SCORE[,which(colnames(RNA_SCORE)==descriptor2)]
  # RNAname <- rownames(RNA_SCORE)
  RNAname <- c("AA","TT","AT","TA","CA","GT","CT","GA","CG","GC","GG","CC","AG","AC","TG","TC")
  MA2 <- MA1
  for(i in 1:length(RNAname)){
    MA2[MA1==RNAname[i]] <- descriptor_score[i]
  }
  MA4 <- matrix(as.numeric(MA2), ncol = ncol(MA2),nrow = nrow(MA2),byrow = F) 
  MA4M <- rowMeans(MA4)
  
  XXXA <- matrix(NA,ncol =LG,nrow = nrow(MA3))
  for(lg in 1:LG){
    for(i in 1:(ncol(MA3)-lg)){
      if(i==1){
        XXX <- (MA3[,i]-MA3M)*(MA4[,(i+lg)]-MA4M)
      }else{
        XXX <- (MA3[,i]-MA3M)*(MA4[,(i+lg)]-MA4M) + XXX
      }
    }
    XXXA[,lg] <- XXX*(1/(ncol(MA3)-lg)) 
  }
  return(XXXA)
}
