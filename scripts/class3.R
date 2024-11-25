PSNPprepare <- function(Data,NI=3,NTYPE="RNA"){
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
  
  MA2 <- matrix(NA,ncol = (ncol(MA)-NI+1),nrow = length(NP))
  rownames(MA2) <- NP
  for(i in 1:ncol(MA2)){
    for(j in 1:length(NP)){
      MA2[j,i] <- length(which(MA1[,i]==NP[j]))
    }
  }
  M3 <- MA2
  for(i in 1:ncol(MA2)){
    M3[,i] <- MA2[,i]/sum(MA2[,i])
  }
  return(M3)
}

PSNP <- function(Data1,Data2,NI=3,NTYPE="RNA"){
  MA1 <- PSNPprepare(Data1,NI,NTYPE)
  MA2 <- PSNPprepare(Data2,NI,NTYPE)
  MAE <- MA1-MA2
  XNAME <- rownames(MAE)
  Data3 <- rbind(NC1(Data1),NC1(Data2))
  MA1 <- matrix(NA,ncol = (ncol(Data3)-NI+1),nrow = nrow(Data3))
  for(i in 1:(ncol(Data3)-NI+1)){
    for(j in 1:NI){
      if(j==1){
        MA1Z <- Data3[,i]
      }else{
        MA1Z <- paste0(MA1Z, Data3[,(i+j-1)])
      }
    }
    MA1[,i] <- MA1Z
  }
  MA <- matrix(NA,ncol = ncol(MA1),nrow = nrow(MA1))
  for(i in 1:ncol(MA)){
    for(j in 1:length(XNAME)){
      MA[which(MA1[,i]==XNAME[j]),i] <- MAE[j,i]
    }
  }
  return(MA)
}

PSNPRE <- function(Data1,Data2,Data5,Data6,NI=3,NTYPE="RNA"){
  MA1 <- PSNPprepare(Data1,NI,NTYPE)
  MA2 <- PSNPprepare(Data2,NI,NTYPE)
  MAE <- MA1-MA2
  XNAME <- rownames(MAE)
  output <- list()
  Data3 <- rbind(NC1(Data1),NC1(Data2))
  MA1 <- matrix(NA,ncol = (ncol(Data3)-NI+1),nrow = nrow(Data3))
  for(i in 1:(ncol(Data3)-NI+1)){
    for(j in 1:NI){
      if(j==1){
        MA1Z <- Data3[,i]
      }else{
        MA1Z <- paste0(MA1Z, Data3[,(i+j-1)])
      }
    }
    MA1[,i] <- MA1Z
  }
  MA <- matrix(NA,ncol = ncol(MA1),nrow = nrow(MA1))
  for(i in 1:ncol(MA)){
    for(j in 1:length(XNAME)){
      MA[which(MA1[,i]==XNAME[j]),i] <- MAE[j,i]
    }
  }
  
  output[[1]] <- MA
  
  Data3 <- rbind(NC1(Data5),NC1(Data6))
  MA1 <- matrix(NA,ncol = (ncol(Data3)-NI+1),nrow = nrow(Data3))
  for(i in 1:(ncol(Data3)-NI+1)){
    for(j in 1:NI){
      if(j==1){
        MA1Z <- Data3[,i]
      }else{
        MA1Z <- paste0(MA1Z, Data3[,(i+j-1)])
      }
    }
    MA1[,i] <- MA1Z
  }
  MA <- matrix(NA,ncol = ncol(MA1),nrow = nrow(MA1))
  for(i in 1:ncol(MA)){
    for(j in 1:length(XNAME)){
      MA[which(MA1[,i]==XNAME[j]),i] <- MAE[j,i]
    }
  }
  
  output[[2]] <- MA
  return(output)
}

# Notmethylated_sample <- readRDS("C:/Users/kunqi.chen/Documents/Robert Dan/PHD/5/PesU/method/exonrealNegativeT.rds")
# mature_sample <- readRDS("C:/Users/kunqi.chen/Documents/Robert Dan/PHD/5/PesU/method/exonmRNA_CorPseU.rds")
# string1 <- as.character(DNAStringSet(Views(Hsapiens,Notmethylated_sample+20)))
# string2 <- as.character(DNAStringSet(Views(Hsapiens,mature_sample+20)))
# 
# PSNP(string1,string2,NTYPE="DNA")
