library(dplyr)
ONE_HOT <- function(Data){
  MA <- NC1(Data)
  MA2 <- matrix(NA,ncol =(4*ncol(MA)),nrow = nrow(MA))
  A <- c(1,0,0,0)%>%matrix(nrow=1)
  C <- c(0,1,0,0)%>%matrix(nrow=1)
  G <- c(0,0,1,0)%>%matrix(nrow=1)
  T <- c(0,0,0,1)%>%matrix(nrow=1)
  U <- c(0,0,0,1)%>%matrix(nrow=1)
  for(i in 1:ncol(MA)){
print(i)
    i1 <- (i-1)*4 +1
    i2 <- (i-1)*4 +2
    i3 <- (i-1)*4 +3
    i4 <- i*4
    for(j in c("A","C","G","T","U")){
      MA2[which(MA[,i]==j),i1] <- eval(parse(text=j))[1]
      MA2[which(MA[,i]==j),i2] <- eval(parse(text=j))[2]
      MA2[which(MA[,i]==j),i3] <- eval(parse(text=j))[3]
      MA2[which(MA[,i]==j),i4] <- eval(parse(text=j))[4]
      # print(eval(parse(text=j)))
    }
  }
  return(MA2)
}

ChemicalProperty <- function(Data){
  MA <- NC1(Data) #将data转成字符串矩阵
  MA2 <- matrix(NA,ncol =(3*ncol(MA)),nrow = nrow(MA)) #生成一个比data多3倍列的矩阵
  A <- c(1,1,1)%>%matrix(nrow=1)# %>%把左侧准备的数据或表达式，传递给右侧的函数调用或表达式进行运行，可以连续操作就像一个链条一样
  C <- c(0,1,0)%>%matrix(nrow=1)
  G <- c(1,0,0)%>%matrix(nrow=1)
  T <- c(0,0,1)%>%matrix(nrow=1)
  U <- c(0,0,1)%>%matrix(nrow=1)
  N <- c(0,0,0)%>%matrix(nrow=1)
  for(i in 1:ncol(MA)){
    i1 <- (i-1)*3 +1 # 1, 4, 7...
    i2 <- (i-1)*3 +2 # 2, 5, 8
    i3 <- i*3 # 3,6,9
    for(j in c("A","C","G","T","U")){
      if(length(which(MA[,i]==j))>0){
        MA2[which(MA[,i]==j),i1] <- eval(parse(text=j))[1] #parse函数解析字符串为表达式；eval函数执行表达式输出结果,等同于 A[1]/C[1]/G[1]
        MA2[which(MA[,i]==j),i2] <- eval(parse(text=j))[2]
        MA2[which(MA[,i]==j),i3] <- eval(parse(text=j))[3] #可以理解为将MA中相应的碱基（ATCGNT）转换为01编码，每个碱基对应3个数字编码
      }
      # print(eval(parse(text=j)))
    }
  }
  return(MA2)
}