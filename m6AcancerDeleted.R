################################################################################
##############################      加载R包      ###############################
################################################################################
library(e1071) #支持向量机（SVM）模型
library(ROCR) #R接收者操作特征（ROC）曲线
library(pROC) #机器学习算准确度要用的包
library(dplyr) #特征选择，数据框操作
library(m6ALogisticModel) #独立开发的机器学习包
library(SummarizedExperiment) #基因组数据处理
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #基因组数据处理
library(BSgenome.Hsapiens.UCSC.hg19) #基因组数据处理
library(MLmetrics) #模型评估
library(caret) #交叉验证
#构建特征要用的包：
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(phastCons100way.UCSC.hg19)
library(fitCons.UCSC.hg19)
#绘制箱线图
library(ggplot2)
################################################################################
###########################     设置工作目录      ##############################
################################################################################
setwd('/Users/jiaming/Desktop/cancer_m6A_predictor/file_m6AcancerDeleted')
################################################################################
############################      数据分析      ################################
################################################################################
#读取原始数据
data <- readRDS('./m6A_hg19.rds')
#取出癌症有关信息列
data <- data@elementMetadata[,-c(1:3)]
#修改非法字符
names(data) <- gsub('-','_',names(data))
#统计每个m6A位点在癌细胞系里出现次数
row_sums <- apply(data, 1, sum)
#给data加上统计次数
data$row_sums <- row_sums
#加回位点信息
dataCopy <- readRDS('./m6A_hg19.rds')
data <- GRanges(seqnames=dataCopy@seqnames,ranges=dataCopy@ranges,strand=dataCopy@strand,data)
#data按照统计次数(患癌指数)重排
data <- data[order(data$row_sums),]
#结果：0 -> 33920 (负样本) 采用 3392 (正样本)
#去除统计次数列
data <- data[,-26]
################################################################################
########################      挑出正负样本数据      ############################
################################################################################
#取data前33920行作为负样本，倒数3392行作为正样本
Nsample <- data[c(1:33920),]#33920行
Psample <- data[c(130692:134083),]#3392行
#去除无用数据
rm(data)
rm(dataCopy)
rm(row_sums)
################################################################################
############################      搭建特征      ################################
################################################################################
#genomic feature
#定义GFgenreation_m6A函数
GFgenreation_m6A <- function(data){
  
  analysis_data <- data
  matureSE <- SummarizedExperiment()
  rowRanges(matureSE) <- analysis_data 
  
  Additional_features_hg19 = list(
    HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
    YTHDC1_TREW = YTHDC1_TREW_gr,
    YTHDF1_TREW = YTHDF1_TREW_gr,
    YTHDF2_TREW = YTHDF2_TREW_gr,
    miR_targeted_genes = miR_targeted_genes_grl,
    TargetScan = TargetScan_hg19_gr,
    Verified_miRtargets = verified_targets_gr,
    METTL3_TREW = METTL3_TREW,
    METTL14_TREW = METTL14_TREW,
    WTAP_TREW = WTAP_TREW,
    METTL16_CLIP = METTL16_CLIP,
    ALKBH5_PARCLIP = ALKBH5_PARCLIP,
    FTO_CLIP = FTO_CLIP,
    FTO_eCLIP = FTO_eCLIP
  )
  
  data_standardized <- predictors_annot(se = matureSE,
                                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        bsgnm = Hsapiens,
                                        fc = fitCons.UCSC.hg19,
                                        pc = phastCons100way.UCSC.hg19,
                                        struct_hybridize = Struc_hg19,
                                        feature_lst = Additional_features_hg19,
                                        hk_genes_list = HK_hg19_eids,
                                        motif = c("DRACH"),
                                        motif_clustering = "A",
                                        isoform_ambiguity_method = "longest_tx",
                                        genes_ambiguity_method = "average",
                                        annot_clustering = matureSE,
                                        standardization = T) 
  
  
  GF <- mcols(data_standardized)
  return(GF)
}

GF_Nsample <- GFgenreation_m6A(Nsample)
GF_Psample <- GFgenreation_m6A(Psample)

saveRDS(GF_Nsample,'GF_Nsample.rds')
saveRDS(GF_Psample,'GF_Psample.rds')

GF_Nsample <- readRDS('./GF_Nsample.rds')
GF_Psample <- readRDS('./GF_Psample.rds')

#sequence feature
SeqFgeneration <- function(data,GF){ #GF多余
  UntestVarSeq <- data
  source("../class1.R")
  source("../class2.R")
  source("../class3.R")
  source("../class4.R")
  source("../class5.R")
  source("../class6.R")
  CP <- ChemicalProperty(UntestVarSeq)
  return(CP)
}
SeqFgeneration2 <- function(data,GF){
  UntestVarSeq <- data
  source("../class1.R")
  source("../class2.R")
  source("../class3.R")
  source("../class4.R")
  source("../class5.R")
  source("../class6.R")
  NF <- sequenceFeatures(UntestVarSeq,NTYPE="DNA")
  return(NF)
}
#用Nsample和GF_Nsample提取序列特征
Nsample$reference_sequence <- DNAStringSet(Views(Hsapiens,Nsample + 20))
sequenceFeature1 <- SeqFgeneration(as.character(Nsample$reference_sequence),GF_Nsample)
sequenceFeature2 <- SeqFgeneration2(as.character(Nsample$reference_sequence),GF_Nsample)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
NFeatureBoth <- cbind(GF_Nsample,sequenceFeature)
#用Psample和GF_Psample提取序列特征
Psample$reference_sequence <- DNAStringSet(Views(Hsapiens,Psample + 20))
sequenceFeature1 <- SeqFgeneration(as.character(Psample$reference_sequence),GF_Psample)
sequenceFeature2 <- SeqFgeneration2(as.character(Psample$reference_sequence),GF_Psample)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
PFeatureBoth <- cbind(GF_Psample,sequenceFeature)

#更新Nsample与Psample
Nsample <- NFeatureBoth
Psample <- PFeatureBoth
Nsample <- Nsample[,-which(colnames(Nsample) %in% c('long_exon','WTAP_TREW','last_exon','last_exon_400bp','exon_stop','constitutive_exon','last_exon_sc400','Stop_codons','UTR3','HK_genes','miR_targeted_genes'))] 
saveRDS(Nsample,'./NFeatureBoth.rds') #feature selection 要用
Psample <- Psample[,-which(colnames(Psample) %in% c('long_exon','WTAP_TREW','last_exon','last_exon_400bp','exon_stop','constitutive_exon','last_exon_sc400','Stop_codons','UTR3','HK_genes','miR_targeted_genes'))] 

################################################################################
#############################    搭建结果表    #################################
################################################################################
final_result <- as.data.frame(matrix(data = NA,10,10))
rownames(final_result) <- c("AUROC_5fcv","Sn_5fcv","Sp_5fcv","ACC_5fcv","MCC_5fcv","AUROC","Sn","Sp","ACC","MCC")
################################################################################
#############################     训练模型     #################################
################################################################################
for(i in 1:10) {
  print(paste0('round_',i,'_starts'))
  #将正样本80% (2714) 作为训练集，20% (678) 作为测试集
  set.seed(i)
  train_P_indx <- sample(1:nrow(Psample),2714)
  train_P <- Psample[train_P_indx,]
  test_P_indx <- setdiff(1:nrow(Psample),train_P_indx)
  test_P <- Psample[test_P_indx,]
  #将负样本取2714作为训练集，678作为测试集
  set.seed(i)
  train_N_indx <- sample(1:nrow(Nsample),2714)
  train_N <- Nsample[train_N_indx,]
  test_N_indx_pool<- setdiff(1:nrow(Nsample),train_N_indx)
  test_N_indx <- sample(test_N_indx_pool,678) 
  test_N <- Nsample[test_N_indx,]
  #将train_P和train_N合并起来，打上标签
  train_data <- rbind(train_P,train_N)
  label_train <- c(rep("positive",nrow(train_P)),rep("unmodified",nrow(train_N)))
  label_train <- factor(label_train,labels=c("positive","unmodified"))
  train_data$label <- label_train
  saveRDS(train_data,'./train_data.rds') #feature selection 要用
  #将test_P和test_N合并起来，打上标签
  test_data <- rbind(test_P,test_N)
  label_test <- c(rep(1,(nrow(test_P))),rep(0,(nrow(test_N))))
  test_data$label <- label_test
  #超参设定
  fitControl <- trainControl(method = "cv",  # 使用的验证类型        
                             number = 5,     # 折叠数 
                             savePred=TRUE,  # 在训练过程中是否应保存预测
                             summaryFunction = twoClassSummary, # 指定了一个应用于总结交叉验证过程结果的函数
                             classProbs = TRUE) # 指定是否应计算并返回类概率
  conservation <- train(label ~ ., data = train_data, 
                        method = "svmRadial", # SVM核函数
                        preProc = c("center", "scale"), # 代表自变量预处理方法,the data are centered and scaled
                        #数据的中心化是指数据集中的各项数据减去数据集的均值
                        #标准化是指中心化之后的数据在除以数据集的标准差，即数据集中的各项数据减去数据集的均值再除以数据集的标准差
                        trControl = fitControl) # trControl：定义函数运行参数的列表
  saveRDS(conservation,file = './conservation.model.rds')
  #评估模型（训练内评估：CV）
  conservation = readRDS('./conservation.model.rds')
  conf_matrix <- confusionMatrix(conservation$pred$pred, conservation$pred$obs) 
  Matt_Coef <- function (conf_matrix){
    TP <- conf_matrix$table[1,1]
    TN <- conf_matrix$table[2,2]
    FP <- conf_matrix$table[1,2]
    FN <- conf_matrix$table[2,1]
    
    mcc_num <- TP*TN - FP*FN
    mcc_den <- as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
    
    mcc_final <- mcc_num/sqrt(mcc_den)
    return(mcc_final)
  }
  mcc <- Matt_Coef(conf_matrix)
  #填入结果表
  final_result[1,i] <- mean(conservation$results$ROC)
  final_result[2,i] <- mean(conservation$results$Sens)
  final_result[3,i] <- mean(conservation$results$Spec)
  final_result[4,i] <- conf_matrix$overall[1]
  final_result[5,i] <- mcc
  #评估模型（测试集评估）
  pred <- predict(conservation, newdata = test_data, type = "prob")
  print(length(label_test) == nrow(pred))
  BIOmotifvsnon_testppred <- prediction(pred$positive,label_test)
  BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc")
  final_result[6,i] <- BIOmotifvsnon_testpppred_auc@y.values[[1]][1]
  
  result <- as.data.frame(matrix(data=NA,2,2))
  result[1,1] <- length(which(pred[1:(length(label_test)/2),1] > 0.5)) #前面一半数据是positive，得到预测positie的样本个数，TP
  result[1,2] <- length(which(pred[1:(length(label_test)/2),1] < 0.5)) #FN
  result[2,1] <- length(which(pred[((length(label_test)/2) + 1):length(label_test),1] > 0.5)) #FP
  result[2,2] <- length(which(pred[((length(label_test)/2) + 1):length(label_test),1] < 0.5)) #TN
  #填入结果表
  final_result[7,i] <- result[1,1]/( result[1,1] + result[1,2] )
  final_result[8,i] <- result[2,2]/( result[2,1] + result[2,2] )
  final_result[9,i] <- ( result[1,1] + result[2,2] )/ nrow(test_data)
  final_result[10,i] <- (result[1,1]*result[2,2]-result[1,2]*result[2,1])/
    (sqrt(as.numeric(result[1,1]+result[2,1])*(result[1,1]+result[1,2])*(result[2,2]+result[2,1])*(result[2,2]+result[1,2])))
  
  print(final_result)
  print(paste0("Round_",i,"_finished"))
}

saveRDS(final_result,'./final_result.rds')
readRDS('./final_result.rds')








