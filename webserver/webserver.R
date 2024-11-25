#制作网站要用的包：
library(jsonlite)

input_json <- commandArgs(trailingOnly = T)
jobID <- input_json[1]

a <- as.data.frame(fromJSON(paste0('/var/www/html/m6A-CAPred/job/',jobID,'/',jobID,'_para.json')))
file <- as.character(a$file)
jobID <- as.character(a$jobID)
cutoff <- as.double(a$cutoff)

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

################################################################################
###########################     设置工作目录      ##############################
################################################################################
target_dir1<-file #阿里云的网址
target_dir2<-'/home/jiaming/webserver_home' #本地服务器
setwd(target_dir2)
################################################################################
###############################      输入       ################################
################################################################################
#data <- read.table(target_dir1)
#data <- makeGRangesFromDataFrame(data, keep.extra.columns = T)
#data <- read.table('/home/jiaming/webserver_ali/dQR7BxnZCk.txt')

data <- read.table(target_dir1, header = TRUE)
#data <- read.table('/home/jiaming/webserver_ali/dQR7BxnZCk.txt', header = TRUE)
data <- GRanges(seqnames = data$seqnames,
                ranges = IRanges(start = data$start,
                                 end = data$end,
                                 width=1),
                strand = data$strand)
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

GF_data <- GFgenreation_m6A(data)
#saveRDS(GF_data,'./GF_data.rds')
#GF_data <- readRDS('./GF_data.rds')

#sequence feature
SeqFgeneration <- function(data,GF){ 
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

#用data和GF_data提取序列特征
data$reference_sequence <- DNAStringSet(Views(Hsapiens,data + 20))
sequenceFeature1 <- SeqFgeneration(as.character(data$reference_sequence),GF_data)
sequenceFeature2 <- SeqFgeneration2(as.character(data$reference_sequence),GF_data)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
dataFeatureBoth <- cbind(GF_data,sequenceFeature)

#更新data
data <- dataFeatureBoth

################################################################################
###############################     打分     ###################################
################################################################################
#评估模型（测试集评估）
conservation <- readRDS(paste0(target_dir2,'/conservation.model.rds')) 
pred <- predict(conservation, newdata = data, type = "prob")

################################################################################
############################     网页表格     ##################################
################################################################################
df <- read.table(target_dir1, header = TRUE)
#df <- read.table('/home/jiaming/webserver_ali/dQR7BxnZCk.txt', header = TRUE)
pred <- data.frame(jobID=jobID,
                   seqnames=as.character(df$seqnames),
                   position=as.character(df$start),
                   strand=as.character(df$strand),
                   probability=pred$positive,
                   cutoff=cutoff,
                   ifCancerRelated=pred$positive>cutoff)
################################################################################
###############################     存为json     ###############################
################################################################################
write_json(pred, paste0('/var/www/html/m6A-CAPred/job/',jobID,'/result/pred.json'))
write.csv(pred, paste0('/var/www/html/m6A-CAPred/job/',jobID,'/result/pred.csv'), row.names = FALSE)























