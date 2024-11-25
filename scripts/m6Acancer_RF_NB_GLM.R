################################################################################
##############################      Load R Packages      ########################
################################################################################
library(e1071) # Support Vector Machine (SVM) Model
library(ROCR) # Receiver Operating Characteristic (ROC) Curve
library(pROC) # Accuracy Metrics for Machine Learning
library(dplyr) # Feature Selection and Dataframe Operations
library(m6ALogisticModel) # Custom Machine Learning Package
library(SummarizedExperiment) # Genomic Data Handling
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # Genomic Data Handling
library(BSgenome.Hsapiens.UCSC.hg19) # Genomic Data Handling
library(MLmetrics) # Model Evaluation Metrics
library(caret) # Cross-Validation
# Packages for Feature Construction:
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(phastCons100way.UCSC.hg19)
library(fitCons.UCSC.hg19)
# For Boxplot Visualization
library(ggplot2)
################################################################################
###########################     Set Working Directory      #####################
################################################################################
setwd('/Users/jiaming/Desktop/cancer_m6A_predictor/file_m6Acancer_RF_NB_GLM')
################################################################################
############################      Data Analysis      ############################
################################################################################
# Load Raw Data
data <- readRDS('./m6A_hg19.rds')
# Extract Cancer-Related Columns
data <- data@elementMetadata[,-c(1:3)]
# Replace Invalid Characters
names(data) <- gsub('-','_',names(data))
# Count the Occurrences of Each m6A Site in Cancer Cell Lines
row_sums <- apply(data, 1, sum)
# Add Count to Data
data$row_sums <- row_sums
# Add Back Genomic Site Information
dataCopy <- readRDS('./m6A_hg19.rds')
data <- GRanges(seqnames=dataCopy@seqnames,ranges=dataCopy@ranges,strand=dataCopy@strand,data)
# Reorder Data by Count (Cancer Index)
data <- data[order(data$row_sums),]
# Result: 0 -> 33920 (Negative Samples), Select 3392 (Positive Samples)
# Remove Count Column
data <- data[,-26]
################################################################################
########################      Extract Positive and Negative Samples      #######
################################################################################
# Select First 33920 Rows as Negative Samples, Last 3392 Rows as Positive Samples
Nsample <- data[c(1:33920),] # 33920 Rows
Psample <- data[c(130692:134083),] # 3392 Rows
# Remove Unused Data
rm(data)
rm(dataCopy)
rm(row_sums)
################################################################################
############################      Feature Construction      ####################
################################################################################
# Genomic Features
# Define Function GFgeneration_m6A
GFgeneration_m6A <- function(data){
  
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

GF_Nsample <- GFgeneration_m6A(Nsample)
GF_Psample <- GFgeneration_m6A(Psample)

saveRDS(GF_Nsample,'GF_Nsample.rds')
saveRDS(GF_Psample,'GF_Psample.rds')

GF_Nsample <- readRDS('./GF_Nsample.rds')
GF_Psample <- readRDS('./GF_Psample.rds')

# Sequence Features
SeqFgeneration <- function(data,GF){ # GF is Redundant
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

# Extract Sequence Features Using Nsample and GF_Nsample
Nsample$reference_sequence <- DNAStringSet(Views(Hsapiens,Nsample + 20))
sequenceFeature1 <- SeqFgeneration(as.character(Nsample$reference_sequence),GF_Nsample)
sequenceFeature2 <- SeqFgeneration2(as.character(Nsample$reference_sequence),GF_Nsample)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
NFeatureBoth <- cbind(GF_Nsample,sequenceFeature)

# Extract Sequence Features Using Psample and GF_Psample
Psample$reference_sequence <- DNAStringSet(Views(Hsapiens,Psample + 20))
sequenceFeature1 <- SeqFgeneration(as.character(Psample$reference_sequence),GF_Psample)
sequenceFeature2 <- SeqFgeneration2(as.character(Psample$reference_sequence),GF_Psample)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
PFeatureBoth <- cbind(GF_Psample,sequenceFeature)

# Update Nsample and Psample
Nsample <- NFeatureBoth
Psample <- PFeatureBoth

Nsample <- Nsample[,-which(colnames(Nsample) %in% c('long_exon','WTAP_TREW','last_exon','last_exon_400bp','exon_stop','constitutive_exon','last_exon_sc400','Stop_codons','UTR3','HK_genes','miR_targeted_genes'))] 
Psample <- Psample[,-which(colnames(Psample) %in% c('long_exon','WTAP_TREW','last_exon','last_exon_400bp','exon_stop','constitutive_exon','last_exon_sc400','Stop_codons','UTR3','HK_genes','miR_targeted_genes'))]

################################################################################
#############################    Build Result Table    #########################
################################################################################
final_result <- as.data.frame(matrix(data = NA,10,10))
rownames(final_result) <- c("AUROC_5fcv","Sn_5fcv","Sp_5fcv","ACC_5fcv","MCC_5fcv","AUROC","Sn","Sp","ACC","MCC")
################################################################################
#############################     Train Model     ##############################
################################################################################
for(i in 1:10) {
  print(paste0('round_',i,'_starts'))
  # Use 80% of Positive Samples (2714) for Training and 20% (678) for Testing
  set.seed(i)
  train_P_indx <- sample(1:nrow(Psample),2714)
  train_P <- Psample[train_P_indx,]
  test_P_indx <- setdiff(1:nrow(Psample),train_P_indx)
  test_P <- Psample[test_P_indx,]
  # Use 2714 Negative Samples for Training and 678 for Testing
  set.seed(i)
  train_N_indx <- sample(1:nrow(Nsample),2714)
  train_N <- Nsample[train_N_indx,]
  test_N_indx_pool<- setdiff(1:nrow(Nsample),train_N_indx)
  test_N_indx <- sample(test_N_indx_pool,678) 
  test_N <- Nsample[test_N_indx,]
  # Combine train_P and train_N and Add Labels
  train_data <- rbind(train_P,train_N)
  label_train <- c(rep("positive",nrow(train_P)),rep("unmodified",nrow(train_N)))
  label_train <- factor(label_train,labels=c("positive","unmodified"))
  train_data$label <- label_train
  # Combine test_P and test_N and Add Labels
  test_data <- rbind(test_P, test_N)
  label_test <- c(rep(1, (nrow(test_P))), rep(0, (nrow(test_N))))
  test_data$label <- label_test
  # Set Hyperparameters
  fitControl <- trainControl(method = "cv",  # Validation Method
                             number = 5,    # Number of Folds
                             savePred = TRUE,  # Save Predictions During Training
                             summaryFunction = twoClassSummary, # Function to Summarize Cross-Validation Results
                             classProbs = TRUE) # Calculate Class Probabilities
  conservation <- train(label ~ ., data = train_data, 
                        method = "nb", # Random Forest: "rf", Naive Bayes: "nb", Generalized Linear Model: "glm"
                        preProc = c("center", "scale"), # Preprocessing Methods: Centering and Scaling
                        trControl = fitControl) # Training Parameters
  saveRDS(conservation, file = './conservation.model.nb.deleted.rds')
  # Evaluate Model (Within-Training Evaluation: CV)
  conservation = readRDS('./conservation.model.nb.deleted.rds')
  conf_matrix <- confusionMatrix(conservation$pred$pred, conservation$pred$obs)
  
  Matt_Coef <- function(conf_matrix){
    TP <- conf_matrix$table[1, 1]
    TN <- conf_matrix$table[2, 2]
    FP <- conf_matrix$table[1, 2]
    FN <- conf_matrix$table[2, 1]
    
    mcc_num <- TP * TN - FP * FN
    mcc_den <- as.double((TP + FP)) * as.double((TP + FN)) * as.double((TN + FP)) * as.double((TN + FN))
    
    mcc_final <- mcc_num / sqrt(mcc_den)
    return(mcc_final)
  }
  mcc <- Matt_Coef(conf_matrix)
  
  # Fill in Result Table
  final_result[1, i] <- mean(conservation$results[2,]$ROC)
  final_result[2, i] <- mean(conservation$results[2,]$Sens)
  final_result[3, i] <- mean(conservation$results[2,]$Spec)
  final_result[4, i] <- conf_matrix$overall[1]
  final_result[5, i] <- mcc
  
  # Evaluate Model (Test Set Evaluation)
  pred <- predict(conservation, newdata = test_data, type = "prob")
  print(length(label_test) == nrow(pred))
  BIOmotifvsnon_testppred <- prediction(pred$positive, label_test)
  BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred, "auc")
  final_result[6, i] <- BIOmotifvsnon_testpppred_auc@y.values[[1]][1]
  
  result <- as.data.frame(matrix(data = NA, 2, 2))
  result[1, 1] <- length(which(pred[1:(length(label_test)/2), 1] > 0.5)) # TP
  result[1, 2] <- length(which(pred[1:(length(label_test)/2), 1] < 0.5)) # FN
  result[2, 1] <- length(which(pred[((length(label_test)/2) + 1):length(label_test), 1] > 0.5)) # FP
  result[2, 2] <- length(which(pred[((length(label_test)/2) + 1):length(label_test), 1] < 0.5)) # TN
  
  # Fill in Result Table
  final_result[7, i] <- result[1, 1] / (result[1, 1] + result[1, 2]) # Sensitivity
  final_result[8, i] <- result[2, 2] / (result[2, 1] + result[2, 2]) # Specificity
  final_result[9, i] <- (result[1, 1] + result[2, 2]) / nrow(test_data) # Accuracy
  final_result[10, i] <- (result[1, 1] * result[2, 2] - result[1, 2] * result[2, 1]) / 
    (sqrt(as.numeric(result[1, 1] + result[2, 1]) * (result[1, 1] + result[1, 2]) * 
            (result[2, 2] + result[2, 1]) * (result[2, 2] + result[1, 2]))) # MCC
  
  print(final_result)
  print(paste0("Round_", i, "_finished"))
}

saveRDS(final_result, './final_result.nb.deleted.rds')
readRDS('./final_result.nb.deleted.rds')

# Available Methods for `train()` in caret:
# - lm: Linear Regression
# - rf: Random Forest
# - svmRadial: Radial Support Vector Machine
# - glm: Generalized Linear Model
# - knn: k-Nearest Neighbors
# - nb: Naive Bayes
# - rpart: Recursive Partitioning and Regression Trees
# - gbm: Gradient Boosting Machine
# - xgbLinear: Extreme Gradient Boosting
# - ctree: Conditional Inference Trees
# - cforest: Conditional Random Forest
# - lda: Linear Discriminant Analysis
# - qda: Quadratic Discriminant Analysis
# - pls: Partial Least Squares
# - pca: Principal Component Analysis
# - deepnet: Deep Neural Networks
# - rfns: Random Forest with Non-Parametric Survival

# Note: Some algorithms require additional R packages (e.g., xgboost requires the xgboost package).

  
