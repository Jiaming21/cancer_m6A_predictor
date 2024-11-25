################################################################################
##############################      Load R Packages      ########################
################################################################################
library(e1071) # Support Vector Machine (SVM) Model
library(ROCR) # Receiver Operating Characteristic (ROC) Curve
library(pROC) # Package for Calculating Machine Learning Accuracy
library(dplyr) # Feature Selection and Dataframe Operations
library(m6ALogisticModel) # Custom Machine Learning Package
library(SummarizedExperiment) # Genomic Data Processing
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # Genomic Data Processing
library(BSgenome.Hsapiens.UCSC.hg19) # Genomic Data Processing
library(MLmetrics) # Model Evaluation
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
setwd('/Users/jiaming/Desktop/cancer_m6A_predictor/file_m6Acancer_data')
################################################################################
############################      Data Analysis      ############################
################################################################################
# Read Raw Data
data <- readRDS('./m6A_hg19.rds')
# Extract Cancer-Related Columns
data <- data@elementMetadata[,-c(1:3)]
# Replace Invalid Characters
names(data) <- gsub('-', '_', names(data))
# Count the Occurrences of Each m6A Site in Cancer Cell Lines
row_sums <- apply(data, 1, sum)
# Add Count to Data
data$row_sums <- row_sums
# Add Back Genomic Site Information
dataCopy <- readRDS('./m6A_hg19.rds')
data <- GRanges(seqnames = dataCopy@seqnames, ranges = dataCopy@ranges, strand = dataCopy@strand, data)
# Reorder Data by Count (Cancer Index)
data <- data[order(data$row_sums),]
# Result: 0 -> 33920 (Negative Samples), Select 3392 (Positive Samples)
# Remove Count Column
data <- data[,-26]
saveRDS(data, file = './data.rds')
rm(dataCopy)
rm(row_sums)
################################################################################
############################      Feature Construction      ####################
################################################################################
# Genomic Features
# Define GFgenreation_m6A Function
GFgenreation_m6A <- function(data) {
  
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
  
  data_standardized <- predictors_annot(
    se = matureSE,
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
    standardization = TRUE
  )
  
  GF <- mcols(data_standardized)
  return(GF)
}

GF_data <- GFgenreation_m6A(data)

saveRDS(GF_data, 'GF_data.rds')

GF_data <- readRDS('GF_data.rds')

# Sequence Features
SeqFgeneration <- function(data, GF) { 
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
SeqFgeneration2 <- function(data, GF) {
  UntestVarSeq <- data
  source("../class1.R")
  source("../class2.R")
  source("../class3.R")
  source("../class4.R")
  source("../class5.R")
  source("../class6.R")
  NF <- sequenceFeatures(UntestVarSeq, NTYPE = "DNA")
  return(NF)
}

# Extract Sequence Features Using data and GF_data
data$reference_sequence <- DNAStringSet(Views(Hsapiens, data + 20))
sequenceFeature1 <- SeqFgeneration(as.character(data$reference_sequence), GF_data)
sequenceFeature2 <- SeqFgeneration2(as.character(data$reference_sequence), GF_data)
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)], sequenceFeature2[,-1])
dataFeatureBoth <- cbind(GF_data, sequenceFeature)

# Update data
data <- dataFeatureBoth

################################################################################
###############################     Scoring     #################################
################################################################################
# Evaluate Model (Test Set Evaluation)
conservation <- readRDS('conservation.model.rds')
pred <- predict(conservation, newdata = data, type = "prob")

# Save and Load Predictions
saveRDS(pred, file = 'pred.rds')
final_predictions <- readRDS('pred.rds')

# Output Prediction Results
print(final_predictions)

# End of Script










