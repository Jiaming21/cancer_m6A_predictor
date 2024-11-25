################################################################################
##############################      Load R Packages      ########################
################################################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(e1071)
library(ROCR)
library(pROC)
library(dplyr)
library(caret)
library(m6ALogisticModel)
library(MLmetrics)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
################################################################################
###########################     Set Working Directory      #####################
################################################################################
setwd('/Users/huangjiaming/Desktop/cancer_m6A_predictor/file_featureSelection')
################################################################################
###########################     Remove Rows with NA      ########################
################################################################################
NsampleALL <- readRDS("./NFeatureBoth.rds") # Genomic and sequence features for negative data

# Identify and remove rows with NA
a <- NA
for (i in 1:length(names(NsampleALL))) {
  b <- which(is.na(NsampleALL[, i]))
  a <- c(a, b)
}
a <- unique(a)
a <- a[-1]

NsampleALL <- NsampleALL[-a, ] # Delete rows containing NA
################################################################################
##########################      Calculate Feature Importance      ##############
################################################################################
# Define Function: F_score_fn
F_score_fn <- function(dataframe) {
  # Convert logical features to numeric
  for (i in 1:(ncol(dataframe))) {
    if (is.logical(dataframe[, i]) == TRUE) {
      dataframe[, i] <- as.numeric(dataframe[, i])
    }
  }
  dataframe <- as.data.frame(dataframe)
  # Variables
  half_index <- nrow(dataframe) / 2
  pos_part <- dataframe[c(1:half_index), ] # Positive samples
  neg_part <- dataframe[-c(1:half_index), ] # Negative samples
  pos_mean <- colMeans(pos_part) # Mean of positive samples
  neg_mean <- colMeans(neg_part) # Mean of negative samples
  total_mean <- colMeans(dataframe) # Overall mean
  df_pos_mean <- matrix(NA, ncol = ncol(pos_part), nrow = nrow(pos_part))
  df_neg_mean <- matrix(NA, ncol = ncol(neg_part), nrow = nrow(neg_part))
  for (i in 1:nrow(pos_part)) {
    df_pos_mean[i, ] <- pos_mean
    df_neg_mean[i, ] <- neg_mean
  }
  # Numerator
  numerator <- (pos_mean - total_mean)^2 + (neg_mean - total_mean)^2
  # Denominator
  denominator <- 1 / (half_index - 1) * (colSums((pos_part - df_pos_mean)^2) + colSums((neg_part - df_neg_mean)^2))
  # F_score
  F_score <- as.data.frame(numerator / denominator)
  return(F_score)
}

# Load training data and calculate feature importance
train_data <- readRDS("./train_data.rds")
train_data <- train_data[, -213] # Remove label column
importance <- F_score_fn(train_data)
indx <- order(-importance$`numerator/denominator`)
indx <- rownames(importance)[indx]

saveRDS(importance, 'importance.rds')
saveRDS(indx, 'featureImportanceOrder.rds')
################################################################################
################################      Plot      ################################
################################################################################
importance <- readRDS("./importance.rds")
indx <- order(-importance$`numerator/denominator`)
rownames <- rownames(importance)
rownames <- as.data.frame(rownames)
rownames_arranged <- rownames[indx, ]

imp <- importance[indx, ]
imp <- as.data.frame(imp)
rownames(imp) <- rownames_arranged

imp$name <- factor(rownames(imp), levels = rownames(imp))
colnames(imp) <- c("importance", "name")
imp$importance2 <- scale(imp$importance, center = TRUE, scale = TRUE) # Normalization

# Plot top 24 features
a <- ggplot(imp[1:24, ], aes(name, importance)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90))

# Save the plot
ggsave('./FeatureImportance.pdf', a, height = 4, width = 5)


