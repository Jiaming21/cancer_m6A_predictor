# SNP density analysis
library(GenomicRanges)
################################################################################
###########################     Set Working Directory      #####################
################################################################################
setwd('/Users/jiaming/Desktop/cancer_m6A_predictor/file_SNP_density')
################################################################################
#############################     Load Files      ##############################
################################################################################
somaticSNP <- readRDS('./somaticSNP_hg19_all.rds') # Load somatic SNP data
pred <- readRDS('./pred.rds') # Load prediction results
data <- readRDS('./data.rds') # Load dataset
################################################################################
#######################     Split Data by Cutoff     ###########################
################################################################################
# Define cutoff values: 0.3 - 0.9
cancerUnrelated <- data[pred$positive < 0.9,] # Data points with low cancer-related scores
length(cancerUnrelated) # Number of unrelated points

cancerRelated <- data[pred$positive >= 0.9,] # Data points with high cancer-related scores
length(cancerRelated) # Number of related points
################################################################################
##########################     Find SNP Positions      #########################
################################################################################
olp1 <- findOverlaps(somaticSNP, cancerUnrelated + 2, ignore.strand = TRUE) # Overlap with unrelated group
qhits1 <- unique(queryHits(olp1)) # Unique SNP positions
length(qhits1) # Count of SNP positions in unrelated group

olp2 <- findOverlaps(somaticSNP, cancerRelated + 2, ignore.strand = TRUE) # Overlap with related group
qhits2 <- unique(queryHits(olp2)) # Unique SNP positions
length(qhits2) # Count of SNP positions in related group
################################################################################
############################     SNP Density      ##############################
################################################################################
length(qhits1) / length(cancerUnrelated) # SNP density in unrelated group

length(qhits2) / length(cancerRelated) # SNP density in related group
################################################################################
###################      Statistical Significance Test     #####################
################################################################################
# Explanation:
# - Left column: Cancer-related status (Yes/No)
# - Top row: Single nucleotide polymorphism (SNP) counts or their complements
# Cancer-unrelated SNP count: length(qhits1)
# Non-SNP count in unrelated region: length(cancerUnrelated) - length(qhits1)
# Cancer-related SNP count: length(qhits2)
# Non-SNP count in related region: length(cancerRelated) - length(qhits2)

# Create a contingency table for chi-squared test
table1 <- matrix(
  c(
    length(qhits1),  # SNP count in cancer-unrelated group
    length(qhits2),  # SNP count in cancer-related group
    length(cancerUnrelated) - length(qhits1), # Non-SNP count in unrelated group
    length(cancerRelated) - length(qhits2)    # Non-SNP count in related group
  ),
  nrow = 2
)

# Perform chi-squared test
chisq_result <- chisq.test(table1)
print(chisq_result)

# Note on `prop.test` usage:
# The prop.test function uses SNP counts (first column) and total counts (third column).
# Syntax: prop.test(x, n, ...)
# x: Vector of SNP counts
# n: Vector of total counts for each group
x <- c(length(qhits1), length(qhits2)) # SNP counts
n <- c(length(cancerUnrelated), length(cancerRelated)) # Total data points in each group

# Perform proportion test
prop_result <- prop.test(x, n)
print(prop_result)

################################################################################
###################      Incorrect Usage for Proportion Test      ##############
################################################################################
# Incorrect example: Passing a matrix directly to `prop.test`
# Explanation: `prop.test` expects two separate vectors (counts and totals),
# not a matrix.

table2 <- matrix(
  c(
    length(qhits1), 
    length(qhits2), 
    length(cancerUnrelated), 
    length(cancerRelated)
  ),
  nrow = 2
)

# Uncommenting the following line will throw an error:
# prop.test(table2)

# Explanation of the error:
# `prop.test` requires:
# - x: A vector of SNP counts
# - n: A vector of total counts
# Passing a matrix as a single argument does not match the expected input format.

################################################################################
###################      Significance Level and P-value      ###################
################################################################################
# Interpreting the p-value:
# A p-value < 0.05 is considered statistically significant.
# This indicates a meaningful difference in SNP proportions between cancer-related 
# and unrelated regions.
