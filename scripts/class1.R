library(BSgenome)
library(matrixStats)

# Function NC1(): Split DNAStringSet (string) into a matrix, each row represents a DNA sequence
NC1 <- function(data) {
  for (i in 1:length(data)) {
    if (i == 1) {
      # Split the first string into characters
      # `strsplit()` splits strings by a specified delimiter (here it's ""), creating character vectors
      # `fixed = TRUE` ensures plain text matching instead of regex
      DCfirst <- unlist(strsplit(data[1], "", fixed = TRUE))
      # Create an empty matrix with rows equal to the number of sequences in `data`
      # and columns equal to the length of the first sequence
      DCsecond <- matrix(NA, nrow = length(data), ncol = length(DCfirst))
      # Assign the first sequence's split characters to the first row
      DCsecond[1, ] <- DCfirst
    } else {
      # For subsequent rows, split the string and assign to the matrix
      DCsecond[i, ] <- unlist(strsplit(data[i], "", fixed = TRUE))
    }
  }
  # Return the matrix where each row corresponds to a DNA sequence split into individual bases
  return(DCsecond)
}


# Function CONPOSITION:
# Calculates nucleotide composition or frequency for sequences.
CONPOSITION <- function(Data, NI = 3, NTYPE = "RNA", Freq = 2) {
  # Set U/T based on RNA or DNA type
  if (NTYPE == "RNA") {
    U <- "U"
  } else {
    U <- "T"
  }
  
  # Generate all possible nucleotide combinations up to `NI` length
  for (i in 1:NI) {
    if (i == 1) {
      NP <- c("A", "G", "C", U)
    } else {
      NP1 <- NULL
      for (j in c("A", "G", "C", U)) {
        for (k in 1:length(NP)) {
          NP1 <- c(NP1, paste0(NP[k], j))
        }
      }
      NP <- NP1
    }
  }
  
  # Initialize result matrix with nucleotide combinations as columns
  MA2 <- matrix(NA, ncol = length(NP), nrow = length(Data))
  colnames(MA2) <- NP
  
  # Count nucleotide combinations in sequences
  if (NTYPE == "RNA") {
    for (i in 1:length(NP)) {
      MA2[, i] <- vcountPattern(NP[i], RNAStringSet(Data))
    }
  } else {
    for (i in 1:length(NP)) {
      MA2[, i] <- vcountPattern(NP[i], DNAStringSet(Data))
    }
  }
  
  # Convert counts to frequencies if Freq == 1
  M3 <- MA2
  if (Freq == 1) {
    for (i in 1:nrow(M3)) {
      M3[i, ] <- MA2[i, ] / sum(MA2[i, ])
    }
  }
  return(M3)
}

# Function sequenceFeatures:
# Calculates cumulative frequency of each nucleotide at each position in sequences.
sequenceFeatures <- function(Data, NTYPE = "RNA") {
  # Convert sequences into a character matrix
  sequences_M <- NC1(Data)
  
  # Set U/T based on RNA or DNA type
  if (NTYPE == "RNA") {
    U <- "U"
  } else {
    U <- "T"
  }
  
  N <- ncol(sequences_M) # Number of columns (sequence length)
  
  # Calculate cumulative counts for each nucleotide
  cumFreq_A <- rowCumsums(matrix(as.numeric(sequences_M == "A"), ncol = N, byrow = FALSE))
  cumFreq_T <- rowCumsums(matrix(as.numeric(sequences_M == U), ncol = N, byrow = FALSE))
  cumFreq_C <- rowCumsums(matrix(as.numeric(sequences_M == "C"), ncol = N, byrow = FALSE))
  cumFreq_G <- rowCumsums(matrix(as.numeric(sequences_M == "G"), ncol = N, byrow = FALSE))
  
  # Combine cumulative counts into a single matrix
  cumFreq_combined <- matrix(0, ncol = N, nrow = length(Data))
  cumFreq_combined[sequences_M == "A"] <- cumFreq_A[sequences_M == "A"]
  cumFreq_combined[sequences_M == U] <- cumFreq_T[sequences_M == U]
  cumFreq_combined[sequences_M == "C"] <- cumFreq_C[sequences_M == "C"]
  cumFreq_combined[sequences_M == "G"] <- cumFreq_G[sequences_M == "G"]
  
  # Normalize cumulative frequencies by position
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined) <- paste0("cumFreq_", seq_len(N))
  
  return(cumFreq_combined)
}

# Example Usage:
# TEST <- c("ACGUCUCUAUCGUACGUACGUAGUG", "CUUCUCGUAACGUAGCAUACGGAUG")
# sequenceFeatures(TEST, "RNA")
# CONPOSITION(TEST)
