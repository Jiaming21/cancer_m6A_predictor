# Predict probabilities for Usample
pred <- predict(conservation, newdata = Usample, type = "prob")

# Load the raw dataset
data <- readRDS('./m6A_hg19.rds')

# Extract cancer-related columns
data <- data@elementMetadata[,-c(1:3)]

# Replace invalid characters in column names
names(data) <- gsub('-', '_', names(data))

# Count occurrences of each m6A site in cancer cell lines
row_sums <- apply(data, 1, sum)

# Add counts to the dataset
data$row_sums <- row_sums

# Add back genomic site information
dataCopy <- readRDS('./m6A_hg19.rds')
data <- GRanges(seqnames = dataCopy@seqnames, ranges = dataCopy@ranges, strand = dataCopy@strand, data)

# Reorder data by cancer index (row_sums)
data <- data[order(data$row_sums),]

# Extract row_sums corresponding to pred
count <- data$row_sums[c(33921:130691)]

# Extract positive predictions
positive <- pred$positive

# Convert count and positive to dataframes
count <- as.data.frame(count)
positive <- as.data.frame(positive)

# Combine count and positive into a single dataframe
sample <- cbind(count, positive)

# Save the combined dataframe
saveRDS(sample, './sample.rds')

# Read the saved sample
sample <- readRDS('./sample.rds') # First column: m6A site counts in cancer cell lines, Second column: positive prediction probabilities

# Convert x-axis variable to factor for plotting
sample$count <- as.factor(sample$count)

################################################################################
###################      ggplot+geom_boxplot() for Boxplots      ################
################################################################################

# Create a boxplot of positive probabilities grouped by counts
p1 <- ggplot(sample, aes(x = count, y = positive)) + geom_boxplot()

# Display the plot
p1

# Save the plot as a PDF
ggsave('./boxplot1.pdf')

################################################################################
###########################      Grouped Boxplots      ##########################
################################################################################

# Convert count to numeric for grouping
sample$count <- as.numeric(sample$count)

# Add a new column for grouped counts
sample$count1 <- '[0-5]'
sample[sample$count > 5, 'count1'] <- '[6-10]'
sample[sample$count > 10, 'count1'] <- '[11-15]'
sample[sample$count > 15, 'count1'] <- '[16-19]'

# Fix the order of levels for the x-axis
sample$count1 <- factor(sample$count1, levels = c('[0-5]', '[6-10]', '[11-15]', '[16-19]'))

# Create a grouped boxplot
p2 <- ggplot(sample, aes(x = count1, y = positive)) + geom_boxplot()

# Display the grouped boxplot
p2

# Save the grouped boxplot as a PDF
ggsave('./boxplot2.pdf')





