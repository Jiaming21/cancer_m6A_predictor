##############################################################################################################
#### SNP density
##############################################################################################################
library(GenomicRanges)
somaticSNP <- readRDS("~/conservation/SNPs/somaticSNP_annotated.rds")
germlineSNP <- readRDS("~/conservation/SNPs/germlineSNP_annotated.rds")

hg19InMm10_m6AConservationInfo <- readRDS("~/conservation/hg19InMm10_m6AConservationInfo.rds")

# take care non-liftover sites
# hg19InMm10_m6AConservationInfo[which(hg19InMm10_m6AConservationInfo$liftoverToMm10 == 'NO')]$PUlearningHg19_score <- 0

hg19InMm10_m6AConservationInfo$score <- hg19InMm10_m6AConservationInfo$br_score + hg19InMm10_m6AConservationInfo$brts_score + hg19InMm10_m6AConservationInfo$nodsHg19_score + hg19InMm10_m6AConservationInfo$nodsmm10_score + hg19InMm10_m6AConservationInfo$sequence_score + hg19InMm10_m6AConservationInfo$PUlearningHg19_score + hg19InMm10_m6AConservationInfo$PUlearningMm10_score + hg19InMm10_m6AConservationInfo$DP2
#na <- na.omit(hg19InMm10_m6AConservationInfo$DP2)
#na2 <- as.numeric(attr(na,"na.action"))
#hg19InMm10_m6AConservationInfo <- hg19InMm10_m6AConservationInfo[-na2]

hg19InMm10_m6AConservationInfo$score <- hg19InMm10_m6AConservationInfo$score/8
hg19InMm10_m6AConservationInfo <- hg19InMm10_m6AConservationInfo[order(-hg19InMm10_m6AConservationInfo$score)]

high <- hg19InMm10_m6AConservationInfo[1:(round(length(hg19InMm10_m6AConservationInfo)*0.3))]
medium <- hg19InMm10_m6AConservationInfo[((round(length(hg19InMm10_m6AConservationInfo)*0.3)) + 1):(round(length(hg19InMm10_m6AConservationInfo)*0.6))]
medium2 <- hg19InMm10_m6AConservationInfo[((round(length(hg19InMm10_m6AConservationInfo)*0.6)) + 1):(round(length(hg19InMm10_m6AConservationInfo)*1))]
#low <- hg19InMm10_m6AConservationInfo[((round(length(hg19InMm10_m6AConservationInfo)*0.75)) + 1):length(hg19InMm10_m6AConservationInfo)]

# high <- hg19InMm10_m6AConservationInfo[which(hg19InMm10_m6AConservationInfo$score >= 5)]
# medium <- hg19InMm10_m6AConservationInfo[which(hg19InMm10_m6AConservationInfo$score >= 2 & hg19InMm10_m6AConservationInfo$score < 5)]
# medium2 <- hg19InMm10_m6AConservationInfo[which(hg19InMm10_m6AConservationInfo$score < 2)]


#################################### somatic
olp1 <- findOverlaps(somaticSNP,low + 2, ignore.strand = T)
qhits1 <- unique(queryHits(olp1))
length(qhits1)/length(low) 

olp2 <- findOverlaps(somaticSNP,medium2 + 2, ignore.strand = T)
qhits2 <- unique(queryHits(olp2))
length(qhits2)
length(qhits2)/length(medium2) # 5.24%

olp3 <- findOverlaps(somaticSNP,medium + 2, ignore.strand = T)
qhits3 <- unique(queryHits(olp3))
length(qhits3)
length(qhits3)/length(medium) # 13.36%

olp4 <- findOverlaps(somaticSNP,high + 2, ignore.strand = T)
qhits4 <- unique(queryHits(olp4))
length(qhits4)
length(qhits4)/length(high) # 16.83%

# significant test
test_number  <- c( length(qhits3) , length(qhits4) )
total_number <- c( length(medium), length(high) )
high_test <- prop.test(test_number, total_number)
high_test

## deleterious level
length(which(somaticSNP[qhits1]$Deleterious_Level > 3))/length(qhits1) 

length(which(somaticSNP[qhits2]$Deleterious_Level > 3))
length(which(somaticSNP[qhits2]$Deleterious_Level > 3))/length(qhits2) # 4.43%

length(which(somaticSNP[qhits3]$Deleterious_Level > 3))
length(which(somaticSNP[qhits3]$Deleterious_Level > 3))/length(qhits3) # 14.59%

length(which(somaticSNP[qhits4]$Deleterious_Level > 3))
length(which(somaticSNP[qhits4]$Deleterious_Level > 3))/length(qhits4) # 18.43%

## disease analysis - ClinVar
length(which(somaticSNP$Clinvar_Number[qhits1] != 0))/length(qhits1) 

length(which(somaticSNP$Clinvar_Number[qhits2] != 0 | somaticSNP$GWAS_number[qhits2] != 0))
length(which(somaticSNP$Clinvar_Number[qhits2] != 0 | somaticSNP$GWAS_number[qhits2] != 0))/length(qhits2) # 0.29%

length(which(somaticSNP$Clinvar_Number[qhits3] != 0 | somaticSNP$GWAS_number[qhits3] != 0))
length(which(somaticSNP$Clinvar_Number[qhits3] != 0 | somaticSNP$GWAS_number[qhits3] != 0))/length(qhits3) # 0.31%

length(which(somaticSNP$Clinvar_Number[qhits4] != 0 | somaticSNP$GWAS_number[qhits4] != 0))
length(which(somaticSNP$Clinvar_Number[qhits4] != 0 | somaticSNP$GWAS_number[qhits4] != 0))/length(qhits4) # 0.47%

## mutType
length(which(somaticSNP[qhits1]$MutType == 'nonsynonymous SNV'))/length(qhits1) 

length(which(somaticSNP[qhits2]$MutType == 'nonsynonymous SNV'))
length(which(somaticSNP[qhits2]$MutType == 'nonsynonymous SNV'))/length(qhits2) # 44.75%

length(which(somaticSNP[qhits3]$MutType == 'nonsynonymous SNV'))
length(which(somaticSNP[qhits3]$MutType == 'nonsynonymous SNV'))/length(qhits3) # 59.69%

length(which(somaticSNP[qhits4]$MutType == 'nonsynonymous SNV'))
length(which(somaticSNP[qhits4]$MutType == 'nonsynonymous SNV'))/length(qhits4) # 61.46%

######################## germline
olp1 <- findOverlaps(germlineSNP,low + 2, ignore.strand = T)
qhits1 <- unique(queryHits(olp1))
length(qhits1)/length(low)

olp2 <- findOverlaps(germlineSNP,medium2 + 2, ignore.strand = T)
qhits2 <- unique(queryHits(olp2))
length(qhits2)
length(qhits2)/length(medium2) # 3.15%

olp3 <- findOverlaps(germlineSNP,medium + 2, ignore.strand = T)
qhits3 <- unique(queryHits(olp3))
length(qhits3)
length(qhits3)/length(medium) # 2.94%

olp4 <- findOverlaps(germlineSNP,high + 2, ignore.strand = T)
qhits4 <- unique(queryHits(olp4))
length(qhits4)
length(qhits4)/length(high) # 2.57%


## deleterious level
length(which(germlineSNP[qhits1]$Deleterious_Level > 3))/length(qhits1) 

length(which(germlineSNP[qhits2]$Deleterious_Level > 3))
length(which(germlineSNP[qhits2]$Deleterious_Level > 3))/length(qhits2) # 0.31%

length(which(germlineSNP[qhits3]$Deleterious_Level > 3))
length(which(germlineSNP[qhits3]$Deleterious_Level > 3))/length(qhits3) # 2.04%

length(which(germlineSNP[qhits4]$Deleterious_Level > 3))
length(which(germlineSNP[qhits4]$Deleterious_Level > 3))/length(qhits4) # 3.28%

## disease analysis - ClinVar
length(which(germlineSNP$Clinvar_Number[qhits1] != 0))/length(qhits1) 

length(which(germlineSNP$Clinvar_Number[qhits2] != 0 | germlineSNP$GWAS_number[qhits2] != 0))
length(which(germlineSNP$Clinvar_Number[qhits2] != 0 | germlineSNP$GWAS_number[qhits2] != 0))/length(qhits2) 

length(which(germlineSNP$Clinvar_Number[qhits3] != 0 | germlineSNP$GWAS_number[qhits3] != 0))
length(which(germlineSNP$Clinvar_Number[qhits3] != 0 | germlineSNP$GWAS_number[qhits3] != 0))/length(qhits3) 

length(which(germlineSNP$Clinvar_Number[qhits4] != 0 | germlineSNP$GWAS_number[qhits4] != 0))
length(which(germlineSNP$Clinvar_Number[qhits4] != 0 | germlineSNP$GWAS_number[qhits4] != 0))/length(qhits4) 

## mutType
length(which(germlineSNP[qhits1]$MutType == 'nonsynonymous SNV'))/length(qhits1) 

length(which(germlineSNP[qhits2]$MutType == 'nonsynonymous SNV'))
length(which(germlineSNP[qhits2]$MutType == 'nonsynonymous SNV'))/length(qhits2) # 13.71%

length(which(germlineSNP[qhits3]$MutType == 'nonsynonymous SNV'))
length(which(germlineSNP[qhits3]$MutType == 'nonsynonymous SNV'))/length(qhits3) # 27.78%

length(which(germlineSNP[qhits4]$MutType == 'nonsynonymous SNV'))
length(which(germlineSNP[qhits4]$MutType == 'nonsynonymous SNV'))/length(qhits4) # 30.52%

##############################################################################################################
#### statistics significant test
##############################################################################################################
total_number <- c( 198, 190 )

high_number  <- c( 73, 59 )

prop.test(high_number, total_number)