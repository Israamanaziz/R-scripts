# run Differential Gene Expression using DESeq2 in R, most basic workflow

#Assuming all the required packages are installed beforehand

#load Libraries

library(DESeq2)
library(tidyverse)
library(airway)

#Data Wrangling: Retrieve Data from airway package

data(airway)
airway

#Convert data into simple Rows as Samples and columns as info

sample_info <- as.data.frame(colData(airway))

# Grab column 2 and 3 only

sample_info <- sample_info[,c(2,3)]

#change labels into more readable format

sample_info$dex <- gsub('trt','treated', sample_info$dex)
sample_info$dex <- gsub('untrt','untreated', sample_info$dex)

# Change column names into more readable format

colnames(sample_info) <- c('CellLines','Dexamethasone')

#Save files in csv format, Make sure your files are getting saved in the required directory

write.table(sample_info,file = 'sample_info.csv',sep = ',', col.names = T, row.names = T, quote = F)

# Extract actual count matrix from airway

countsData <- assay(airway)


# Save files in csv format

write.table(countsData, file = 'Counts_data.csv', sep = ',', col.names = T,row.names = T,quote = F)


# RUN DESeq2 Workflow
# Read Counts Data

counts_data <- read.csv('Counts_data.csv')
head(counts_data)


#read Sample info

columnData <- read.csv('sample_info.csv')


#Now we have samples as column in counts_data, and samples as rows in columnData
#Make sure they match each other

all(colnames(counts_data) %in% rownames(columnData))

#Make sure they are same in order

all(colnames(counts_data) == rownames(columnData))


#Construct the DESeqDataset(dds) object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = columnData,
                              design = ~ Dexamethasone)

dds


#prefiltering low Quality reads

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds


#Set factor level

dds$Dexamethasone <- relevel(dds$Dexamethasone, ref = "untreated")

#RUN DESeq2 

dds <- DESeq(dds)
res <- results(dds)


#Explore results

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# Construct MA plot

plotMA(res)

