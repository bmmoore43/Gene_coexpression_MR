#!/usr/bin/env Rscript

## 20190921
## Jen Wisecaver
## Normalize and transform an RNAseq counts matrix
##
## See transform_count jupyter notebook for walkthrough

library('getopt')
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help'        , 'h', 0, "logical", '',
  'abundances'  , 'a', 1, "character", 'Enter list of abundance estimates from Kallisto (see jupyter notebook tutorial)',
  'lookup'      , 'l', 1, "character", 'Enter transcript to gene lookup file (see jupyter notebook tutorial)',
  'samples'     , 's', 1, "character", 'Enter experimental design file (see jupyter notebook tutorial)',
  'seqtype'     , 't', 1, "character", 'Full-length or 3 prime tagged RNAseq? Enter [Full/full/F/f or Tag/tag/T/t]',
  'bltype'      , 'b', 1, "character", 'Blind vst transformation to exp design? Enter [FALSE/False/false/F/f or TRUE/True/true/T/t]',
  'out'         , 'o', 1, "character", 'Enter base name for all output files'
), byrow=TRUE, ncol=5)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( is.null(opt$abundances) || is.null(opt$lookup) || is.null(opt$samples) || is.null(opt$out)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

filesfile <- opt$abundances
#cat(system((paste('head -n 5', filesfile)), intern=TRUE), sep='\n')

library(stringr)
files <- scan(filesfile, what = 'list')
file_names <- str_split(files, '/', simplify = TRUE)
file_names <- file_names[,ncol(file_names)-1]
names(files) <- file_names
#head(files)

t2gfile <- opt$lookup
#cat(system((paste('head -n 5', t2gfile)), intern=TRUE), sep='\n')

tx2gene <- read.table(t2gfile, sep="\t", header=TRUE)
#head(tx2gene)

samplefile <- opt$samples
#cat(system((paste('head -n 5', samplefile)), intern=TRUE), sep='\n')

sampleTable <- read.table(samplefile, sep="\t", header=TRUE)
#head(sampleTable)

CFA <- 'scaledTPM'
if ( !is.null(opt$seqtype) ) {
	seqtype <- opt$seqtype

	if (seqtype == 'full' || seqtype == 'Full' || seqtype == 'F' || seqtype == 'f'){
		CFA <- 'scaledTPM'
	} 

	if (seqtype == 'tag' || seqtype == 'Tag' || seqtype == 'T' || seqtype == 't'){
		CFA <- 'no'
	} 
}

library(tximport)
library(rhdf5)
paste('tximport options: type = kallisto, tx2gene = tx2gene, countsFromAbundance = ', CFA)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = CFA)
#names(txi)


library("DESeq2")
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- estimateSizeFactors(dds)
#head(counts(dds))

#nrow(dds)
dds_filt <- dds[ rowSums(counts(dds)) > 0, ]
#nrow(dds_filt)

normalized_matrix <- counts(dds_filt, normalized = TRUE)
#head(normalized_matrix)

BLINDED <- 'FALSE'
if ( !is.null(opt$bltype) ) {
	bltype <- opt$bltype

	if (bltype == 'true' || bltype == 'True' || bltype == 'TRUE' || bltype == 'T' || bltype == 't'){
		BLINDED <- 'TRUE'
	} 

	if (bltype == 'false' || bltype == 'False' || bltype == 'FALSE' || bltype == 'F' || bltype == 'f'){
		BLINDED <- 'FALSE'
	} 
}

paste('vst options: blind = ', BLINDED)
vsd <- vst(dds_filt, blind = BLINDED)
normalized_vst_matrix <- assay(vsd)
#head(normalized_vst_matrix)

# Tally number of genes with expression greater than zero
df_binary <- normalized_matrix
df_binary[df_binary == 0] <- 0
df_binary[df_binary > 0] <- 1
stats_df <- as.data.frame(rowSums(df_binary))
names(stats_df) <- c('no. samples with count > 0')

# Tally number of genes with expression greater than or equal to five
df_binary <- normalized_matrix
df_binary[df_binary < 5] <- 0
df_binary[df_binary > 0] <- 1
stats_df$'no. samples with count >= 5' = rowSums(df_binary)

# Tally number of genes with expression greater than or equal to ten
df_binary <- normalized_matrix
df_binary[df_binary < 10] <- 0
df_binary[df_binary > 0] <- 1
stats_df$'no. samples with count >= 10' = rowSums(df_binary)

# Calculate median expression genes in the normalized matrix
median_expression  <- apply(normalized_matrix, 1, median)
stats_df$'pre-VST median expression' = median_expression

# Calculate mean expression genes in the normalized matrix
mean_expression  <- apply(normalized_matrix, 1, mean)
stats_df$'pre-VST mean expression' = mean_expression

# Calculate expression variance genes in the normalized matrix
exp_variance <- apply(normalized_matrix, 1, var)
stats_df$'pre-VST expression variance' = exp_variance

# Calculate median expression genes in the VST matrix
median_expression  <- apply(normalized_vst_matrix, 1, median)
stats_df$'post-VST median expression' = median_expression

# Calculate mean expression genes in the VST matrix
mean_expression  <- apply(normalized_vst_matrix, 1, mean)
stats_df$'post-VST mean expression' = mean_expression

# Calculate expression variance genes in the VST matrix
exp_variance <- apply(normalized_vst_matrix, 1, var)
stats_df$'post-VST expression variance' = exp_variance
#head(stats_df)

norm_outfil = paste(opt$out, '_normalized.matrix', sep = '')
write.table(normalized_matrix, norm_outfil, sep="\t", quote = FALSE)

normVst_outfil = paste(opt$out, '_normalized_vst_transformed.matrix', sep = '')
write.table(normalized_vst_matrix, normVst_outfil, sep="\t", quote = FALSE)

stats_outfil = paste(opt$out, '_statistics.txt', sep = '')
write.table(stats_df, stats_outfil, sep="\t", quote = FALSE)
