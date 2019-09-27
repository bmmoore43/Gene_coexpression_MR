## 20190921
## Jen Wisecaver
## Normalize and transform an RNAseq counts matrix
##
## See transform_count jupyter notebook for walkthrough
## Usage: Rscript transform_counts.R counts.matrix conditions.txt

library('getopt')
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help'        , 'h', 0, "logical", '',
  'matrix'      , 'm', 1, "character", 'Enter path to preformated expression matrix for DEseq',
  'samples'     , 's', 1, "character", 'Enter experimental design file (see jupyter notebook tutorial)',
  'out'         , 'o', 1, "character", 'Enter base name for all output files'
), byrow=TRUE, ncol=5)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( is.null(opt$matrix) || is.null(opt$samples) || is.null(opt$out)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

matrixfile <- opt$matrix
countmatrix <- read.table(matrixfile, sep="\t", header=TRUE)
countmatrix <- round(countmatrix,0)
#head(countmatrix)

samplefile <- opt$samples
sampleTable <- read.table(samplefile, sep="\t", header=TRUE)
#head(sampleTable)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countmatrix, colData = sampleTable, design = ~ condition)
dds <- estimateSizeFactors(dds)
#head(counts(dds))

#nrow(dds)
dds_filt <- dds[ rowSums(counts(dds)) > 0, ]
#nrow(dds_filt)

normalized_matrix <- counts(dds_filt, normalized = TRUE)
#head(normalized_matrix)

vsd <- vst(dds_filt, blind = FALSE)
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
