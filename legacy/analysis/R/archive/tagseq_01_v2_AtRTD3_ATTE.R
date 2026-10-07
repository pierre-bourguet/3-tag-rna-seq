# import libraries ####
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbreak)
library(patchwork)
library(RColorBrewer)
library("pheatmap")
library(svglite)
library(ComplexHeatmap)
library(agricolae)
library(ggbeeswarm)

# import functions and parameters ####

source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/DEG_heatmap.R")
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/graphical_parameters.R")

#
# define outliers to remove ####
outliers <- c("empty_R1", # these two have low coverage
              "ddm1_a_long_b_2_R1", 
              "ab_2_R2", # this one is a clear outlier on euclidean distance heatmaps
              "b_2_R1", # these ones looks contaminated: have reads on up TEs
              "a_long_2_R2", 
              "a_long_1_R3", 
              "a_long_2_R1" # outlier on PCA of TEG reads
)

# remove mom1 and F2 segregants samples, not needed
outlier_patterns <- c("mom1", "ddm1_mom1", "F2_WT", "F2_a_2") # here one can define patterns: all samples matching will be discarded
#
# import and format counts ####

setwd("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE/")

# import counts, we don't use tximport here to avoid gene length normalization since we are only sequencing 3' ends
counts_sense <- as_tibble(read.delim("counts.tsv", header=T, sep="\t"))
counts_AS <- as_tibble(read.delim("counts_AS.tsv", header=T, sep="\t"))

# add suffix to Geneids and remove "_AS" suffix from column names
counts_AS$Geneid <- paste0(counts_AS$Geneid, "_AS")
names(counts_AS) <- gsub("_AS", "", names(counts_AS))

# remove outliers
counts_sense <- counts_sense %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns)))))
counts_AS <- counts_AS %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns)))))

# verify that column names match for sense and antisense count dataframes
if (!all(names(counts_sense) == names(counts_AS))) {
  stop("Column names do not match between sense and antisense count dataframes")
}

# rbind the two dataframes
counts <- rbind(counts_sense, counts_AS)

# Import annotations to filter counts only at TAIR10 ATTEs and PCGs
TEs <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/02_input/TAIR10_Transposable_Elements.txt", head=T, sep="\t", quote="", comment.char="", col.names=c("Geneid", "sense", "start", "end", "family", "superfamily")))
PCGs_pseudogenes <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/02_genome_files/03_output/TAIR10_GFF_PCG_pseudogene.tsv", head=T, sep="\t", quote="", comment.char=""))

# merge the counts at the gene level and remove transcript fusions
counts_merged <- counts %>%
  dplyr::mutate(Geneid = gsub("\\.\\d+", "", Geneid)) %>%  # Remove the .1, .2, etc. at the end of Geneid
  dplyr::group_by(Geneid) %>%
  dplyr::summarise(across(everything(), sum)) %>% # Sum the counts of all isoforms for each gene
  dplyr::filter(!str_detect(Geneid, "-")) %>% # Remove transcript fusions
  dplyr::filter(str_detect(Geneid, "^AT[1-5]")) %>% # only keep geneids from chromosome 1 to 5
  dplyr::filter(Geneid %in% TEs$Geneid | Geneid %in% paste0(TEs$Geneid, "_AS") | Geneid %in% PCGs_pseudogenes$Geneid | Geneid %in% paste0(PCGs_pseudogenes$Geneid, "_AS")) # only keep geneids from TEs and PCGs

# import salmon TPM normalized samples and merge by annotation ####
TPM_sense <- as_tibble(read.delim("TPM_salmon.tsv", header=T, sep="\t"))
TPM_AS <- as_tibble(read.delim("TPM_salmon_AS.tsv", header=T, sep="\t"))

# add suffix to geneids and remove "_AS" suffix from column names
TPM_AS$Geneid <- paste0(TPM_AS$Geneid, "_AS")
names(TPM_AS) <- gsub("_AS", "", names(TPM_AS))

# remove outliers
TPM_sense <- TPM_sense %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns)))))
TPM_AS <- TPM_AS %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns)))))

# verify that column names match for sense and antisense count dataframes
if (!all(names(TPM_sense) == names(TPM_AS))) {
  stop("Column names do not match between sense and antisense count dataframes")
}

# rbind the two dataframes
TPM <- rbind(TPM_sense, TPM_AS)

# merge the counts at the gene level, remove transcript fusions and geneids not on chromosome 1 to 5
TPM_merged <- TPM %>%
  dplyr::mutate(Geneid = gsub("\\.\\d+", "", Geneid)) %>%  # Remove the .1, .2, etc. at the end of Geneid
  dplyr::group_by(Geneid) %>%
  dplyr::summarise(across(everything(), sum)) %>% # Sum the TPM of all isoforms for each gene
  dplyr::filter(!str_detect(Geneid, "-")) %>% # Remove transcript fusions
  dplyr::filter(str_detect(Geneid, "^AT[1-5]")) # %>% # only keep geneids from chromosome 1 to 5
#filter(Geneid %in% TEs$GeneId | Geneid %in% paste0(TEs$GeneId, "_AS") | Geneid %in% PCGs$GeneId | Geneid %in% paste0(PCGs$GeneId, "_AS")) # only keep geneids from TEs and PCGs
#
# average the replicates ####

process_TPM <- function(df) {
  df %>%
    pivot_longer(
      cols = -Geneid,
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(genotype = str_sub(sample, 1, -4)) %>%
    group_by(Geneid, genotype) %>%
    summarise(avg_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(
      names_from = genotype,
      values_from = avg_value
    )
}

# Apply the function
TPM_merged_avg <- process_TPM(TPM_merged)

#
# import sample information ####
samples <- read.delim("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/03_sample_lists/DESeq2_tagseq_01_cdca7.tsv"
                      , header=T, sep=',', quote="", dec=".", comment.char="")[,c(2,3)]
samples <- samples[order(samples[,1], samples[,2]),] # reorder samples by alphabetical names

# prepare a full name to remove outliers
samples$full_name <- paste0(samples$condition, "_", samples$sample)
# remove outliers
samples <- samples %>%
  dplyr::filter(
    !full_name %in% outliers &
      !str_detect(full_name, paste(outlier_patterns, collapse = "|"))
  )

sample_columns <- which(names(counts_merged) %in% paste(samples$condition, samples$sample, sep="_")) # retrieve columns that contain sample names

# find all conditions different from reference
reference <- "Col_0"
conditions <- unique(samples$condition)[unique(samples$condition) != reference]

#
# DESeq2 function ####
library("DESeq2")
DESeq2_function <- function(x) { # x should be cts_summary.tsv file (summary of raw counts)
  cts <- ceiling(x[,sample_columns]) ; row.names(cts) <- x$Geneid # ceiling is to round up
  # defining metadata
  coldata <- data.frame(
    condition=gsub("_R.$", "", names(x)[sample_columns]),
    type=rep("single-strand", nrow(samples))
  )
  row.names(coldata) <- names(x)[sample_columns]
  if (!all(rownames(coldata) == colnames(cts))) {
    print("the names in count_file and the sample table do not match")
  } # IF FALSE YOU ARE IN DEEP SHIT MY MAN, NOT GONNA WORK
  dds <<- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ condition)
  # adding meta data to the dataframe
  mcols(dds) <- DataFrame(mcols(dds))
  # pre-filtering, here keeping only genes with at least 10 reads in at least 3 samples
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  # setting the reference treatment
  dds$condition <- relevel(dds$condition, ref = reference)
  # differential analysis
  return(DESeq(dds))
}
dds <- DESeq2_function(counts_merged)

# save the deseq2 object
save(dds, file = "DEGs_batch/deseq2_object.RData")
load("DEGs_batch/deseq2_object.RData")

# create output directory for plots and tables
output_dir <- paste0("DEGs_batch/")
ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
#
# quality controls: PCAs ####

#### PCAs

# vst transformation
vsd <- vst(dds, blind=FALSE)

# extract PCA data to filter samples easily
ntop_variable_features <- 1000 # number of most variable features for PCA
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE, ntop=ntop_variable_features)
pcaData <- pcaData %>%
  mutate(ddm1 = str_detect(name, "ddm1"))
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Define a range of shapes to use
set.seed(123)
shape_list <- sample(15:25, length(unique(pcaData$condition)), replace = TRUE) # Using a set of distinct shapes

## all samples
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable features")) +
  coord_fixed() +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal()

# without ddm1
ggplot(pcaData %>% filter(ddm1 == FALSE), aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable features")) +
  coord_fixed() +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal()

# ddm1 samples only
ggplot(pcaData %>% filter(ddm1 == TRUE), aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable features")) +
  coord_fixed() +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal() # Optional: use a minimal theme for a cleaner look

# plotting PCA on transposable elements #### 

## default plotting (useful as it shows the % of PC1 & PC2)
# subset TEs
counts_merged_subset <- counts_merged %>%
  filter(Geneid %in% TEs$Geneid)

vsd_TEG <- vst(DESeq2_function(counts_merged_subset), blind=FALSE, nsub=100) # you might have to lower nsub to lower than default (1000) if there are too few reads at TEs
# extract PCA data & create a new column to plot only some samples
ntop_variable_features <- 500
pcaData <- plotPCA(vsd_TEG, intgroup=c("condition", "type"), returnData=TRUE, ntop = ntop_variable_features) %>%
  mutate(ddm1 = str_detect(name, "ddm1"))
# extract % variance for each PC
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Define a range of shapes to use
shape_list <- sample(15:25, length(unique(pcaData$condition)), replace = TRUE) # Using a set of distinct shapes

# Plot all samples
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable TEs")) +
  coord_fixed() +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal() # Optional: use a minimal theme for a cleaner look

# exclude ddm1 samples
ggplot(pcaData %>% filter(ddm1 == F), aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable TEs")) +
  coord_fixed() +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal() # Optional: use a minimal theme for a cleaner look

#pdf(paste0(output_dir, "PCA_TEG.pdf"))
#plotPCA(vsd_TEG, intgroup=c("condition"))
#graphics.off()

# with extra color and styling: all TEs
#pdf(paste0(output_dir, "PCA_TE_colored.pdf"))
ggplot(pcaData, aes(x=PC1, y=PC2, fill=condition)) +
  geom_point(pch=21, size=3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable TEs")) +
  scale_fill_manual(values=many_colors) +
  theme_minimal()
#graphics.off()

## now make a barplot of PC1 values

# Create an ordered dataframe
data <- as_tibble(plotPCA(vsd_TEG, intgroup = c("condition"), returnData = TRUE)) %>%
  mutate(condition = as.factor(condition)) %>% # Make sure condition is a factor
  arrange(desc(PC1)) # Arrange by PC1 to get the order

# barplot
ggplot(data
       , aes(x = PC1, y = reorder(name, PC1), fill = condition)) +
  geom_bar(stat = "identity", orientation = "y", colour="black") +
  labs(title = "Ordered PC1 Values", x = "PC1", y = "Condition") +
  scale_fill_manual(values=c(col_vibrant, col_high_contrast, col_bright[-7], col_muted)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) # Ensure y-axis labels are readable

# barplots of PC1 for mutant lines
# ggplot(data %>% filter(str_detect(condition, "cdca7_ab_"))
#        , aes(x = PC1, y = reorder(name, PC1), fill = condition)) +
#   geom_bar(stat = "identity", orientation = "y", colour="black") +
#   labs(title = "Ordered PC1 Values", x = "PC1", y = "Condition") +
#   scale_fill_manual(values=c(col_vibrant, col_high_contrast, col_bright[-7], col_muted)) +
#   theme_minimal() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1)) # Ensure y-axis labels are readable


# quality controls: heatmaps ####

# all samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- names(vsd$sizeFactor)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# show samples without ddm1 mutations
sampleDists <- dist(t(
  assay(vsd)[, grep("ddm1", colnames(assay(vsd)), invert = TRUE) ]
))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- grep("ddm1", names(vsd$sizeFactor), invert = T, value = T)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# keep samples with ddm1 mutations
sampleDists <- dist(t(
  assay(vsd)[, grep("ddm1", colnames(assay(vsd))) ]
))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- grep("ddm1", names(vsd$sizeFactor), value = T)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# same but removing mom1 ddm1 which are actually just mom1: doesn't change the interpretation
# remove samples without ddm1 mutations
x <- assay(vsd)[, grep("ddm1_[^mom1]", colnames(assay(vsd)))] # this is a negative lookahead to exclude ddm1_mom1
sampleDists <- dist(t(x))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- grep("ddm1_[^mom1]", names(vsd$sizeFactor), value = T)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# import plate well position to see if position correlates with batch effects
sample_wells <- read_tsv("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/03_sample_lists/sample_list_tagseq_01_cdca7.tsv"
                         , col_names = FALSE, show_col_types = FALSE)

# Extract the 'well' information
sample_wells <- sample_wells %>%
  mutate(
    sample = X2,
    well = sub(".*/([^/]+)_(.*)\\..*", "\\1", X1)
  ) %>%
  dplyr::select(sample, well) %>%
  mutate(row = substr(well, 1, 1),
         column = as.numeric(substr(well, 2, nchar(well))))

#### plot sample to sample correlation with well position in the plate to control for batch effects
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/plot_heatmap_plate_batch_effect.R")

# sample distance matrix with all samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- names(vsd$sizeFactor)
colnames(sampleDistMatrix) <- NULL
# heatmap
pdf(paste0(output_dir, "euclidean_distance_heatmap_all_samples.pdf"), width=17, height=15)
create_heatmap_with_annotations(sampleDistMatrix, sample_wells)
dev.off()

# without ddm1 mutations
sampleDists <- dist(t(
  assay(vsd)[, grep("ddm1", colnames(assay(vsd)), invert = TRUE) ]
))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- grep("ddm1", names(vsd$sizeFactor), invert = T, value = T)
colnames(sampleDistMatrix) <- NULL
# heatmap
pdf(paste0(output_dir, "euclidean_distance_heatmap_no_ddm1_samples.pdf"), width=17, height=15)
create_heatmap_with_annotations(sampleDistMatrix, sample_wells)
dev.off()

# with ddm1 samples only
x <- assay(vsd)[, grep("ddm1", colnames(assay(vsd)))]
sampleDists <- dist(t(x))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- grep("ddm1", names(vsd$sizeFactor), value = T)
colnames(sampleDistMatrix) <- NULL
# heatmap
pdf(paste0(output_dir, "euclidean_distance_heatmap_ddm1_samples.pdf"), width=17, height=15)
create_heatmap_with_annotations(sampleDistMatrix, sample_wells)
dev.off()

#
# compare all mutants with control ####

f <- function(aa, bb) { # this creates a res_mutant dataframe comparing mutant vs control
  eval(substitute( a <- results(dds, contrast=c("condition",as.character(b), reference))
                   , list(a = aa, b = bb))) # this second argument provides an environment for substitute, defining variables used in previous line
}
all_res <- Map(f, paste0("res_", conditions), as.list(conditions)) # list of all mutant vs WT comparisons
#
# creates list of all up & down TEs & PCGs (across all mutant conditions) #### 
up <- function(x, y) {
  geneids <- c(eval(parse(text=paste0(y, "$Geneid"))), paste0(eval(parse(text=paste0(y, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  return( row.names(x[row.names(x) %in% geneids & x$log2FoldChange >=1 & !is.na(x$padj) & x$padj < 0.1,]) )
}
down <- function(x, y) {
  geneids <- c(eval(parse(text=paste0(y, "$Geneid"))), paste0(eval(parse(text=paste0(y, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  return( row.names(x[row.names(x) %in% geneids & x$log2FoldChange <=-1 & !is.na(x$padj) & x$padj < 0.1,]) )
}
DEGs <- list(
  upTEs = unique(unlist(lapply(FUN = up, X = all_res, y="TEs"))),
  upPCGs = unique(unlist(lapply(FUN = up, X = all_res, y="PCGs_pseudogenes"))),
  downTEs = unique(unlist(lapply(FUN = down, X = all_res, y="TEs"))),
  downPCGs = unique(unlist(lapply(FUN = down, X = all_res, y="PCGs_pseudogenes")))
)

lapply(DEGs, length)
#
# rlog normalization ####
rld <- rlog(dds, blind=FALSE)
rld_df <- as_tibble(assay(rld))
rld_df$Geneid <- row.names(assay(dds))

# long format
long_rld <- rld_df %>%
  pivot_longer(cols = -Geneid, names_to = "sample", values_to = "expression") %>%
  mutate(condition = str_replace(sample, "_R[123]", ""))
# average expression
long_rld_avg <- long_rld %>%
  group_by(Geneid, condition) %>%
  summarise(average_expression = mean(expression), .groups = 'drop')
# wide format with average to export
wide_rld_avg <- long_rld_avg %>%
  pivot_wider(names_from = condition, values_from = average_expression)

# export
DEG_write_rlog <- function(x, y) {
  dir.create("DEGs_batch/normalized_counts/", showWarnings = F, recursive = T)
  write.table(subset(wide_rld_avg, subset=wide_rld_avg$Geneid %in% x), file=paste0("DEGs_batch/normalized_counts/batch_", y, "_rlog.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
}
mapply(FUN = DEG_write_rlog, x=DEGs, y=names(DEGs))

# write tables
write.table(x=rld_df, file=paste0("DEGs_batch/normalized_counts/rlog.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=wide_rld_avg, file=paste0("DEGs_batch/normalized_counts/rlog_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# write DEG tables with TPM normalized counts ####
DEG_write_TPM <- function(x, y) {
  dir.create("DEGs_batch/normalized_counts", showWarnings = F, recursive = T)
  write.table(subset(TPM_merged, subset=TPM_merged$Geneid %in% x), file=paste0("DEGs_batch/normalized_counts/batch_", y, "_TPM.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
}
mapply(FUN = DEG_write_TPM, x=DEGs, y=names(DEGs))
# number of batch DEGs
write.table(x=t(as.data.frame(lapply(DEGs, FUN=length))), file=paste0(output_dir, "number_of_DEGs_batch.tsv"), quote = F, sep="\t", row.names=T, col.names=F)

#
# plot number of up TEs in cdca7 mutants: fig1 ####

#### up TEs

up_TEs <- lapply(FUN = length, X = (lapply(FUN = up, X = all_res, y="TEs")))
up_TEs_df <- data.frame(
  Category = gsub("res_", "", names(up_TEs)),
  Count = unlist(up_TEs)
)
ggplot(up_TEs_df, aes(y = Category, x = Count)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  theme_minimal()

# filter out samples
up_TEs_df <- up_TEs_df %>%
  filter(!str_detect(Category, "mom1"))
up_TEs_df_no_ddm1 <- up_TEs_df %>% filter(!str_detect(Category, "ddm1_[^2]|F2")) %>%
  filter(!str_detect(Category, "ddm1_2_G5"))

# reorder levels
up_TEs_df_no_ddm1$Category <- factor(up_TEs_df_no_ddm1$Category, levels = rev(c("a_1", "a_2", "b_1", "b_2", "a_long_1", "a_long_2", "a_long_3", "ab_1", "ab_2", "a_long_b", "ddm1_2_G2")))

# create output directory
dir.create("DEGs_batch/figures", showWarnings = F, recursive = T)

# plot
fig1_up_TEs_hist <- ggplot(up_TEs_df_no_ddm1 %>% dplyr::filter(!str_detect(Category, "long"))
                          , aes(y = Category, x = Count, fill = Category)) +
  geom_col() +
  labs(x = "Count", y = "Mutant") +
  scale_fill_manual(values = rev(col_muted_2_replicates[c(1:4,7,8,5)])) +
  scale_x_continuous(breaks = c(0,100,200,300), limits=c(0,1510), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(350, 1350), ticklabels=c(1400, 1500), scales=0.3, space = 0.1, expand = F) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", linewidth = pt_0.5_to_mm),
    axis.text.y = element_markdown(size = 6),
    axis.title.x = element_text(margin = margin(t = -5)), # reduce space between axis titles and axis labels
    axis.title.y = element_text(margin = margin(r = -5)),
    plot.margin = margin(t = -5, r = -5)
    , legend.position = "none"
    , panel.grid.minor.x = element_blank()
    , panel.grid.major.y = element_blank()
  ) ; fig1_up_TEs_hist

# export in svg, works fine but fonts are not exported properly (open the svg file in a text editor to see).
svglite::svglite(filename = "DEGs_batch/figures/up_TEs_cdca7_svglite.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
fig1_up_TEs_hist + theme_horizontal_nature
dev.off()

#### hacking the fonts by modifying the svg file
source("DEGs_batch/figures/hack_svg_fonts.R")
input_file <- "DEGs_batch/figures/up_TEs_cdca7_svglite.svg"
output_file <- modify_svg_fonts(
  input_file,
  new_style = "ArialMT",
  replacements = list(
    list(c("α", "β"), "SymbolMT"),
    list(c("cdca7", "ddm1", "-1", "-2"), "Arial-ItalicMT")
  ),
  resize_dashes = TRUE  # Set to TRUE or FALSE as needed
)

# add a WT category
fig1_up_TEs_df <- up_TEs_df_no_ddm1 %>%
  filter(!str_detect(Category, "long")) %>%
  bind_rows(tibble(Category = "WT", Count = NA))

fig1_up_TEs_df$Category <- factor(fig1_up_TEs_df$Category, levels = rev(c("WT", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))

fig1_up_TEs_hist <- ggplot(fig1_up_TEs_df, aes(y = Category, x = Count, fill = Category)) +
  geom_col() +
  labs(x = "Count", y = "Genotype") +
  scale_y_discrete(labels = c(italicized_labels_DEGs, "WT"), position = "left") +
  scale_fill_manual(values = rev(col_muted_2_replicates[c(19,1:4,7,8,5)])) +
  xlim(0,950) +
  scale_x_continuous(breaks = c(0,100,200), limits=c(0,930), expand = c(0, 0)) +
  scale_x_break(c(230, 870), ticklabels=c(800, 900), scales=0.3, space = 0.1, expand = F) +
  theme_horizontal_nature +
  theme(
    panel.border = element_rect(colour = "black", linewidth = pt_0.5_to_mm),
    axis.text.y = element_markdown(size = 6),
    axis.title.x = element_text(margin = margin(t = -5)), # reduce space between axis titles and axis labels
    axis.title.y = element_text(margin = margin(r = -5)),
    plot.margin = margin(t = -5, r = -5),
    legend.position = "none"
  ) ; fig1_up_TEs_hist

# export in svg, works fine but fonts are not exported properly (open the svg file in a text editor to see).
svglite::svglite(filename = "DEGs_batch/figures/up_TEs_cdca7_histogram.svg", width = 60*mm_to_inches, height = 50*mm_to_inches)
fig1_up_TEs_hist
dev.off()

#### hacking the fonts by modifying the svg file
source("DEGs_batch/figures/hack_svg_fonts.R")
input_file <- "DEGs_batch/figures/up_TEs_cdca7_histogram.svg"
output_file <- modify_svg_fonts(
  input_file,
  new_style = "ArialMT",
  replacements = list(
    list(c("α", "β"), "SymbolMT"),
    list(c("cdca7", "ddm1", "-1", "-2"), "Arial-ItalicMT")
  ),
  resize_dashes = TRUE  # Set to TRUE or FALSE as needed
)

#
# plot number of up PCGs in cdca7 mutants ####

up_PCGs <- lapply(FUN = length, X = (lapply(FUN = up, X = all_res, y="PCGs")))
up_PCGs_df <- data.frame(
  Category = gsub("res_", "", names(up_PCGs)),
  Count = unlist(up_PCGs)
)
ggplot(up_PCGs_df, aes(y = Category, x = Count)) +
  geom_col() +
  labs(title = "Number of upregulated PCGs in each genotype", x = "Count", y = "Genotype") +
  theme_minimal()

# remove mom1
up_PCGs_df <- up_PCGs_df %>%
  filter(!str_detect(Category, "mom1"))

up_PCGs_df_no_ddm1 <- up_PCGs_df %>% filter(!str_detect(Category, "ddm1_[^2]|F2")) %>%
  filter(!str_detect(Category, "ddm1_2_G5"))
up_PCGs_df_no_ddm1$Category <- factor(up_PCGs_df_no_ddm1$Category, levels = rev(c("a_1", "a_2", "b_1", "b_2", "a_long_1", "a_long_2", "a_long_3", "ab_1", "ab_2", "a_long_b", "ddm1_2_G2")))

ggplot(up_PCGs_df_no_ddm1 %>% filter(!str_detect(Category, "ddm1|long")), aes(y = Category, x = Count)) +
  geom_col() +
  labs(title = "Number of upregulated PCGs in each genotype", x = "Genotype", y = "Count") +
  scale_y_discrete(labels = italicized_labels_DEGs) +
  theme_minimal()

# plot number of up TEs in ddm1 cdca7 mutants ####

up_TEs <- lapply(FUN = length, X = (lapply(FUN = up, X = all_res, y="TEs")))
up_TEs_df <- data.frame(
  condition = gsub("res_", "", names(up_TEs)),
  Count = unlist(up_TEs)
)

# filter samples
up_TEs_df_w_ddm1 <- up_TEs_df %>%
  # filter(str_detect(condition, "ddm1|ab")) %>% # use this one include cdca7ab controls
  filter(str_detect(condition, "ddm1")) %>%
  filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
  filter(!condition %in% "F2_ddm1") %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
up_TEs_df_w_ddm1$condition <- factor(up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_ddm1 <- ggplot(up_TEs_df_w_ddm1, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,4,3,5)], rep(col_muted[5], 5)))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_ddm1

#
# cdca7 long mutants: number of up TEs, heatmap and superplots ####

up_TEs <- lapply(FUN = length, X = (lapply(FUN = up, X = all_res, y="TEs")))
up_TEs_df <- data.frame(
  condition = gsub("res_", "", names(up_TEs)),
  Count = unlist(up_TEs)
)

# filter samples
up_TEs_df_a_long_no_ddm1 <- up_TEs_df %>% filter(str_detect(condition, "long|a_1|a_2|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

up_TEs_df_a_long_only_ab <- up_TEs_df %>% filter(str_detect(condition, "a_long_b|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0)) %>%
  mutate(condition = factor(condition, levels = rev(c("Col_0", "ab_1", "ab_2", "a_long_b")))) %>%
  arrange(condition)
  
up_TEs_df_a_long_w_ddm1 <- up_TEs_df %>% filter(str_detect(condition, "ddm1") & str_detect(condition, "G5|long|ab")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
# up_TEs_df_w_ddm1$condition <- factor(up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_a_long <- ggplot(rbind(up_TEs_df_a_long_no_ddm1, up_TEs_df_a_long_w_ddm1), aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  #scale_fill_manual(values = rev(c(col_muted[c(10,4,4,3,5)], rep(col_muted[5], 5)))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long

barplot_nb_up_TEs_a_long_b <- ggplot(up_TEs_df_a_long_only_ab, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,4,6)]))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long_b

## heatmap at TEs up in cdca7a/b
cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7a_long_mutants"
            , z = rld_df %>% dplyr::select(matches("long_b|ab|Geneid") & -matches("ddm1"))
            , n = "rlog")
wide_rld_avg

# superplot
#cdca7_ab_up_TEs_union <- union(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))

a_long_b_rld <- rld_df %>%
  dplyr::filter(Geneid %in% cdca7_ab_up_TEs) %>%
  dplyr::select(matches("long_b|ab|Geneid|Col_0") & -matches("ddm1")) %>%
  tidyr::pivot_longer(cols = -Geneid, names_to = "sample", values_to = "value") %>%
  dplyr::mutate(condition=stringr::str_replace(sample, "_R[123]", ""))

superplot_a_long_b <- superplot_w_boxplot(data = a_long_b_rld, condition_order = rev(c("Col_0", "ab_1", "ab_2", "a_long_b"))
                                          , colors = rev(c("grey55", col_muted[c(4,4,6)]))) + 
  coord_flip(ylim=c(1.2,9)) +
  xlab("Genotype") + ylab("Transcript levels (rlog)") + labs(title="upregulated TEs")

# write svg output
svglite::svglite(filename = "DEGs_batch/figures/superplot_a_long_b.svg", width = 3, height = 1.5)
barplot_nb_up_TEs_a_long_b + superplot_a_long_b + plot_layout(axes = "collect", widths = c(2,3)) #& theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- a_long_b_rld %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
# heatmap of upregulated TEs ####
upTEs_TPM <- as_tibble(read.delim("DEGs_batch/normalized_counts/batch_upTEs_TPM.tsv", header=T, sep="\t"))
upTEs_rlog <- as_tibble(read.delim("DEGs_batch/normalized_counts/batch_upTEs_rlog.tsv", header=T, sep="\t"))

# filter out samples with ddm1 mutations, mom1 mutants, F2 segregants
upTEs_TPM_no_ddm1 <- upTEs_TPM[, grep("ddm1|F2|mom1", colnames(upTEs_TPM), invert = T)]
upTEs_rlog_no_ddm1 <- upTEs_rlog[, grep("ddm1|F2|mom1", colnames(upTEs_rlog), invert = T)]

# heatmap of upregulated TEs
normalization <- "TPM"
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=TPM_merged, n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=rld_df, n = "rlog") ) # DOESN'T WORK FOR SOME REASON

### TEs up regulated in cdca7a/b or ddm1

cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
ddm1_G2_up_TEs <- up(all_res$res_ddm1_2_G2, "TEs")

## using log2 (TPM + 1) & rlog at cdca7 upTEs

# with all samples
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs", z = TPM_merged, n = normalization)

# only cdca7 mutants
# log2(TPM+1)
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants", z = TPM_merged %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = normalization)
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants_avg", z = TPM_merged_avg %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = normalization)
# rlog
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants", z = rld_df %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = 'rlog')
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants_avg", z = wide_rld_avg %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = 'rlog')

## heatmap of rlog at cdca7-upTEs for fig1
# Specify the normalization method
normalization_method <- "rlog"  # "TPM" or "rlog"

# Prepare the matrices and other variables based on the normalization method
if (normalization_method == "TPM") {
  z <- TPM_merged_avg %>%
    dplyr::select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
  # Select up TEs & log2 transform & remove Geneid
  m <- as.matrix(log2(z[z$Geneid %in% cdca7_ab_up_TEs, -1] + 1))
  color_scale <- colorRamp2(seq(from = 0, to = 9, by = 1), scico(n = 10, direction = -1, palette = "lajolla"))
  heatmap_name <- "log2\n(TPM+1)"
} else if (normalization_method == "rlog") {
  z <- wide_rld_avg %>%
    dplyr::select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
  # Select up TEs & remove Geneid
  m <- as.matrix(z[z$Geneid %in% cdca7_ab_up_TEs, -1])
  color_scale <- colorRamp2(seq(from = 3, to = 10, by = 1), scico(n = 8, direction = -1, palette = "lajolla"))
  heatmap_name <- "rlog"
}

# See data spread to define the color scale
summary(m)

# Prepare labels
italicized_labels_heatmap <- c(
  expression(plain('WT')),
  expression(italic('cdca7α-1')),
  expression(italic('cdca7α-2')),
  expression(italic('cdca7β-1')),
  expression(italic('cdca7β-2')),
  expression(italic('cdca7α/β-1')),
  expression(italic('cdca7α/β-2')),
  expression(italic('ddm1'))
)

# Plot the heatmap
Heatmap(m,
        name = heatmap_name,
        col = color_scale,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width = 15, height = 15,
        column_title = paste0("cdca7α/β upTEs\nn=", length(cdca7_ab_up_TEs)),
        column_labels = italicized_labels_heatmap
)

# remove clustering & sort the heatmap by row sums (TEs with highest expression across genotypes at the top)
m <- m[order(rowSums(m), decreasing = T),] 
fig1_up_TEs_heatmap <- Heatmap(t(m),
                                name = heatmap_name,
                                col = color_scale,
                                cluster_rows = F,
                                cluster_columns = F,
                                width=15, height=15,
                                column_title=NULL,
                                row_labels = italicized_labels_heatmap,
                                row_names_side = "left",
                                row_names_gp = grid::gpar(fontsize = 6),
                                use_raster = T,
                                border= T,
                                heatmap_legend_param = list(
                                  title = "log2\n(TPM+1)", at = c(0, 5, 10), 
                                  labels = c("0", "5", "10"),
                                  legend_height = unit(1, "cm"),
                                  legend_width = unit(1, "cm"),
                                  labels_gp = gpar(fontsize = 6),
                                  title_gp = gpar(fontsize = 6)
                                )
) ; fig1_up_TEs_heatmap

svglite::svglite(filename = "DEGs_batch/figures/up_TEs_cdca7_heatmap.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
svglite::svglite(filename = "DEGs_batch/figures/up_TEs_cdca7_heatmap_rlog.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
draw(fig1_up_TEs_heatmap)
dev.off()

## using log2 (TPM + 1) at ddm1 upTEs

# select samples of interest
z <- TPM_merged_avg %>%
  select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
# log2 transform
m <- as.matrix( log2(z[z$Geneid %in% ddm1_G2_up_TEs,-1] + 1))
# see data spread to scale the heatmap
summary(m)
col_log2_TPM <- colorRamp2(seq(from=0, to=9, by=1), scico(n=10, direction=-1, palette="lajolla"))
Heatmap(m,
        name = "log2\n(TPM+1)",
        col = col_log2_TPM,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15,
        column_title=paste0("ddm1-2 upTEs\nn=", length(ddm1_G2_up_TEs))
)

## Z-score heatmap, scaling on TPM
z <- TPM_merged_avg %>%
  select(-matches("ddm1|F2|mom1|long"))
col_order <- colnames(m)[c(1:3,6,7,4,5)]

m <- as.matrix(z[z$Geneid %in% cdca7_ab_up_TEs,-1])
m <- t(scale(t(m)))
summary(m)
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
Heatmap(m,
        name = "z-score",
        col = col_zscore,
        column_order = col_order, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15
)

# z-score scaling after log2 + 1 transformation
m <- as.matrix( log2(z[z$Geneid %in% cdca7_ab_up_TEs,-1] + 1))
m <- t(scale(t(m)))
summary(m)
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
Heatmap(m,
        name = "z-score",
        col = col_zscore,
        column_order = col_order, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15
)

# Save the heatmap to a PDF file
ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
file_name <- paste0(output_dir, n, "/heatmap_", y, ".pdf")
pdf(file_name, width = 15, height = 7)
draw(ht, heatmap_legend_side = "right")
dev.off()

# heatmap at ddm1 up TEs
ddm1_G5_up_TEs <- up(all_res$res_ddm1_2_G5, "TEs")
# with all samples
DEG_heatmap(ddm1_G5_up_TEs, "ddm1_G5", z = TPM_merged, n = normalization)

#
# superplot at TEs up in any sample ####

## superplot using log2 (TPM + 1) or rlog at TEs up in any sample

# Specify the normalization method
normalization_method <- "rlog"  # "TPM" or "rlog"

# Prepare the matrices and other variables based on the normalization method
if (normalization_method == "TPM") {
  z <- TPM_merged %>%
    tidyr::pivot_longer(-Geneid, names_to = "sample", values_to = "expression") %>%
    dplyr::mutate(condition = gsub("_R[0-9]", "", sample)) %>%
    dplyr::filter(str_detect(condition, "ddm1|Col_0|ab")) %>%
    dplyr::filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
    dplyr::filter(!condition %in% "F2_ddm1")
  # Select up TEs & log2 transform
  m <- z %>%
    dplyr::filter(Geneid %in% DEGs$upTEs) %>%
    dplyr::mutate(expression = log2(expression + 1))
} else if (normalization_method == "rlog") {
  z <- long_rld %>%
    #dplyr::filter(str_detect(condition, "ddm1|Col_0|ab")) %>% # use this one to include Col-0 & ab controls
    dplyr::filter(str_detect(condition, "ddm1|Col_0")) %>% # use this one to include Col-0 only
    dplyr::filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
    dplyr::filter(!condition %in% "F2_ddm1")
  # Select up TEs
  m <- z[z$Geneid %in% DEGs$upTEs,]
}

# import the superplot function
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/superplot_w_boxplot.R")

# Define the specific order vector
order_vector <- c("Col_0", "ab_1", "ab_2" , "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")
order_vector <- c("Col_0", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")

# Define the color palette
my_colors <- c("grey55", col_muted[c(4,4,3)], rep(col_muted[5], 6)) # to include cdca7ab controls
my_colors <- c("grey55", col_muted[c(3)], rep(col_muted[5], 6)) # without
# rename the dataframe for compatibility with the function
names(m)[3] <- "value"

# run the function and customize the plot
superplot_upTEs_ddm1 <- superplot_w_boxplot(m, rev(order_vector), rev(my_colors)) + 
  coord_flip(ylim=c(0,10)) +
  xlab("Genotype") + ylab("Transcript levels (rlog)") + labs(title="upregulated TEs")
superplot_upTEs_ddm1

barplot_nb_up_TEs_ddm1 + superplot_upTEs_ddm1 + patchwork::plot_layout(axes = 'collect')

svglite::svglite(filename = "DEGs_batch/figures/fig_sup_TEs_ddm1_cdca7.svg", width = 80*mm_to_inches, height = 50*mm_to_inches)
barplot_nb_up_TEs_ddm1 + superplot_upTEs_ddm1 + plot_layout(axes = 'collect') & theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- m %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

library(agricolae)
# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)


#
# looking at TEs ####

# all upregulated TEs at all samples

d <- TPM %>%
  filter(Geneid %in% DEGs$upTEs) %>%
  pivot_longer(cols=-Geneid) %>%
  mutate(condition=substr(name, 1, nchar(name)-3))  %>%
  mutate(log2_TPM = log2(value + 1)) %>% # distinguish ddm1 conditions from others
  mutate(ddm1 = if_else(str_detect(condition, "ddm1"), TRUE, FALSE))

ggplot(d %>% filter(ddm1 == T), aes(x = name, y = log2_TPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = paste0("quantification of TE transcripts\nn = ", length(unique(d$Geneid))), y = "Log2 (TPM + 1)", x = "Genotype") +
  coord_flip() +
  theme(legend.position = "none")

ggplot(d %>% filter(ddm1 == F), aes(x = name, y = log2_TPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = paste0("quantification of TE transcripts\nn = ", length(unique(d$Geneid))), y = "Log2 (TPM + 1)", x = "Genotype") +
  coord_flip() +
  theme(legend.position = "none") +
  ylim(0,8)

# Calculating the median of log2_TPM for each name and merge with the transgene expression levels

medians <- d %>%
  group_by(name, condition) %>%
  summarise(median_log2_TPM = median(log2_TPM, na.rm = TRUE), .groups = 'drop')

# merge with transgene expression levels
mTurq_data <- TPM %>%
  filter(Geneid %in% "mTurq_3xcMyc") %>%
  pivot_longer(cols=-Geneid) %>%
  mutate(condition=substr(name, 1, nchar(name)-3)) %>%
  mutate(log2_TPM = log2(value + 1)) %>%
  select(name, mTurq_log2_TPM = log2_TPM)

medians <- left_join(medians, mTurq_data, by = c("name"))

# plot all data
ggplot(medians, aes(y = condition, x = median_log2_TPM, fill = condition, size = mTurq_log2_TPM)) +
  geom_jitter(width = 0.01, height = 0.2, shape = 21, show.legend = F) +
  theme_minimal(base_family = "Arial") +  
  theme(text = element_text(family = "Arial")) +
  labs(title = paste0("upregulated TE transcripts\nn = ", length(DEGs$upTEs)),
       y = "", x = "Median log2 (TPM + 1)", fill="genotype")# +
scale_fill_manual(values=col_muted[c(4,2,1,10)]) +
  scale_y_discrete(labels= rev(italicized_labels_tagseq)) +
  theme(axis.title.y = element_blank())


# all upregulated TEs in samples of interest ####

samples <- c("WT", "cdca7_a", "cdca7_b", "cdca7_ab")

dd <- d %>%
  filter(condition %in% samples)

ggplot(dd, aes(x = name, y = log2_TPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = "quantification of TE transcripts", y = "log2 (TPM + 1)", x = "genotype") +
  coord_flip()

# now plotting the median of replicates
# Filtering and ordering data based on specific samples
dd <- d %>%
  filter(condition %in% samples) %>%
  mutate(condition = factor(condition, levels = rev(samples)))

# Calculating the median of log2_TPM for each name
medians <- dd %>%
  group_by(name, condition) %>%
  summarise(median_log2_TPM = median(log2_TPM, na.rm = TRUE), .groups = 'drop')

# Creating the plot with ordered conditions
italicized_labels_tagseq <- c(
  expression(plain('WT')),
  expression(italic('cdca7-α')),
  expression(italic('cdca7-β')),
  expression(italic('cdca7-α/β'))
)
grDevices::cairo_pdf("figures/CDCA7-ab_log2_TPM_at_TEs.pdf", width = 4, height = 2)
set.seed(130)
p2 <- ggplot(medians, aes(y = condition, x = median_log2_TPM, fill = condition)) +
  geom_jitter(width = 0.01, height = 0.2, shape = 21, size = 3, show.legend = F) +
  theme_minimal(base_family = "Arial") +  
  theme(text = element_text(family = "Arial")) +
  labs(title = paste0("upregulated TE transcripts\nn = ", length(DEGs$upTEs)),
       y = "", x = "Median log2 (TPM + 1)", fill="genotype") +
  scale_fill_manual(values=col_muted[c(4,2,1,10)]) +
  scale_y_discrete(labels= rev(italicized_labels_tagseq)) +
  theme(axis.title.y = element_blank()) ; p2
dev.off()

library(patchwork)
setwd("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/tagseq_03_cdca7/")
set.seed(130)
magnifier <- 1.5 ; grDevices::cairo_pdf("figures/poster_fig2A_2B.pdf", width = 5*magnifier, height = 4*magnifier)
p1 + p2 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold"),
                                                    text = element_text(size=14, family="Arial"),
                                                    axis.text.y = element_text(size=16)
)
dev.off()



# make linear correlations of rlogs between samples ####

### between cdca7-ab mutants
cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
cdca7_ab_up_PCGs <- intersect(up(all_res$res_ab_1, "PCGs"), up(all_res$res_ab_2, "PCGs"))
ddm1_up_PCGs <- up(all_res$res_ddm1_2_G2, "PCGs")
cdca7_ab_rld <- wide_rld_avg %>%
  dplyr::filter(Geneid %in% cdca7_ab_up_TEs)

# linear correlation between samples
cor(cdca7_ab_rld %>% dplyr::select(-Geneid), method = "pearson")
# scatterplot of cdca7-ab mutants
ggplot(cdca7_ab_rld, aes(x = ab_1, y = a_long_b)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  geom_abline() +
  labs(title = "cdca7-ab mutants", x = "rlog(ab_1)", y = "rlog(ab_2)") +
  coord_fixed(ylim=c(2,11), xlim=c(2,11))

data <- cdca7_ab_rld %>% dplyr::select(c("Col_0", "ab_1", "a_long_b", "ab_2", "ddm1_2_G2", "ddm1_2_G5"))
data <- cdca7_ab_rld %>% dplyr::select(c("Col_0", "a_1", "a_2", "a_long_1", "a_long_2", "a_long_3"))
pairs(data)
library(GGally)
ggpairs(data, title = "Scatter Plot Matrix for mtcars Dataset", axisLabels = "show") +
  coord_fixed(ylim=c(1,11), xlim=c(1,11)) +
  geom_abline(intercept = 0, slope =1)


#
### NOT UPDATED from there ####

grDevices::cairo_pdf("figures/CDCA7-ab_log2_TPM_at_TEs.pdf", width = 4, height = 2)
set.seed(130)

# plot heatmaps ####
suppressPackageStartupMessages(library("RColorBrewer")) ; suppressPackageStartupMessages(library("pheatmap"))
DEG_heatmap <- function(x, y, z, n) {
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      if (n=="rlog") {
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\n", n)
      }
      else {
        sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x)+1))
        title <- paste0(y, "\nn=", length(x),"\nlog2(", n, "+1)")
      }
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(output_dir, n, "/heatmap_", y, ".pdf"), main=title, cluster_cols = F)
    }
  }
}

mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,9:(9+length(conditions))], n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_norm[,9:(8+nrow(samples))], n = normalization) )


# write file with TPM / RPM values for all samples, averaged over replicates ####
write.table(x=cts_summary_norm_mean, file=paste0(args[1], normalization, strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# function that adds annotations to a df using its row.names as Geneid
annotate_df <- function(df) { 
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}
# median of ratios (from DESeq2, median of ratios to geometric mean). Write tables and heatmaps ####
cts_MoR <- annotate_df( as.data.frame(counts(dds, normalized=T)) )
cts_MoR_mean <- average_replicates(cts_MoR)
write.table(x=cts_MoR, file=paste0(args[1], "ESF", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=cts_MoR_mean, file=paste0(args[1], "ESF", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_MoR[,9:(8+nrow(samples))], n = "MoR") )
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_MoR_mean[,9:(9+length(conditions))], n = "MoR") )
# exploring differences between median of ratios normalization by DESeq and RPM / TPM ####
if (FALSE==TRUE) { # just protecting this code so it's not executed when i run the script
  test <- cts_summary[,sample_columns]
  par(pty = "s")
  plot(x=colSums(test) / colSums(test)[6], y=sizeFactors(dds), xlim=c(0.5,1.7), ylim=c(0.5,1.7), xlab=normalization)
  abline(a = 0, b=1)
  text(x=colSums(test) / colSums(test)[6], sizeFactors(dds), labels=names(sizeFactors(dds)), cex= 0.5, pos=1)
}
