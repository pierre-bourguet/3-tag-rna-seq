# import DESeq2 environment
load("../06_DESeq2/DESeq2_environment.RData")
# quality controls: heatmaps of euclidean distances between samples ####

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

