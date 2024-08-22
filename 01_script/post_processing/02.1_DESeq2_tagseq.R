#!/usr/bin/env Rscript

# import environment and arguments and prepare output directory ####

# Source environment file for paths, constants and features
# modify the environment file to analyze different sets of features

print("loading environment")
source("02.2_DESeq2_environment.R")

print("R: importing arguments")

# these are arguments for troubleshooting, they are over-ridden below by user-provided inputs when the whole script is executed
args <- c(
  "../../04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min50bp/02_counts/",
  "empty_R1,ddm1_a_long_b_2_R1,F2_ddm1_R1,F2_ddm1_a_2_R1,ddm1_ab_2_R1,ddm1_a_1_R3,a_long_b_R1,b_2_R1,ab_2_R2,ab_1_R1,a_2_R3",
  "empty_R1,ddm1_a_long_b_2_R1,F2_ddm1_a_2_R1,ddm1_ab_2_R1,ddm1_a_1_R3,a_long_b_R1,b_2_R1,ab_2_R2,ab_1_R1,a_2_R3",
  "mom1,F2_WT,F2_a_2",
  "../../03_sample_lists/sample_list_tagseq_01_cdca7.tsv",
  "Col_0"
)

args <- c(
  "../../04_output/tagseq_03_cdca7_complementation_AtRTD3_ATTE_150bp_3M_min_50bp/02_counts/",
  "WT_R6,cdca7_ab_R1,cdca7_ab_R6,cdca7_ab_dCter_R4",
  "none",
  "../../03_sample_lists/sample_list_tagseq_03_cdca7.tsv",
  "WT"
)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]
outliers <- strsplit(args[2], ",")[[1]]
outlier_patterns <- strsplit(args[3], ",")[[1]]
sample_info_file <- args[4]
reference_condition <- args[5]

# Define the arguments for printing
args_list <- list(
  "count directory" = base_dir,
  "outliers" = outliers,
  "outlier patterns" = outlier_patterns,
  "sample info file" = sample_info_file,
  "reference condition" = reference_condition
)

# save current directory
current_dir <- getwd()
# set working directory
setwd(base_dir)

# Create output directory for plots and tables
output_dir <- "../06_DESeq2/"
dir.create(paste0(output_dir, "default_plots/"), showWarnings = F, recursive = T)

# prepare output log file
sink(paste0("../06_DESeq2/DESeq2_log.txt"))

# Print each argument with an additional newline between each
for (arg_name in names(args_list)) {
  cat(paste(arg_name, ": ", args_list[[arg_name]], "\n\n", sep = ""))
}

sink()

#
# import counts ####
print("importing counts")

# Import and format counts
counts_sense <- as_tibble(read.delim("star_counts.tsv", header = TRUE, sep = "\t"))
counts_AS <- as_tibble(read.delim("star_counts_AS.tsv", header = TRUE, sep = "\t"))

# Add suffix to Geneids and remove "_AS" suffix from column names
counts_AS$Geneid <- paste0(counts_AS$Geneid, "_AS")
names(counts_AS) <- gsub("_AS", "", names(counts_AS))

# Remove outliers & reorder columns alphabetically
counts_sense <- counts_sense %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns))))) %>%
  dplyr::select(order(names(.)))

counts_AS <- counts_AS %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns))))) %>%
  dplyr::select(order(names(.)))

# Verify that column names match for sense and antisense count dataframes
if (!all(names(counts_sense) == names(counts_AS))) {
  stop("Column names do not match between sense and antisense count dataframes")
}

# Combine the two dataframes
counts <- rbind(counts_sense, counts_AS)

# Merge the counts at the gene level and remove transcript fusions
counts_merged <- counts %>%
  dplyr::mutate(Geneid = gsub("\\.\\d+", "", Geneid)) %>%  # Remove the .1, .2, etc. at the end of Geneid
  dplyr::group_by(Geneid) %>%
  dplyr::summarise(across(everything(), sum)) %>%  # Sum the counts of all isoforms for each gene
  dplyr::filter(!str_detect(Geneid, "-")) %>%  # Remove transcript fusions
  dplyr::filter(str_detect(Geneid, "(^AT[1-5]|mTurq_3xcMyc|Venus)")) %>%  # Only keep geneids from chromosome 1 to 5, or transgenic genes containing mTurq_3xcMyc or Venus
  dplyr::filter(Geneid %in% TEs$Geneid | Geneid %in% paste0(TEs$Geneid, "_AS") | Geneid %in% features$Geneid | Geneid %in% paste0(features$Geneid, "_AS") | Geneid %in% c("mTurq_3xcMyc", "mTurq_3xcMyc_AS", "pAlli_Venus", "pAlli_Venus_AS"))  # Only keep geneids from TEs and selected features

# import RPM values, merge transcript variants & replicates, write merged files ####
print("importing RPM values and merging replicates")

RPM_sense <- as_tibble(read.delim("normalized_counts/no_filter_no_transcript_merge/star_RPM.tsv", header = TRUE, sep = "\t"))
RPM_AS <- as_tibble(read.delim("normalized_counts/no_filter_no_transcript_merge/star_RPM_AS.tsv", header = TRUE, sep = "\t"))

# Add suffix to geneids and remove "_AS" suffix from column names
RPM_AS$Geneid <- paste0(RPM_AS$Geneid, "_AS")
names(RPM_AS) <- gsub("_AS", "", names(RPM_AS))

# Remove outliers & reorder columns alphabetically
RPM_sense <- RPM_sense %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns))))) %>%
  dplyr::select(order(names(.)))

RPM_AS <- RPM_AS %>%
  dplyr::select(all_of(names(.) %>%
                         setdiff(outliers) %>%
                         discard(~ any(str_detect(.x, outlier_patterns))))) %>%
  dplyr::select(order(names(.)))

# Verify that column names match for sense and antisense count dataframes
if (!all(names(RPM_sense) == names(RPM_AS))) {
  stop("Column names do not match between sense and antisense count dataframes")
}

# Combine the two dataframes
RPM <- rbind(RPM_sense, RPM_AS)

# Merge the counts at the gene level, remove transcript fusions and genes not on chromosome 1 to 5
RPM_merged <- RPM %>%
  dplyr::mutate(Geneid = gsub("\\.\\d+", "", Geneid)) %>%  # Remove the .1, .2, etc. at the end of Geneid
  dplyr::group_by(Geneid) %>%
  dplyr::summarise(across(everything(), sum)) %>%  # Sum the RPM of all isoforms for each gene
  dplyr::filter(!str_detect(Geneid, "-")) %>%  # Remove transcript fusions
  dplyr::filter(str_detect(Geneid, "^AT[1-5]"))  # Only keep geneids from chromosome 1 to 5

# Average the replicates
process_RPM <- function(df) {
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
RPM_merged_avg <- process_RPM(RPM_merged)

# Write the RPM merged dataframe
write.table(RPM_merged, file = "normalized_counts/RPM.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(RPM_merged_avg, file = "normalized_counts/RPM_averaged.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

#
# import sample information ####
print("importing sample information")

# import sample info
samples <- read_tsv(paste0(current_dir, "/", sample_info_file), col_names = c("path", "full_name"), show_col_types = FALSE)[,2] %>%
  mutate(sample = str_extract(full_name, "R\\d+$")) %>%
  mutate(condition = str_remove(full_name, "_R\\d+$")) %>%
  # re order by alphabetical order
  dplyr::arrange(condition, sample, .locale = "en")

# Remove outliers
samples <- samples %>%
  dplyr::filter(
    !full_name %in% outliers &
      !str_detect(full_name, paste(outlier_patterns, collapse = "|"))
  )

sample_columns <- which(names(counts_merged) %in% paste(samples$condition, samples$sample, sep = "_"))  # Retrieve columns that contain sample names

#
# run DESeq2 ####
print("running DESeq2")

# find all conditions different from reference
conditions <- unique(samples$condition)[unique(samples$condition) != reference_condition]

# DESeq2 function
DESeq2_function <- function(x) {  # x should be cts_summary.tsv file (summary of raw counts)
  
  # round up in case you have fractional counts
  cts <- round(x[, sample_columns])
  row.names(cts) <- x$Geneid
  
  # Defining metadata
  coldata <- data.frame(
    condition = gsub("_R\\d+$", "", names(x)[sample_columns]),
    type = rep("single-strand", nrow(samples))
  )
  row.names(coldata) <- names(x)[sample_columns]
  if (!all(rownames(coldata) == colnames(cts))) {
    stop("the names in count_file and the sample table do not match")
  }
  dds <<- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ condition)
  
  # Adding meta data to the dataframe
  mcols(dds) <- DataFrame(mcols(dds))
  
  # Pre-filtering, here keeping only genes with at least 10 reads in at least 3 samples
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  
  # Setting the reference treatment
  dds$condition <- relevel(dds$condition, ref = reference_condition)
  
  # Differential analysis
  return(DESeq(dds))
  
}

dds <- DESeq2_function(counts_merged)

# Save the DESeq2 object
save(dds, file = "../06_DESeq2/DESeq2_object.RData")

#
# PCA using VST ####

print("plotting PCA using vst transformation")

# vst transformation
vsd <- vst(dds, blind=FALSE)

#### PCAs

# extract PCA data to filter samples easily
ntop_variable_features <- 1000 # number of most variable features for PCA
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE, ntop = ntop_variable_features)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Define a range of shapes to use
set.seed(12)
shape_list <- sample(15:25, length(unique(pcaData$condition)), replace = TRUE) # Using a set of distinct shapes

## all samples
PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"),
       y=paste0("PC2: ", percentVar[2], "% variance"),
       title=paste0("PCA with ", ntop_variable_features, " most variable features")) +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal()

# save the plot
ggsave(plot = PCA_plot, paste0(output_dir, "default_plots/PCA.pdf"), width = 10, height = 6)

## barplot of PC1 values to find outliers if needed

# Create an ordered dataframe
data <- as_tibble(pcaData) %>%
  mutate(condition = as.factor(condition)) %>% # Make sure condition is a factor
  arrange(desc(PC1)) # Arrange by PC1 to get the order

# barplot
barplot_PC1 <- ggplot(data
       , aes(x = PC1, y = reorder(name, PC1), fill = condition)) +
  geom_bar(stat = "identity", orientation = "y", colour="black") +
  labs(title = paste0("Ordered PC1 Values (", ntop_variable_features, " most variable features)"), x = "PC1", y = "Condition") +
  scale_fill_manual(values=c(col_vibrant, col_high_contrast, col_bright[-7], col_muted)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) # Ensure y-axis labels are readable

ggsave(plot = barplot_PC1, paste0(output_dir, "default_plots/PCA_barplot_of_PC1.pdf"), width = 8, height = 10)

#
# heatmaps of euclidean distances between samples using VST ####

# import plate well position to control if plate position correlates with batch effects
sample_wells <- read_tsv(paste0(current_dir, "/", sample_info_file), col_names = FALSE, show_col_types = FALSE)

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

# sample distance matrix with all samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- names(vsd$sizeFactor)
colnames(sampleDistMatrix) <- NULL
# heatmap
pdf(paste0(output_dir, "/default_plots/euclidean_distance_plate_position_effect.pdf"), width=17, height=15)
create_heatmap_with_annotations(sampleDistMatrix, sample_wells)
dev.off()

#
# normalization with DESeq2 estimated size factors (ESF, median of ratios to geometric mean): write tables ####
print("normalizing counts with estimated size factors")

# median of ratios
ESF <- as.data.frame(counts(dds, normalized=T)) %>%
  mutate(Geneid = row.names(.)) %>%
  dplyr::select(Geneid, everything())

# average replicates
ESF_avg <- ESF %>%
  pivot_longer(cols = -Geneid, names_to = "sample", values_to = "expression") %>%
  mutate(condition = str_replace(sample, "_R\\d+$", "")) %>% # remove the _R + number suffix
  group_by(Geneid, condition) %>%
  summarise(average_expression = mean(expression), .groups = 'drop') %>% # average expression
  pivot_wider(names_from = condition, values_from = average_expression, names_sort = T)

# write table of library-size normalized counts
write.table(x=ESF, file=paste0("normalized_counts/ESF.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=ESF_avg, file=paste0("normalized_counts/ESF_averaged.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

#
# variance-stabilizing transformation (VST): write tables ####

print("export variance-stabilized counts")

vsd <- as_tibble(assay(vsd)) %>%
  mutate(Geneid = row.names(assay(dds))) %>%
  dplyr::select(Geneid, everything())

# average replicates
vsd_avg <- vsd %>%
  pivot_longer(cols = -Geneid, names_to = "sample", values_to = "expression") %>%
  mutate(condition = str_replace(sample, "_R\\d+$", "")) %>% # remove the _R + number suffix
  group_by(Geneid, condition) %>%
  summarise(average_expression = mean(expression), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = average_expression, names_sort = T)

# write tables
write.table(x=vsd, file=paste0("normalized_counts/vst.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=vsd_avg, file=paste0("normalized_counts/vst_averaged.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

#
# regularized-log transformation (rlog): write tables ####
print("rlog transformation")

rld <- rlog(dds, blind=FALSE)
rld <- as_tibble(assay(rld)) %>%
  mutate(Geneid = row.names(assay(dds))) %>%
  dplyr::select(Geneid, everything())

# average replicates
rld_avg <- rld %>%
  pivot_longer(cols = -Geneid, names_to = "sample", values_to = "expression") %>%
  mutate(condition = str_replace(sample, "_R\\d+$", "")) %>% # remove the _R + number suffix
  group_by(Geneid, condition) %>%
  summarise(average_expression = mean(expression), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = average_expression, names_sort = T)

# write tables
write.table(x=rld, file=paste0("normalized_counts/rlog.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=rld_avg, file=paste0("normalized_counts/rlog_averaged.tsv"), quote = F, sep="\t", row.names=F, col.names=T)

#
# functions to extract pairwise comparisons and DEGs ####

# function to output a list of all treatment vs control comparisons
f <- function(aa, bb) { # this creates a res_mutant dataframe comparing mutant vs control
  eval(substitute( a <- results(dds, contrast=c("condition", as.character(b), reference_condition))
                   , list(a = aa, b = bb))) # this second argument provides an environment for substitute, defining variables used in previous line
}


# functions to extract up or down DEGs for annotations of interest
up <- function(res, annotations) {
  
  # res is the input dataframe
  # annotations matches a feature dataframe to subset
  
  # retrieve geneids
  geneids <- c(eval(parse(text=paste0(annotations, "$Geneid"))), paste0(eval(parse(text=paste0(annotations, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  
  # return geneids for significant DEGs
  return( row.names(res[row.names(res) %in% geneids & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.1,]) )
}
down <- function(res, annotations) {
  
  # res is the input dataframe
  # annotations matches a feature dataframe to subset
  
  # retrieve geneids
  geneids <- c(eval(parse(text=paste0(annotations, "$Geneid"))), paste0(eval(parse(text=paste0(annotations, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  
  # return geneids for significant DEGs
  return( row.names(res[row.names(res) %in% geneids & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.1,]) )
  
}

# these functions filter out TEs:
  # up or down in sense overlapping a PCG in the same orientation
  # up or down in antisense overlapping a PCG in the opposite orientation

up_TEs_intersect_filter <- function(res, TE_df, feature_df) {
  
  # this function filters out TEs:
    # up in sense overlapping a PCG in the same orientation
    # up in antisense overlapping a PCG in the opposite orientation
  
  # this function DOES NOT filter out cases where PCG antisense transcription would cause artifactual TE upregulation in sense or antisense
  
  # retrieve geneids
  feature_geneids <- c(eval(parse(text=paste0(feature_df, "$Geneid"))), paste0(eval(parse(text=paste0(feature_df, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  TE_geneids <- c(eval(parse(text=paste0(TE_df, "$Geneid"))), paste0(eval(parse(text=paste0(TE_df, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  
  # return geneids for significant DEGs
  up_TEs <- row.names(res[row.names(res) %in% TE_geneids & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.1,])
  up_features <- row.names(res[row.names(res) %in% feature_geneids & res$log2FoldChange >=1 & !is.na(res$padj) & res$padj < 0.1,])
  
  # filter out:
  # - TEs up in sense that have 1 bp overlap or more with a PCG in the same orientation
  # - TEs up in antisense that have 1 bp overlap or more with a PCG in the inverse orientation
  up_TEs <- up_TEs[!up_TEs %in% c(
    TE_PCG_intersect_sense$Geneid_TE[TE_PCG_intersect_sense$overlap > 0],
    TE_PCG_intersect_antisense$Geneid_TE[TE_PCG_intersect_sense$overlap > 0]
  )
  ]
  
  # with at least 1 bp overlap, identify up TEs that intersect up PCG in the same and inverse orientation
  # now this is redundant with the previous filter, has no effect
  sense_upTE_overlap <- TE_PCG_intersect_sense %>%
    filter(Geneid_PCG %in% up_features & Geneid_TE %in% up_TEs) %>%
    pull(Geneid_TE)
  antisense_upTE_overlap <- TE_PCG_intersect_antisense %>%
    filter(Geneid_PCG %in% up_features & Geneid_TE %in% up_TEs) %>%
    pull(Geneid_TE)
  
  # filter out these TEs
  up_TEs <- up_TEs[!up_TEs %in% c(sense_upTE_overlap, antisense_upTE_overlap)]
  
  return(up_TEs)
  
}
down_TEs_intersect_filter <- function(res, TE_df, feature_df) {
  
  # this function filters out TEs:
  # down in sense overlapping a PCG in the same orientation
  # down in antisense overlapping a PCG in the opposite orientation
  
  # this function DOES NOT filter out cases where PCG antisense transcription would cause artifactual TE downregulation in sense or antisense
  
  # retrieve geneids
  feature_geneids <- c(eval(parse(text=paste0(feature_df, "$Geneid"))), paste0(eval(parse(text=paste0(feature_df, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  TE_geneids <- c(eval(parse(text=paste0(TE_df, "$Geneid"))), paste0(eval(parse(text=paste0(TE_df, "$Geneid"))), "_AS")) # this is to include both sense and antisense quantifications
  
  # return geneids for significant DEGs
  down_TEs <- row.names(res[row.names(res) %in% TE_geneids & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.1,])
  down_features <- row.names(res[row.names(res) %in% feature_geneids & res$log2FoldChange <=-1 & !is.na(res$padj) & res$padj < 0.1,])
  
  # filter out:
  # - TEs down in sense that have 1 bp overlap or more with a PCG in the same orientation
  # - TEs down in antisense that have 1 bp overlap or more with a PCG in the inverse orientation
  down_TEs <- down_TEs[!down_TEs %in% c(
    TE_PCG_intersect_sense$Geneid_TE[TE_PCG_intersect_sense$overlap > 0],
    TE_PCG_intersect_antisense$Geneid_TE[TE_PCG_intersect_sense$overlap > 0]
  )
  ]
  
  # with at least 1 bp overlap, identify down TEs that intersect down PCG in the same and inverse orientation
  # now this is redundant with the previous filter, has no effect
  sense_downTE_overlap <- TE_PCG_intersect_sense %>%
    filter(Geneid_PCG %in% down_features & Geneid_TE %in% down_TEs) %>%
    pull(Geneid_TE)
  antisense_downTE_overlap <- TE_PCG_intersect_antisense %>%
    filter(Geneid_PCG %in% down_features & Geneid_TE %in% down_TEs) %>%
    pull(Geneid_TE)
  
  # filter out these TEs
  down_TEs <- down_TEs[!down_TEs %in% c(sense_downTE_overlap, antisense_downTE_overlap)]
  
  return(down_TEs)
}

#
# DEGs found in batch mode (including all conditions): export tables and heatmaps #### 
print("DEGs found in batch mode: export tables and heatmaps")

# list of all treatment vs control comparisons
all_res <- Map(f, paste0("res_", conditions), as.list(conditions)) 

## find up and down annotations including all conditions
DEGs <- list(
  upTEs = unique(unlist(lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features"))),
  upfeatures = unique(unlist(lapply(FUN = up, X = all_res, annotations="features"))),
  downTEs = unique(unlist(lapply(FUN = down_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features"))),
  downfeatures = unique(unlist(lapply(FUN = down, X = all_res, annotations="features")))
)

## write DEG tables with normalized counts
DEG_write_tables <- function(x, y) {
  
  # Create output directory
  output_dir <- paste0("../06_DESeq2/batch_DEGs/", y)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper function to write tables
  write_table <- function(data, suffix) {
    df <- data %>%
      filter(Geneid %in% x) %>%
      left_join(Araport11_annotations, by = "Geneid")
    write.table(df, file = paste0(output_dir, "/batch_", y, suffix), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  # Write all required tables
  write_table(RPM_merged, "_RPM.tsv")
  write_table(RPM_merged_avg, "_RPM_averaged.tsv")
  write_table(ESF, "_ESF.tsv")
  write_table(ESF_avg, "_ESF_averaged.tsv")
  write_table(vsd, "_VST.tsv")
  write_table(vsd_avg, "_VST_averaged.tsv")
  write_table(rld, "_rlog.tsv")
  write_table(rld_avg, "_rlog_averaged.tsv")
}

# Apply the function to the list of DEGs
mapply(FUN = DEG_write_tables, x = DEGs, y = names(DEGs))


## write number of batch DEGs
write.table(x=t(as.data.frame(lapply(DEGs, FUN=length))), file=paste0(output_dir, "batch_DEGs/number_of_batch_DEGs.tsv"), quote = F, sep="\t", row.names=T, col.names=F)

## write heatmaps for batch DEGs

# Helper function to generate heatmaps
generate_heatmaps <- function(DEG_list, suffix, data, select_cols, method) {
  print(paste("DEGs found in batch mode:", method))
  
  mapply(
    FUN = DEG_heatmap,
    x = DEG_list,
    y = paste0(names(DEG_list), suffix),
    output_dir = paste0("../06_DESeq2/batch_DEGs/", names(DEG_list)),
    MoreArgs = list(z = data %>% dplyr::select(Geneid, one_of(select_cols)), n = method)
  )
}

# RPM
generate_heatmaps(DEGs, "_RPM.pdf", RPM_merged, samples$full_name, "RPM")
generate_heatmaps(DEGs, "_RPM_averaged.pdf", RPM_merged_avg, c(reference_condition, conditions), "RPM")

# ESF
generate_heatmaps(DEGs, "_ESF.pdf", ESF, samples$full_name, "ESF")
generate_heatmaps(DEGs, "_ESF_averaged.pdf", ESF_avg, c(reference_condition, conditions), "ESF")

# VST
generate_heatmaps(DEGs, "_VST.pdf", vsd, samples$full_name, "VST")
generate_heatmaps(DEGs, "_VST_averaged.pdf", vsd_avg, c(reference_condition, conditions), "VST")

# rlog
generate_heatmaps(DEGs, "_rlog.pdf", rld, samples$full_name, "rlog")
generate_heatmaps(DEGs, "_rlog_averaged.pdf", rld_avg, c(reference_condition, conditions), "rlog")

#
# all pairwise comparisons of treated vs control: export tables and heatmaps ####
print("all pairwise comparisons: export tables and heatmaps")

## export the entire DESeq2 table for each pairwise comparison

export_all_pairwise <- function(x, y) {
  
  # x is the input dataframe
  # y is the sample name
  
  # create output directory
  output_dir_pairwise <- paste0(output_dir, "pairwise_comparisons/", y, "_vs_", reference_condition)
  dir.create(output_dir_pairwise, showWarnings = F, recursive = T)
  
  # prepare the dataframe
  x <- as.data.frame(x) %>%
    mutate(Geneid = row.names(x)) %>%  # add Geneid column
    dplyr::select(Geneid, everything()) %>% # Geneid as the first column
    left_join(Araport11_annotations, by="Geneid") # add annotations
  
  # export
  write.table(x, file=paste0(output_dir_pairwise, "/", y, "_vs_", reference_condition, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
}

mapply(x = all_res
       , y = gsub("res_", "", names(all_res))
       , FUN = export_all_pairwise)

## number of DEGs for each pairwise comparison

# create a dataframe with the number of DEGs for each pairwise comparison
DEG_nb <- t(data.frame(
  up_TEs = unlist(lapply(FUN = length, X = lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features")))
  , up_features = unlist(lapply(FUN = length, X = lapply(FUN = up, X = all_res, annotations="features")))
  , down_TEs = unlist(lapply(FUN = length, X = lapply(FUN = down_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features")))
  , down_features = unlist(lapply(FUN = length, X = lapply(FUN = down, X = all_res, annotations="features")))
  )) %>%
  as.data.frame() %>%
  # rename columns to remove "res_" prefix
  setNames(gsub("res_", "", names(all_res))) %>%
  mutate(DEG = row.names(.))

# format for plotting  
DEG_nb_long <- DEG_nb %>%
  pivot_longer(cols = -DEG, names_to = "treatment", values_to = "number_of_DEGs") %>%
  # order the treatments as in the sample sheet
  mutate(treatment = factor(treatment, levels = unique(samples$condition)))

# plot the data as barplots using ggplot
barplot_DEGs <- ggplot(DEG_nb_long, aes(x = treatment, y = number_of_DEGs, fill = DEG)) +
geom_bar(stat = "identity") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Number of DEGs for each pairwise comparison", x = "treatment", y = "Number of DEGs") +
theme(
  legend.position = "none"
      , panel.grid.major.x = element_blank() ) +
facet_wrap(~DEG, scales = "free_y", ncol=1)

# save the plot
plot_width <- 2 + 0.15 * length(unique(samples$condition))
ggsave(plot = barplot_DEGs, filename = paste0(output_dir, "default_plots/number_of_DEGs.pdf"), width = plot_width, height = 8)

## export up and down TEs / genes

# function to export tables of normalized counts and heatmaps for a given set of DEGs
export_DEGs_pairwise <- function(x, name, DEG_list, DEG_type) {
  
  # Create output directory
  output_dir_pairwise <- file.path(output_dir, "pairwise_comparisons", paste0(name, "_vs_", reference_condition))
  dir.create(output_dir_pairwise, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir_pairwise, "all_replicates"), showWarnings = FALSE, recursive = TRUE)
  
  # Prepare the dataframe with Geneid and annotations
  z <- as.data.frame(x) %>%
    mutate(Geneid = rownames(x)) %>%
    select(Geneid, everything()) %>%
    filter(Geneid %in% DEG_list) %>%
    left_join(Araport11_annotations, by = "Geneid")
  
  # Helper function to export tables
  export_table <- function(data, suffix, directory) {
    write.table(data, file = file.path(directory, paste0(name, "_vs_", reference_condition, "_", DEG_type, suffix)), 
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  # Helper function to generate heatmaps
  generate_heatmap <- function(data, file_suffix, normalization, directory) {
    DEG_heatmap(
      x = DEG_list,
      y = paste0(name, "_vs_", reference_condition, "_", DEG_type, file_suffix),
      z = data, 
      n = normalization, 
      output_dir = directory
    )
  }
  
  # Export normalized counts and generate heatmaps for averaged data
  export_table(z %>% left_join(ESF_avg, by = "Geneid"), "_normalized_counts.tsv", output_dir_pairwise)
  export_table(z %>% left_join(vsd_avg, by = "Geneid"), "_VST.tsv", output_dir_pairwise)
  export_table(z %>% left_join(rld_avg, by = "Geneid"), "_rlog.tsv", output_dir_pairwise)
  
  generate_heatmap(ESF_avg, "_normalized_counts_heatmap.pdf", "ESF", output_dir_pairwise)
  generate_heatmap(vsd_avg, "_VST_heatmap.pdf", "VST", output_dir_pairwise)
  generate_heatmap(rld_avg, "_rlog_heatmap.pdf", "rlog", output_dir_pairwise)
  
  # Export normalized counts and generate heatmaps for all replicates
  export_table(z %>% left_join(ESF, by = "Geneid"), "_normalized_counts_all_replicates.tsv", file.path(output_dir_pairwise, "all_replicates"))
  export_table(z %>% left_join(vsd, by = "Geneid"), "_VST_all_replicates.tsv", file.path(output_dir_pairwise, "all_replicates"))
  export_table(z %>% left_join(rld, by = "Geneid"), "_rlog_all_replicates.tsv", file.path(output_dir_pairwise, "all_replicates"))
  
  generate_heatmap(ESF, "_normalized_counts_all_replicates_heatmap.pdf", "ESF", file.path(output_dir_pairwise, "all_replicates"))
  generate_heatmap(vsd, "_VST_all_replicates_heatmap.pdf", "VST", file.path(output_dir_pairwise, "all_replicates"))
  generate_heatmap(rld, "_rlog_all_replicates_heatmap.pdf", "rlog", file.path(output_dir_pairwise, "all_replicates"))
}

# apply the function to all pairwise comparisons

# up TEs
mapply(FUN = export_DEGs_pairwise
       , x = all_res # list of all pairwise comparisons
       , name = gsub("res_", "", names(all_res))
       , DEG_list = lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features") #
       , MoreArgs = list(DEG_type = "up_TEs"))

# down TEs
mapply(FUN = export_DEGs_pairwise
       , x = all_res # list of all pairwise comparisons
       , name = gsub("res_", "", names(all_res))
       , DEG_list = lapply(FUN = down_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features") #
       , MoreArgs = list(DEG_type = "down_TEs"))

# up features
mapply(FUN = export_DEGs_pairwise
       , x = all_res # list of all pairwise comparisons
       , name = gsub("res_", "", names(all_res))
       , DEG_list = lapply(FUN = up, X = all_res, annotations="features") #
       , MoreArgs = list(DEG_type = "up_genes"))

# down features
mapply(FUN = export_DEGs_pairwise
       , x = all_res # list of all pairwise comparisons
       , name = gsub("res_", "", names(all_res))
       , DEG_list = lapply(FUN = down, X = all_res, annotations="features") #
       , MoreArgs = list(DEG_type = "down_genes"))

#

# save the environment and write a R script file for downstream analysis ####
print("saving environment")
save.image(file = paste0(output_dir, "DESeq2_environment.RData"))

# Define the content of the new R script
script_content <- '# Load the DESeq2 environment
load("../06_DESeq2/DESeq2_environment.RData")
# your analysis ####
'
# Write the content to the file
file_path <- "../07_analysis/your_analysis.R"
dir.create(dirname(file_path), recursive = TRUE, showWarnings = F)
writeLines(script_content, file_path)

#
