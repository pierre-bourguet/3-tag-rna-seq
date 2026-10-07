#!/usr/bin/env Rscript
# DESeq2 analysis of 3' tag-seq counts (STAR -> salmon counts written by the pipeline in 02_counts/).
#
# Usage: Rscript DESeq2_tagseq.R <count_dir> <outliers> <outlier_patterns> <sample_sheet> <reference_condition> [--key=value ...]
#   count_dir            <run>/02_counts (star_counts.tsv, star_counts_AS.tsv, normalized_counts/no_filter_no_transcript_merge/)
#   outliers             comma-separated sample names to exclude, or "none"
#   outlier_patterns     comma-separated patterns; samples whose name contains one are excluded, or "none"
#   sample_sheet         samplesheet of the run (CSV sample,fastq[,well,condition,replicate] or legacy TSV fastq<TAB>sample)
#   reference_condition  condition every other condition is compared to
# Options:
#   --outdir=<dir>       default <count_dir>/../06_DESeq2; when set, normalized tables go to <outdir>/normalized_counts
#   --manifest=<file>    reference_manifest.tsv; default <count_dir>/../pipeline_info/reference_manifest.tsv, else the
#                        shared ATTE reference (runs made before the manifest existed)
#   --lfc=1 --padj=0.1   DEG thresholds (|log2FC| >= lfc and padj < padj)
#   --min_count=10 --min_samples=3   pre-filter: keep features with >= min_count reads in >= min_samples samples

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))

all_args <- commandArgs(trailingOnly = TRUE)
opt_args <- grep("^--[A-Za-z_]+=", all_args, value = TRUE)
args <- all_args[!all_args %in% opt_args]
if (length(args) != 5) stop("expected 5 positional arguments, see the header of DESeq2_tagseq.R")
opts <- setNames(sub("^--[A-Za-z_]+=", "", opt_args), sub("^--([A-Za-z_]+)=.*", "\\1", opt_args))
opt <- function(key, default) if (key %in% names(opts)) opts[[key]] else default

as_list <- function(x) if (identical(x, "none") || x == "") character(0) else strsplit(x, ",")[[1]]
base_dir <- normalizePath(args[1])
outliers <- as_list(args[2])
outlier_patterns <- as_list(args[3])
sample_info_file <- normalizePath(args[4])
reference_condition <- args[5]
lfc_threshold <- as.numeric(opt("lfc", 1))
padj_threshold <- as.numeric(opt("padj", 0.1))
min_count <- as.numeric(opt("min_count", 10))
min_samples <- as.numeric(opt("min_samples", 3))

print("loading environment")
source(file.path(script_dir, "DESeq2_environment.R"))

run_manifest <- file.path(base_dir, "..", "pipeline_info", "reference_manifest.tsv")
default_manifest <- file.path(script_dir, "..", "..", "..", "resources", "genomes", "a_thaliana", "AtRTD3", "processed", "tagseq_ATTE",
                              "reference_manifest.tsv")
manifest_file <- normalizePath(opt("manifest", if (file.exists(run_manifest)) run_manifest else default_manifest))
manifest <- read_manifest(manifest_file)
load_annotations(manifest)
transgene_ids <- if (!is.null(manifest[["unique_count_ids"]]) && nzchar(manifest[["unique_count_ids"]])) {
  strsplit(manifest[["unique_count_ids"]], ",")[[1]]
} else {
  c("mTurq_3xcMyc", "pAlli_Venus")
}

setwd(base_dir)
output_dir <- paste0(normalizePath(opt("outdir", file.path(base_dir, "..", "06_DESeq2")), mustWork = FALSE), "/")
norm_dir <- if ("outdir" %in% names(opts)) file.path(output_dir, "normalized_counts") else "normalized_counts"
dir.create(paste0(output_dir, "default_plots/"), showWarnings = FALSE, recursive = TRUE)
dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(c(
  paste0("count directory: ", base_dir),
  paste0("outliers: ", paste(outliers, collapse = ",")),
  paste0("outlier patterns: ", paste(outlier_patterns, collapse = ",")),
  paste0("sample info file: ", sample_info_file),
  paste0("reference condition: ", reference_condition),
  paste0("reference manifest: ", manifest_file, " (", manifest[["te_annotation"]], ")"),
  paste0("DEG thresholds: |log2FC| >= ", lfc_threshold, ", padj < ", padj_threshold),
  paste0("pre-filter: >= ", min_count, " reads in >= ", min_samples, " samples")
), paste0(output_dir, "DESeq2_log.txt"))

is_outlier <- function(x) {
  if (length(outlier_patterns) == 0) return(x %in% outliers)
  x %in% outliers | str_detect(x, paste(outlier_patterns, collapse = "|"))
}

# sample information ####
print("importing sample information")
samples <- read_samples(sample_info_file) %>%
  dplyr::arrange(condition, replicate, .locale = "en") %>%
  dplyr::filter(!is_outlier(full_name))
condition_of <- setNames(samples$condition, samples$full_name)
to_condition <- function(sample) ifelse(sample %in% names(condition_of), condition_of[sample], str_remove(sample, "_R\\d+$"))

# import counts ####
print("importing counts")

# sense and antisense tables, outliers removed, columns in alphabetical order; antisense features get an _AS suffix
read_strands <- function(sense_file, antisense_file) {
  sense <- as_tibble(read.delim(sense_file, header = TRUE, sep = "\t", check.names = FALSE))
  AS <- as_tibble(read.delim(antisense_file, header = TRUE, sep = "\t", check.names = FALSE))
  AS$Geneid <- paste0(AS$Geneid, "_AS")
  names(AS) <- sub("_AS$", "", names(AS))
  keep <- function(df) df %>% dplyr::select(Geneid, all_of(sort(setdiff(names(df), "Geneid")[!is_outlier(setdiff(names(df), "Geneid"))])))
  sense <- keep(sense)
  AS <- keep(AS)
  if (!all(names(sense) == names(AS))) stop("Column names do not match between sense and antisense tables")
  rbind(sense, AS)
}

# isoforms summed per gene, transcript fusions removed, chromosomes 1-5 and transgenes only
merge_isoforms <- function(df) {
  df %>%
    dplyr::mutate(Geneid = gsub("\\.\\d+", "", Geneid)) %>%
    dplyr::group_by(Geneid) %>%
    dplyr::summarise(across(everything(), sum)) %>%
    dplyr::filter(!str_detect(Geneid, "-"))
}

counts <- read_strands("star_counts.tsv", "star_counts_AS.tsv")
counts_merged <- merge_isoforms(counts) %>%
  dplyr::filter(str_detect(Geneid, "^AT[1-5]") | sub("_AS$", "", Geneid) %in% transgene_ids) %>%
  dplyr::filter(sub("_AS$", "", Geneid) %in% c(TEs$Geneid, features$Geneid, transgene_ids))

# RPM values (salmon TPM with --noLengthCorrection), merged per gene, averaged per condition ####
print("importing RPM values and merging replicates")
RPM <- read_strands("normalized_counts/no_filter_no_transcript_merge/star_RPM.tsv",
                    "normalized_counts/no_filter_no_transcript_merge/star_RPM_AS.tsv")
RPM_merged <- merge_isoforms(RPM) %>% dplyr::filter(str_detect(Geneid, "^AT[1-5]"))

average_conditions <- function(df, sort_columns = TRUE) {
  df %>%
    pivot_longer(cols = -Geneid, names_to = "sample", values_to = "value") %>%
    mutate(condition = to_condition(sample)) %>%
    group_by(Geneid, condition) %>%
    summarise(average = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = average, names_sort = sort_columns)
}

RPM_merged_avg <- average_conditions(RPM_merged, sort_columns = FALSE)
write.table(RPM_merged, file = file.path(norm_dir, "RPM.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(RPM_merged_avg, file = file.path(norm_dir, "RPM_averaged.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

missing <- setdiff(samples$full_name, names(counts_merged))
if (length(missing) > 0) stop("samples of the sample sheet missing from the count table: ", paste(missing, collapse = ", "))
sample_columns <- which(names(counts_merged) %in% samples$full_name)

# run DESeq2 ####
print("running DESeq2")
conditions <- unique(samples$condition)[unique(samples$condition) != reference_condition]

DESeq2_function <- function(x) {
  cts <- round(x[, sample_columns])
  row.names(cts) <- x$Geneid
  coldata <- data.frame(condition = condition_of[names(x)[sample_columns]], type = "single-strand")
  row.names(coldata) <- names(x)[sample_columns]
  if (!all(rownames(coldata) == colnames(cts))) stop("the names in count_file and the sample table do not match")
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  mcols(dds) <- DataFrame(mcols(dds))
  dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]
  dds$condition <- relevel(dds$condition, ref = reference_condition)
  DESeq(dds)
}

dds <- DESeq2_function(counts_merged)
save(dds, file = paste0(output_dir, "DESeq2_object.RData"))

# PCA using VST ####
print("plotting PCA using vst transformation")
vsd <- vst(dds, blind = FALSE)
ntop_variable_features <- 1000
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE, ntop = ntop_variable_features)
percentVar <- round(100 * attr(pcaData, "percentVar"))
set.seed(12)
shape_list <- sample(15:25, length(unique(pcaData$condition)), replace = TRUE)

PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = paste0("PCA with ", ntop_variable_features, " most variable features")) +
  scale_color_manual(values = many_colors) +
  scale_shape_manual(values = shape_list) +
  theme_minimal()
ggsave(plot = PCA_plot, paste0(output_dir, "default_plots/PCA.pdf"), width = 10, height = 6)

# PC1 per sample, to spot outliers
data <- as_tibble(pcaData) %>% mutate(condition = as.factor(condition)) %>% arrange(desc(PC1))
barplot_PC1 <- ggplot(data, aes(x = PC1, y = reorder(name, PC1), fill = condition)) +
  geom_bar(stat = "identity", orientation = "y", colour = "black") +
  labs(title = paste0("Ordered PC1 Values (", ntop_variable_features, " most variable features)"), x = "PC1", y = "Condition") +
  scale_fill_manual(values = c(col_vibrant, col_high_contrast, col_bright[-7], col_muted)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))
ggsave(plot = barplot_PC1, paste0(output_dir, "default_plots/PCA_barplot_of_PC1.pdf"), width = 8, height = 10)

# sample-to-sample euclidean distances (VST), annotated with plate row and column ####
sample_wells <- samples %>%
  transmute(sample = full_name, well) %>%
  mutate(row = substr(well, 1, 1), column = as.numeric(substr(well, 2, nchar(well))))
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
rownames(sampleDistMatrix) <- names(vsd$sizeFactor)
colnames(sampleDistMatrix) <- NULL
pdf(paste0(output_dir, "/default_plots/euclidean_distance_plate_position_effect.pdf"), width = 17, height = 15)
create_heatmap_with_annotations(sampleDistMatrix, sample_wells)
dev.off()

# normalized tables: size-factor normalized counts (ESF, median of ratios), VST, rlog; per sample and averaged ####
print("normalizing counts with estimated size factors")
ESF <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  mutate(Geneid = row.names(.)) %>%
  dplyr::select(Geneid, everything())
ESF_avg <- average_conditions(ESF)
write.table(x = ESF, file = file.path(norm_dir, "ESF.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = ESF_avg, file = file.path(norm_dir, "ESF_averaged.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("export variance-stabilized counts")
vsd <- as_tibble(assay(vsd)) %>%
  mutate(Geneid = row.names(assay(dds))) %>%
  dplyr::select(Geneid, everything())
vsd_avg <- average_conditions(vsd)
write.table(x = vsd, file = file.path(norm_dir, "vst.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = vsd_avg, file = file.path(norm_dir, "vst_averaged.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("rlog transformation")
rld <- rlog(dds, blind = FALSE)
rld <- as_tibble(assay(rld)) %>%
  mutate(Geneid = row.names(assay(dds))) %>%
  dplyr::select(Geneid, everything())
rld_avg <- average_conditions(rld)
write.table(x = rld, file = file.path(norm_dir, "rlog.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = rld_avg, file = file.path(norm_dir, "rlog_averaged.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# pairwise comparisons and DEG selection ####

# results of one condition vs the reference
f <- function(aa, bb) {
  eval(substitute(a <- results(dds, contrast = c("condition", as.character(b), reference_condition)), list(a = aa, b = bb)))
}

# sense and antisense IDs of an annotation set
with_AS <- function(ids) c(ids, paste0(ids, "_AS"))
is_up <- function(res) res$log2FoldChange >= lfc_threshold & !is.na(res$padj) & res$padj < padj_threshold
is_down <- function(res) res$log2FoldChange <= -lfc_threshold & !is.na(res$padj) & res$padj < padj_threshold

up <- function(res, annotations) {
  geneids <- with_AS(get(annotations)$Geneid)
  row.names(res[row.names(res) %in% geneids & is_up(res), ])
}
down <- function(res, annotations) {
  geneids <- with_AS(get(annotations)$Geneid)
  row.names(res[row.names(res) %in% geneids & is_down(res), ])
}

# Differential TEs, excluding those that overlap a protein-coding gene (>= 1 bp) in the orientation where gene
# transcription would be counted on the TE: sense TEs on a same-strand gene, antisense TEs on an opposite-strand
# gene. TEs whose signal could come from antisense transcription of a gene are not excluded.
TE_overlapping_PCG <- function() {
  c(TE_PCG_intersect_sense$Geneid_TE[TE_PCG_intersect_sense$overlap > 0],
    TE_PCG_intersect_antisense$Geneid_TE[TE_PCG_intersect_antisense$overlap > 0])
}
up_TEs_intersect_filter <- function(res, TE_df, feature_df) {
  ids <- row.names(res[row.names(res) %in% with_AS(get(TE_df)$Geneid) & is_up(res), ])
  ids[!ids %in% TE_overlapping_PCG()]
}
down_TEs_intersect_filter <- function(res, TE_df, feature_df) {
  ids <- row.names(res[row.names(res) %in% with_AS(get(TE_df)$Geneid) & is_down(res), ])
  ids[!ids %in% TE_overlapping_PCG()]
}

# DEGs in any comparison ("batch" DEGs): tables and heatmaps ####
print("DEGs found in batch mode: export tables and heatmaps")
all_res <- Map(f, paste0("res_", conditions), as.list(conditions))

DEGs <- list(
  upTEs = unique(unlist(lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df = "TEs", feature_df = "features"))),
  upfeatures = unique(unlist(lapply(FUN = up, X = all_res, annotations = "features"))),
  downTEs = unique(unlist(lapply(FUN = down_TEs_intersect_filter, X = all_res, TE_df = "TEs", feature_df = "features"))),
  downfeatures = unique(unlist(lapply(FUN = down, X = all_res, annotations = "features")))
)

DEG_write_tables <- function(x, y) {
  dir <- paste0(output_dir, "batch_DEGs/", y)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  write_table <- function(data, suffix) {
    df <- data %>% filter(Geneid %in% x) %>% left_join(Araport11_annotations, by = "Geneid")
    write.table(df, file = paste0(dir, "/batch_", y, suffix), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  write_table(RPM_merged, "_RPM.tsv")
  write_table(RPM_merged_avg, "_RPM_averaged.tsv")
  write_table(ESF, "_ESF.tsv")
  write_table(ESF_avg, "_ESF_averaged.tsv")
  write_table(vsd, "_VST.tsv")
  write_table(vsd_avg, "_VST_averaged.tsv")
  write_table(rld, "_rlog.tsv")
  write_table(rld_avg, "_rlog_averaged.tsv")
}
invisible(mapply(FUN = DEG_write_tables, x = DEGs, y = names(DEGs)))
write.table(x = t(as.data.frame(lapply(DEGs, FUN = length))), file = paste0(output_dir, "batch_DEGs/number_of_batch_DEGs.tsv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

generate_heatmaps <- function(DEG_list, suffix, data, select_cols, method) {
  print(paste("DEGs found in batch mode:", method))
  invisible(mapply(
    FUN = DEG_heatmap,
    x = DEG_list,
    y = paste0(names(DEG_list), suffix),
    output_dir = paste0(output_dir, "batch_DEGs/", names(DEG_list)),
    MoreArgs = list(z = data %>% dplyr::select(Geneid, one_of(select_cols)), n = method)
  ))
}
generate_heatmaps(DEGs, "_RPM.pdf", RPM_merged, samples$full_name, "RPM")
generate_heatmaps(DEGs, "_RPM_averaged.pdf", RPM_merged_avg, c(reference_condition, conditions), "RPM")
generate_heatmaps(DEGs, "_ESF.pdf", ESF, samples$full_name, "ESF")
generate_heatmaps(DEGs, "_ESF_averaged.pdf", ESF_avg, c(reference_condition, conditions), "ESF")
generate_heatmaps(DEGs, "_VST.pdf", vsd, samples$full_name, "VST")
generate_heatmaps(DEGs, "_VST_averaged.pdf", vsd_avg, c(reference_condition, conditions), "VST")
generate_heatmaps(DEGs, "_rlog.pdf", rld, samples$full_name, "rlog")
generate_heatmaps(DEGs, "_rlog_averaged.pdf", rld_avg, c(reference_condition, conditions), "rlog")

# each condition vs reference: full results, log2FC table, DEG numbers, DEG tables and heatmaps ####
print("all pairwise comparisons: export tables and heatmaps")

export_all_pairwise <- function(x, y) {
  dir <- paste0(output_dir, "pairwise_comparisons/", y, "_vs_", reference_condition)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  x <- as.data.frame(x) %>%
    mutate(Geneid = row.names(x)) %>%
    dplyr::select(Geneid, everything()) %>%
    left_join(Araport11_annotations, by = "Geneid")
  write.table(x, file = paste0(dir, "/", y, "_vs_", reference_condition, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
invisible(mapply(x = all_res, y = gsub("res_", "", names(all_res)), FUN = export_all_pairwise))

log2FC_df <- as_tibble(map_dfc(all_res, ~ as.data.frame(.x) %>% dplyr::select(log2FoldChange))) %>%
  rename_with(~ gsub("res_", "", names(all_res))) %>%
  mutate(Geneid = vsd$Geneid) %>%
  dplyr::select(Geneid, everything())
write.table(log2FC_df, file = file.path(norm_dir, "log2FC.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

DEG_nb <- t(data.frame(
  up_TEs = unlist(lapply(FUN = length, X = lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df = "TEs", feature_df = "features"))),
  up_features = unlist(lapply(FUN = length, X = lapply(FUN = up, X = all_res, annotations = "features"))),
  down_TEs = unlist(lapply(FUN = length, X = lapply(FUN = down_TEs_intersect_filter, X = all_res, TE_df = "TEs", feature_df = "features"))),
  down_features = unlist(lapply(FUN = length, X = lapply(FUN = down, X = all_res, annotations = "features")))
)) %>%
  as.data.frame() %>%
  setNames(gsub("res_", "", names(all_res))) %>%
  mutate(DEG = row.names(.))

DEG_nb_long <- DEG_nb %>%
  pivot_longer(cols = -DEG, names_to = "treatment", values_to = "number_of_DEGs") %>%
  mutate(treatment = factor(treatment, levels = unique(samples$condition)))

barplot_DEGs <- ggplot(DEG_nb_long, aes(x = treatment, y = number_of_DEGs, fill = DEG)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of DEGs for each\npairwise comparison", x = "treatment", y = "Number of DEGs") +
  theme(legend.position = "none", panel.grid.major.x = element_blank()) +
  facet_wrap(~DEG, scales = "free_y", ncol = 1)
plot_width <- 2 + 0.15 * length(unique(samples$condition))
ggsave(plot = barplot_DEGs, filename = paste0(output_dir, "default_plots/number_of_DEGs.pdf"), width = plot_width, height = 8)

export_DEGs_pairwise <- function(x, name, DEG_list, DEG_type) {
  dir <- file.path(output_dir, "pairwise_comparisons", paste0(name, "_vs_", reference_condition))
  dir.create(file.path(dir, "all_replicates"), showWarnings = FALSE, recursive = TRUE)
  z <- as.data.frame(x) %>%
    mutate(Geneid = rownames(x)) %>%
    dplyr::select(Geneid, everything()) %>%
    filter(Geneid %in% DEG_list) %>%
    left_join(Araport11_annotations, by = "Geneid")
  export_table <- function(data, suffix, directory) {
    write.table(data, file = file.path(directory, paste0(name, "_vs_", reference_condition, "_", DEG_type, suffix)),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  generate_heatmap <- function(data, file_suffix, normalization, directory) {
    DEG_heatmap(x = DEG_list, y = paste0(name, "_vs_", reference_condition, "_", DEG_type, file_suffix), z = data, n = normalization,
                output_dir = directory)
  }
  export_table(z %>% left_join(ESF_avg, by = "Geneid"), "_normalized_counts.tsv", dir)
  export_table(z %>% left_join(vsd_avg, by = "Geneid"), "_VST.tsv", dir)
  export_table(z %>% left_join(rld_avg, by = "Geneid"), "_rlog.tsv", dir)
  generate_heatmap(ESF_avg, "_normalized_counts_heatmap.pdf", "ESF", dir)
  generate_heatmap(vsd_avg, "_VST_heatmap.pdf", "VST", dir)
  generate_heatmap(rld_avg, "_rlog_heatmap.pdf", "rlog", dir)
  rep_dir <- file.path(dir, "all_replicates")
  export_table(z %>% left_join(ESF, by = "Geneid"), "_normalized_counts_all_replicates.tsv", rep_dir)
  export_table(z %>% left_join(vsd, by = "Geneid"), "_VST_all_replicates.tsv", rep_dir)
  export_table(z %>% left_join(rld, by = "Geneid"), "_rlog_all_replicates.tsv", rep_dir)
  generate_heatmap(ESF, "_normalized_counts_all_replicates_heatmap.pdf", "ESF", rep_dir)
  generate_heatmap(vsd, "_VST_all_replicates_heatmap.pdf", "VST", rep_dir)
  generate_heatmap(rld, "_rlog_all_replicates_heatmap.pdf", "rlog", rep_dir)
}

comparison_names <- gsub("res_", "", names(all_res))
for (type in list(list("up_TEs", up_TEs_intersect_filter, "TEs"), list("down_TEs", down_TEs_intersect_filter, "TEs"),
                  list("up_genes", up, "features"), list("down_genes", down, "features"))) {
  DEG_lists <- if (type[[3]] == "TEs") {
    lapply(FUN = type[[2]], X = all_res, TE_df = "TEs", feature_df = "features")
  } else {
    lapply(FUN = type[[2]], X = all_res, annotations = "features")
  }
  invisible(mapply(FUN = export_DEGs_pairwise, x = all_res, name = comparison_names, DEG_list = DEG_lists, MoreArgs = list(DEG_type = type[[1]])))
}

# environment for downstream analyses ####
print("saving environment")
save.image(file = paste0(output_dir, "DESeq2_environment.RData"))
analysis_script <- file.path(dirname(normalizePath(output_dir)), "07_analysis", "your_analysis.R")
dir.create(dirname(analysis_script), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(analysis_script)) {
  writeLines(c("# Load the DESeq2 environment", paste0('load("', output_dir, 'DESeq2_environment.RData")'), "# your analysis ####"), analysis_script)
}
