#!/usr/bin/env Rscript
# QC plots from all_read_counts_summarized.txt and pipeline_statistics.tsv (written by pipeline_stats.py).
# Usage: plot_qc.R <qc_dir> [min_reads_warn]
#   min_reads_warn: UMI-collapsed read number below which a sample is flagged (default 2.5e6)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)
qc_dir <- args[1]
min_reads_warn <- if (length(args) >= 2) as.numeric(args[2]) else 2.5e6

read_stats <- read_tsv(file.path(qc_dir, "all_read_counts_summarized.txt"), show_col_types = FALSE) %>%
  mutate(Sample = factor(Sample, levels = Sample))
n_samples <- nrow(read_stats)
sample_axis <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot_width <- max(6, 2 + 0.12 * n_samples)

# read numbers after each step
p <- read_stats %>%
  pivot_longer(c(Raw, Trimmed, Umi_Collapsed), names_to = "step", values_to = "reads") %>%
  mutate(step = factor(step, levels = c("Raw", "Trimmed", "Umi_Collapsed"))) %>%
  ggplot(aes(Sample, reads / 1e6, fill = step)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  scale_fill_manual(values = c(Raw = "black", Trimmed = "grey60", Umi_Collapsed = "green4")) +
  labs(x = NULL, y = "Reads (millions)", fill = NULL, title = "Reads after each step", subtitle = paste0("n = ", n_samples, " samples")) +
  theme_classic() + sample_axis
ggsave(file.path(qc_dir, "all_read_counts_summarized.pdf"), p, width = plot_width, height = 5)

# fraction of trimmed reads with a unique UMI-sequence combination
p <- read_stats %>%
  mutate(percent_unique = 100 * Umi_Collapsed / Trimmed) %>%
  ggplot(aes(Sample, percent_unique)) +
  geom_col() +
  labs(x = NULL, y = "Unique UMIs (% of trimmed reads)", title = "UMI complexity", subtitle = paste0("n = ", n_samples, " samples")) +
  theme_classic() + sample_axis
ggsave(file.path(qc_dir, "percent_unique_UMIs.pdf"), p, width = plot_width, height = 5)

# UMI-collapsed reads, low-depth samples flagged
p <- read_stats %>%
  mutate(low = Umi_Collapsed < min_reads_warn) %>%
  ggplot(aes(Sample, Umi_Collapsed / 1e6, fill = low)) +
  geom_col() +
  geom_hline(yintercept = min_reads_warn / 1e6, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "green4"), guide = "none") +
  labs(x = NULL, y = "UMI-collapsed reads (millions)", title = "Library depth after UMI collapsing",
       subtitle = paste0("n = ", n_samples, " samples; red: < ", min_reads_warn / 1e6, "M reads (",
                         sum(read_stats$Umi_Collapsed < min_reads_warn), " samples)")) +
  theme_classic() + sample_axis
ggsave(file.path(qc_dir, "number_unique_UMI_reads.pdf"), p, width = plot_width, height = 5)

# mapping statistics, one point per sample
stats <- read_tsv(file.path(qc_dir, "pipeline_statistics.tsv"), show_col_types = FALSE) %>%
  pivot_longer(-parameter, names_to = "sample", values_to = "value") %>%
  mutate(value = suppressWarnings(as.numeric(value)))
wrap_labels <- function(x) str_wrap(gsub(":", ": ", gsub("_", " ", x)), width = 12)
plot_stats <- function(df, ylab, file) {
  set.seed(1)
  p <- ggplot(df, aes(parameter, value)) +
    geom_quasirandom(shape = 16, stroke = 0, size = 1.5, alpha = 0.6) +
    scale_x_discrete(labels = wrap_labels) +
    labs(x = NULL, y = ylab, title = "Mapping statistics", subtitle = paste0("n = ", n_distinct(df$sample), " samples")) +
    theme_classic()
  ggsave(file.path(qc_dir, file), p, width = 1.2 + 0.9 * n_distinct(df$parameter), height = 4)
}
plot_stats(filter(stats, str_detect(parameter, "%")), "Percentage", "alignment_statistics_percent.pdf")
plot_stats(filter(stats, !str_detect(parameter, "%")), "Reads", "alignment_statistics_reads.pdf")
