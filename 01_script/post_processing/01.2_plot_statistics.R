#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(ggbeeswarm)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_folder <- args[1]

# Create output folder
output_folder <- file.path(input_folder, "01_QC")
dir.create(output_folder, showWarnings = FALSE)

# Read statistics
read_stats <- as_tibble(read.delim(file.path(input_folder, "01_QC/all_read_counts_summarized.txt"), header=TRUE, sep="\t"))

# Plot percent_unique_UMIs
pdf(file.path(output_folder, "percent_unique_UMIs.pdf"), width = 15, height = 6)
ggplot(read_stats %>%
         mutate(percent_unique_UMIs = (Umi_Collapsed / Trimmed) * 100) %>%
         pivot_longer(cols = c(percent_unique_UMIs), names_to = "read", values_to = "percentage"),
       aes(x = Sample, y = percentage)) +
  geom_col(position = "identity") +
  labs(x = "Sample", y = "percentage of unique UMIs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Plot millions of unique_UMIs
pdf(file.path(output_folder, "number_unique_UMI_reads.pdf"), width = 15, height = 6)
ggplot(read_stats %>%
         pivot_longer(cols = c(Umi_Collapsed), values_to = "number_of_reads", names_to = "read") %>%
         mutate(
           number_of_reads = number_of_reads / 1e6,
           highlight = if_else(number_of_reads < 2.5, "Below 2.5M", "Above or Equal 2.5M")
         ),
       aes(x = Sample, y = number_of_reads, fill = highlight)) +
  geom_col(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("Below 2.5M" = "red", "Above or Equal 2.5M" = "green4")) +
  labs(x = "Sample", y = "Number of UMI-collapsed reads (in millions)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
dev.off()

# Plot all read types in a single histogram
pdf(file.path(output_folder, "all_read_types_histogram.pdf"), width = 15, height = 6)
ggplot(read_stats %>%
         pivot_longer(cols = c(Raw, Trimmed, Umi_Collapsed), values_to = "number_of_reads", names_to = "read") %>%
         mutate(number_of_reads = number_of_reads / 1e6),
       aes(x = Sample, y = number_of_reads, fill = read)) +
  geom_col(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("Raw" = "pink", "Trimmed" = "lightblue", "Umi_Collapsed" = "green4")) +
  labs(x = "Sample", y = "Number of Reads (in millions)", fill = "Read Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Read alignment statistics
alignment_stats <- as_tibble(read.delim(file.path(input_folder, "01_QC/pipeline_statistics.tsv"), header=TRUE, sep="\t"))

alignment_long_percent <- alignment_stats %>%
  filter(str_detect(parameter, '%')) %>%
  pivot_longer(
    cols = -parameter,
    names_to = "sample",
    values_to = "value"
  )

# Plot alignment statistics
pdf(file.path(output_folder, "alignment_statistics_percent.pdf"), width = 5, height = 4)
set.seed(1)
ggplot(alignment_long_percent, aes(x = parameter, y = value)) +
  geom_quasirandom(alpha=0.5) +
  scale_x_discrete(labels = ~ str_wrap(
    gsub(":",": ",gsub("_", " ", as.character(.x)))
    , width = 10)
  ) +  labs(title = "Distribution of alignment statistics", x = "", y = "Percentage") +
  ylim(0,100) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())
dev.off()

alignment_long_reads <- alignment_stats %>%
  filter(!str_detect(parameter, '%')) %>%
  pivot_longer(
    cols = -parameter,
    names_to = "sample",
    values_to = "value"
  )

# Plot alignment statistics
pdf(file.path(output_folder, "alignment_statistics_reads.pdf"), width = 5, height = 4)
set.seed(1)
ggplot(alignment_long_reads, aes(x = parameter, y = value)) +
  geom_quasirandom(alpha=0.5) +
  scale_x_discrete(labels = ~ str_wrap(
    gsub(":",": ",gsub("_", " ", as.character(.x)))
    , width = 10)
    ) +
  labs(title = "Distribution of alignment statistics", x = "", y = "Reads") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())
dev.off()

