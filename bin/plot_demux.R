#!/usr/bin/env Rscript
# Reads per well after demultiplexing (demultiplex_summary.tsv from make_samplesheet.py)
# Usage: plot_demux.R demultiplex_summary.tsv

suppressPackageStartupMessages(library(tidyverse))

d <- read_tsv(commandArgs(trailingOnly = TRUE)[1], show_col_types = FALSE)
unmatched <- d %>% filter(well == "no_barcode_match")
wells <- d %>%
  filter(well != "no_barcode_match") %>%
  mutate(label = paste0(well, " ", sample),
         label = factor(label, levels = label[order(reads)]),
         assigned = sample != "unassigned")

p <- ggplot(wells, aes(reads / 1e6, label, fill = assigned)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "grey40", `FALSE` = "darkorange"), labels = c(`TRUE` = "in well sheet", `FALSE` = "unassigned")) +
  labs(x = "Reads (millions)", y = NULL, fill = NULL, title = "Reads per well",
       subtitle = paste0("n = ", nrow(wells), " wells; no barcode match: ", round(unmatched$percent_of_total, 1), "% of reads")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6))
ggsave("demultiplex_summary.pdf", p, width = 7, height = 2 + 0.12 * nrow(wells))
