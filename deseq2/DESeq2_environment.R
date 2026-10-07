# Libraries, helper functions and annotation tables for DESeq2_tagseq.R

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

source(file.path(script_dir, "R_functions/DEG_heatmap.R"))
source(file.path(script_dir, "R_functions/graphical_parameters.R"))
source(file.path(script_dir, "R_functions/plot_heatmap_plate_batch_effect.R"))

ht_opt$message <- FALSE

# reference_manifest.tsv: key<TAB>value, '#' comments
read_manifest <- function(path) {
  m <- read.delim(path, header = FALSE, comment.char = "#", quote = "", colClasses = "character", col.names = c("key", "value"))
  setNames(m$value, m$key)
}

# Annotation tables listed in the reference manifest:
#   TEs: TE set of the reference (TAIR10 ATTEs, or TE genes with the family of the TE they derive from)
#   features: TAIR10 protein-coding genes
#   Araport11_annotations: functional annotation (genes) or family/superfamily (TEs), duplicated with _AS suffixes
#   TE_PCG_intersect_sense / _antisense: TEs overlapping a protein-coding gene on the same / opposite strand;
#     overlap as a fraction of TE length; antisense TE IDs carry the _AS suffix
load_annotations <- function(manifest) {
  rd <- function(key, ...) as_tibble(read.delim(manifest[[key]], sep = "\t", quote = "", comment.char = "", ...))

  TEs <<- rd("deseq2_te", header = TRUE, col.names = c("Geneid", "sense", "start", "end", "family", "superfamily"))
  features <<- rd("deseq2_pcg", header = TRUE)

  ann <- rd("deseq2_annotations", header = FALSE, col.names = c("chr", "start", "end", "Geneid", "strand", "type", "comment_1", "comment_2"))
  Araport11_annotations <<- bind_rows(ann, mutate(ann, Geneid = paste0(Geneid, "_AS")))

  column_names <- c("chr", "start", "end", "Geneid", "type", "strand")
  cols <- c(paste0(column_names, "_TE"), paste0(column_names, "_PCG"), "overlap")
  TE_PCG_intersect_sense <<- rd("deseq2_te_pcg_same", header = FALSE, col.names = cols) %>%
    mutate(overlap = overlap / (end_TE - start_TE))
  TE_PCG_intersect_antisense <<- rd("deseq2_te_pcg_opposite", header = FALSE, col.names = cols) %>%
    mutate(overlap = overlap / (end_TE - start_TE), Geneid_TE = paste0(Geneid_TE, "_AS"))
}

# Sample sheet: pipeline CSV (sample,fastq[,well,condition,replicate]) or legacy headerless TSV (fastq<TAB>sample).
# condition defaults to the sample name without its _R<n> suffix; well defaults to the fastq name prefix (A1_...).
read_samples <- function(path) {
  first <- readLines(path, n = 1)
  if (grepl("^sample,|,sample,|,sample$", first) && grepl("fastq", first)) {
    s <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else {
    s <- read_tsv(path, col_names = c("fastq", "sample"), show_col_types = FALSE, col_types = cols(.default = "c"))
  }
  if (!"condition" %in% names(s)) s$condition <- str_remove(s$sample, "_R\\d+$")
  s$condition[is.na(s$condition) | s$condition == ""] <- str_remove(s$sample, "_R\\d+$")[is.na(s$condition) | s$condition == ""]
  if (!"replicate" %in% names(s)) s$replicate <- str_extract(s$sample, "R\\d+$")
  if (!"well" %in% names(s)) s$well <- sub(".*/([^/]+)_(.*)\\..*", "\\1", s$fastq)
  s %>% transmute(full_name = sample, condition, replicate, well)
}
