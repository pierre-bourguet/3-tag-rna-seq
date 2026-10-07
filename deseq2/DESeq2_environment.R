# environment_DESeq2.R

# libraries

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

# Paths to additional scripts
source("R_functions/DEG_heatmap.R")
source("R_functions/graphical_parameters.R")
source("R_functions/plot_heatmap_plate_batch_effect.R")

# remove error messages from complexheatmap
ht_opt$message = FALSE

# Annotation files ####################################################################################################

# TAIR10 ATTEs
TEs <- as_tibble(read.delim("../../02_genome_files/02_input/TAIR10_Transposable_Elements.txt"
                            , header = TRUE, sep = "\t", quote = "", comment.char = "", col.names = c("Geneid", "sense", "start", "end", "family", "superfamily")))

# TAIR10 PCGs + pseudogenes
# features <- as_tibble(read.delim("../../02_genome_files/03_output/TAIR10_GFF_PCG_pseudogene.tsv", header = TRUE, sep = "\t", quote = "", comment.char = ""))
features <- as_tibble(read.delim("../../02_genome_files/03_output/TAIR10_GFF_PCG.tsv", header = TRUE, sep = "\t", quote = "", comment.char = ""))


# gene annotations
# this contains two types of information, depending on the annotation type:
# for non-TEs: the functional information from Araport11, and the predicted subcellular compartment of the gene product for protein coding genes, using the ".1" isoform, which is canonical in most but not all cases
# for TEs: the family and superfamily of the TE

Araport11_annotations <- as_tibble(read.delim("../../02_genome_files/03_output/Araport11_gene_ATTE_annotations.tsv"
                                              , header = FALSE, sep = "\t", quote = "", comment.char = "", col.names = c("chr", "start", "end", "Geneid", "strand", "type", "comment_1", "comment_2")))

# create duplicate annotations to match with antisense DEGs
Araport11_annotations_AS <- Araport11_annotations %>%
  mutate(Geneid = paste0(Geneid, "_AS"))
Araport11_annotations <- rbind (Araport11_annotations, Araport11_annotations_AS)


# TAIR10 ATTEs intersecting with TAIR10 PCGs ####################################################################################################

# import data
column_names <- c("chr", "start", "end", "Geneid", "type", "strand")
TE_PCG_intersect_sense <- as_tibble(read.delim("../../02_genome_files/03_output/TAIR10_ATTE_PCG_intersect_same_orientation.tsv"
                            , header = F, sep = "\t", quote = "", comment.char = "", col.names = c(paste0(column_names, "_TE"), paste0(column_names, "_PCG"), "overlap")))
TE_PCG_intersect_antisense <- as_tibble(read.delim("../../02_genome_files/03_output/TAIR10_ATTE_PCG_intersect_opposite_orientation.tsv"
                                               , header = F, sep = "\t", quote = "", comment.char = "", col.names = c(paste0(column_names, "_TE"), paste0(column_names, "_PCG"), "overlap")))

# convert overlap to fraction overlap
TE_PCG_intersect_sense <- TE_PCG_intersect_sense %>%
  mutate(overlap = overlap / (end_TE - start_TE))
TE_PCG_intersect_antisense <- TE_PCG_intersect_antisense %>%
  mutate(overlap = overlap / (end_TE - start_TE))

# add a "_AS" suffix to TEs in antisense orientation
TE_PCG_intersect_antisense <- TE_PCG_intersect_antisense %>%
  mutate(Geneid_TE = paste0(Geneid_TE, "_AS"))
