# function to plot a heatmap showing the row and columns of the input plate, to control for local batch effect
create_heatmap_with_annotations <- function(sampleDistMatrix, sample_wells) {
  
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  suppressPackageStartupMessages(require(RColorBrewer))
  suppressPackageStartupMessages(require(dplyr))
  
  # Ensure row names in sampleDistMatrix match with sample names in sample_wells
  sample_wells <- sample_wells[match(rownames(sampleDistMatrix), sample_wells$sample), ]
  
  # Define color palettes for annotations
  row_colors <- colorRampPalette(rev(brewer.pal(9, "Paired")))(length(unique(sample_wells$row)))
  column_colors <- colorRampPalette(rev(brewer.pal(12, "Set3")))(length(unique(sample_wells$column)))
  
  # Create color mapping functions
  row_col_fun <- structure(row_colors, names = unique(sample_wells$row))
  column_col_fun <- structure(column_colors, names = as.character(unique(sample_wells$column)))
  
  # Create row and column annotations
  row_annotation <- rowAnnotation(
    row = sample_wells$row,
    col = list(row = row_col_fun)
  )
  
  column_annotation <- HeatmapAnnotation(
    column = as.numeric(sample_wells$column),
    col = list(column = column_col_fun)
  )
  
  # Define the color function for the heatmap
  breaks <- seq(min(sampleDistMatrix), max(sampleDistMatrix), length.out = 255)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  col_fun <- colorRamp2(breaks, colors)
  
  # Create the heatmap
  Heatmap(
    sampleDistMatrix,
    name = "distance",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    left_annotation = row_annotation,
    top_annotation = column_annotation,
    width=15, height=15
  )
}

# used for analysis of 3' tagseq data, but is suitable for any other data
#
# Usage:
# sample_well should have these columns: sample well row column
# sampleDistMatrix should be a distance matrix with row.names containing sample names