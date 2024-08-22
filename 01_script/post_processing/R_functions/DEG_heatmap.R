suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))

DEG_heatmap <- function(x, y, z, n, output_dir) {
  
  # x is a list of annotations
  # y is the name of the output file
  # z is the data frame with the expression values
  # n is the normalization method
  # output_dir is the output directory
  
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      if (n %in% c("rlog", "VST")) {
        sampleDistMatrix <- as.matrix(as.data.frame(z[z$Geneid %in% x,-which(names(z) == "Geneid")]))
        title <- paste0(y, "\nn=", length(x),"\n", n)
      }
      else {
        n <- paste0("log2 (", n,"+1)")
        sampleDistMatrix <- as.matrix(log2(z[z$Geneid %in% x,-which(names(z) == "Geneid") ]+1))
        title <- paste0(y, "\nn=", length(x),"\n", n)
      }
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      
      # output directory and filename
      ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
      file_path <- paste0(output_dir, "/", y)
      
      # Create the heatmap
      ht <- Heatmap(
        sampleDistMatrix,
        name = n,
        col = colors,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_title = title,
        width=15, height=15
      )
      
      # Save the heatmap to a PDF file
      plot_height <- min(  0.3 * ncol(sampleDistMatrix), 8)
      
      pdf(file_path, width = 0.3 * ncol(sampleDistMatrix), height = plot_height)
      draw(ht, heatmap_legend_side = "right")
      dev.off()
    }
  }
}
