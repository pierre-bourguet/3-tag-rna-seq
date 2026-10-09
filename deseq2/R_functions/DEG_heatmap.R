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

      # heatmap body of fixed size (0.25 in per column, 5 in high); the page is sized to the drawn heatmap
      ht <- Heatmap(
        sampleDistMatrix,
        name = n,
        col = colors,
        cluster_rows = as.dendrogram(hclust(dist(sampleDistMatrix))),  # clustered once for both draws below
        row_dend_reorder = TRUE,
        cluster_columns = FALSE,
        column_title = title,
        column_title_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        width = unit(max(0.25 * ncol(sampleDistMatrix), 1.5), "in"),
        height = unit(5, "in")
      )

      pdf(NULL)
      ht_drawn <- draw(ht, heatmap_legend_side = "right")
      plot_width <- convertWidth(ComplexHeatmap:::width(ht_drawn), "in", valueOnly = TRUE)
      plot_height <- convertHeight(ComplexHeatmap:::height(ht_drawn), "in", valueOnly = TRUE)
      # the title can be wider than the heatmap
      title_width <- max(strwidth(strsplit(title, "\n")[[1]], units = "inches", cex = 10 / par("ps")))
      invisible(dev.off())

      pdf(file_path, width = max(plot_width, title_width) + 0.5, height = plot_height + 0.3)
      draw(ht, heatmap_legend_side = "right")
      dev.off()
    }
  }
}
