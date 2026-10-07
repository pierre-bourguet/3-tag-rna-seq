# colors ####

suppressPackageStartupMessages(library("scico"))
suppressPackageStartupMessages(library("circlize"))
col_muted <- c("#CC6677", "#88CCEE", "#DDCC77", "#117733", "#332288", "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD") # source is https://personal.sron.nl/~pault/
col_vibrant <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB')
col_high_contrast <- c("#FFFFFF", '#004488', '#DDAA33', '#BB5566', '#000000')
col_bright <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
col_rlog <- colorRamp2(seq(from=1, to=10, by=1), scico(n=10, direction=-1, palette="lajolla"))
many_colors <- c(col_vibrant, col_high_contrast[-1], col_bright[-7], col_muted)


# Function to generate similar colors for a given color
generate_similar_colors_x3 <- function(color) {
  similar_colors <- c(color, adjustcolor(color, 0.75), adjustcolor(color, 1.25))
  return(similar_colors)
}
generate_similar_colors_x2 <- function(color) {
  similar_colors <- c(color, adjustcolor(color, 0.75))
  return(similar_colors)
}

# Generate the new vector of colors
col_muted_3_replicates <- c()
for (i in 1:length(col_muted)) {
  new_group <- generate_similar_colors_x3(col_muted[i])
  col_muted_3_replicates <- c(col_muted_3_replicates, new_group)
}
col_muted_2_replicates <- c()
for (i in 1:length(col_muted)) {
  new_group <- generate_similar_colors_x2(col_muted[i])
  col_muted_2_replicates <- c(col_muted_2_replicates, new_group)
}

#
# graphical parameters ####
pt_0.5_to_mm <- 0.176389 # use this in mm to get 0.5 pt
mm_to_inches <- 0.0393701 # multiply mm by this to get inches

## .pt <- 72.27 / 25.4
## .stroke <- 96 / 25.4

## line width exact 1 pt 
ggplot_line_width_1pt <- 1 / (72.27 / 25.4) / (72.27 / 96)

## stroke width exact 1 pt 
ggplot_stroke_width_1pt <- 1 / (96 / 25.4) / (72.27 / 96) * 2




# ggplots themes ####

theme_horizontal_nature <- list(
  theme(
    panel.grid = element_blank()
    , axis.line.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.line.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , line = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    #, rect = element_rect(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")  # this gives a black rectangle around the plot
    , axis.text.x.top = element_blank()
    , axis.ticks.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , text = element_text(size = 6, family = "Arial") # change font size of all text
    , axis.text.x = element_text(size = 6, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    #, legend.title = element_text(size = 6, family = "Arial") # change font size of legend title
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent", color = NA)
    , legend.background = element_blank()
    , legend.key.size = unit(0.25, "cm")
    , legend.position = "inside"
    , legend.position.inside = c(0.75, 0.75)
    , legend.title = element_blank()
    # strip options are about labels that come with facet_wrap or facet_grid
    , strip.background = element_blank() # to remove rectangles around titles
    , strip.text = element_text(size = 6, family = "Arial")  )
)


theme_vertical_nature <- list(
  theme(
    panel.grid = element_blank()
    , panel.border = element_blank()
    , axis.line.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.line.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , line = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    #, rect = element_rect(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D") # this gives a black rectangle around the plot
    , axis.text.x.top = element_blank()
    , axis.ticks.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , text = element_text(size = 6, family = "Arial") # change font size of all text
    , axis.text.y = element_text(size = 6, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    #, legend.title = element_text(size = 6, family = "Arial") # change font size of legend title
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent", color = NA)
    , legend.background = element_blank()
    , legend.key.size = unit(0.25, "cm")
    , legend.position = "inside"
    , legend.position.inside = c(0.75, 0.75)
    , legend.title = element_blank()
    # strip options are about labels that come with facet_wrap or facet_grid
    , strip.background = element_blank() # to remove rectangles around titles
    , strip.text = element_text(size = 6, family = "Arial")    
  )
)

theme_scatterplot_nature <- list(
  theme(
    panel.grid = element_blank()
    , axis.line.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.line.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.x = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.ticks.y = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , line = element_line(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , rect = element_rect(linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D")
    , axis.text.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , text = element_text(size = 6, family = "Arial") # change font size of all text
    , axis.text = element_text(size = 6, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent", color = NA)
    , legend.background = element_blank()
    , legend.key.size = unit(0.25, "cm")
    , legend.position = "inside"
    , legend.position.inside = c(0.75, 0.75)
    , legend.title = element_blank()
    # strip options are about labels that come with facet_wrap or facet_grid
    , strip.background = element_blank() # to remove rectangles around titles
    , strip.text = element_text(size = 6, family = "Arial")
  )
)
#
