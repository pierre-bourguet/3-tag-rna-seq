# colors ####

suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("scico"))
col_mC <- colorRampPalette(c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506'))(100)
col_mC_diff <- colorRamp2(c(-1,0,1), c("#2166AC", "#F7F7F7", "#B2182B")) 
col_mCHG_diff <- colorRamp2(c(-0.5,0,0.5), c("#2166AC", "#F7F7F7", "#B2182B")) 
col_mCHH_diff <- colorRamp2(c(-0.3,0,0.3), c("#2166AC", "#F7F7F7", "#B2182B")) 
col_zscore <- colorRamp2(c(-3,0,3), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
col_zscore <- colorRamp2(c(-3,0,3), c("cornflowerblue", "#FFFFF0", "brown1"))
col_log2_TPM <- colorRamp2(seq(from=6, to=11, by=1), scico(n=6, direction=-1, palette="lajolla"))
col_log10_TPM <- colorRamp2(seq(from=1, to=4, by=0.2), scico(n=16, direction=-1, palette="lajolla"))
col_ChIP_diff <- c('#762A83', '#9970AB', '#C2A5CF', '#E7D4E8', '#F7F7F7', '#D9F0D3', '#ACD39E', '#5AAE61', '#1B7837') # an alternative from https://personal.sron.nl/~pault/
col_ChIP_diff <- colorRamp2(c(-2,0,2), c("#762A83", "#F7F7F7", "#1B7837"))
col_deciles <- c('white', '#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506', '#888888') # source is "YlOrBr" on https://personal.sron.nl/~pault/
col_muted <- c("#CC6677", "#88CCEE", "#DDCC77", "#117733", "#332288", "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD") # source is https://personal.sron.nl/~pault/
col_vibrant <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB')
col_high_contrast <- c("#FFFFFF", '#004488', '#DDAA33', '#BB5566', '#000000')
col_bright <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
col_pale <- c('#BBCCEE', '#CCEEFF', '#CCDDAA', '#EEEEBB', '#FFCCCC', '#DDDDDD') # for text, source is https://personal.sron.nl/~pault/
col_dark <- c('#222255', '#225555', '#225522', '#666633', '#663333', '#555555') # for text background, source is https://personal.sron.nl/~pault/
col_sunset <- c('#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B', '#DD3D2D', '#A50026', '#FFFFFF')
colors_h2aw_scatterplot <- c("#AA4499", "grey20")
colors_h2aw_barplot <- c("black", "#AA4499")
colors_h2aw_boxplot <- c("grey70", "#AA4499")
col_TE_families <- colorRampPalette(c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506'))(100)
col_rlog <- colorRamp2(seq(from=1, to=10, by=1), scico(n=10, direction=-1, palette="lajolla"))
col_linear_pale_to_brown <- c('#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29', '#EC7014', '#CC4C02', '#993404', '#662506', '#888888')
col_TE_size <- colorRampPalette(c('#FEFBE9', '#FCF7D5', '#F5F3C1', '#EAF0B5', '#DDECBF', '#D0E7CA', '#C2E3D2', '#B5DDD8', '#A8D8DC', '#9BD2E1', '#8DCBE4', '#81C4E7', '#7BBCE7', '#7EB2E4', '#88A5DD', '#9398D2', '#9B8AC4', '#9D7DB2', '#9A709E', '#906388', '#805770', '#684957', '#46353A', '#999999'))(100)
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

# graphical functions from the discontinued colortools package https://github.com/gastonstat/colortools/tree/master ####

#' @title Pizza color wheel
#' 
#' @description
#' This function displays a color wheel with specified colors
#' 
#' @details
#' This function is based on the \code{\link{pie}} function
#' 
#' @param colors a vector with R color names of colors in hexadecimal notation
#' @param bg background color of the plot. Default \code{"gray95"}
#' @param border color of the border separating the pizza slices
#' @param init.angle integer value indicating the start angle (in degrees) for
#' the slices
#' @param cex numeric value indicating the character expansion of the labels
#' @param lty argument passed to \code{\link{polygon}} which draws each slice
#' @param labcol color for the labels (i.e. names of the colors)
#' @param \dots graphical parameters (\code{\link{par}}) can be given as
#' argument to \code{pizza}
#' @author Gaston Sanchez
#' @seealso \code{\link{wheel}}
#' @export
#' @examples
#'
#' # pizza color wheel for rainbow colors
#' pizza(rainbow(7))
#' 
#' # pizza color wheel for tomato (18 colors)
#' pizza(setColors("tomato", 18), bg = "gray20", cex = 0.7)
#'
pizza <-
  function(colors, bg = "gray95", border = NA, 
           init.angle = 105, cex = 0.8, lty = 1, labcol = NULL, ...)
  {
    n <- length(colors)
    x <- rep(1, n)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    # set colors
    labels = colors
    # 
    if (is.null(labcol))
    {
      if (mean(col2rgb(bg)) > 127)
        labcol = rep("black", n)
      if (mean(col2rgb(bg)) <= 127)
        labcol = rep("white", n)
    }
    # prepare plot window
    par(bg = bg)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
      xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    # get ready to plot
    border <- rep(border, length.out = nx)
    if (is.null(border[1]))
      border <- rep(bg, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(45, length.out = nx)
    radius = seq(1, 0, by=-1/n)[1:n]
    twopi <- -2 * pi
    t2xy <- function(t, rad) {
      t2p <- twopi * t + init.angle * pi/180
      list(x = rad * cos(t2p), y = rad * sin(t2p))
    }
    # plot colored segments
    for (i in 1L:nx)
    {
      n <- max(2, floor(200 * dx[i]))
      P <- t2xy(seq.int(x[i], x[i + 1], length.out = n), rad=radius[1])
      polygon(c(P$x, 0), c(P$y, 0), angle = angle[i], 
              border = border[i], col = colors[i], lty = lty[i])
      P <- t2xy(mean(x[i + 0:1]), rad=radius[1])
      lab <- labels[i]
      if (!is.na(lab) && nzchar(lab)) {
        adjs = 0.5
        if (P$x > 1e-08) adjs <- 0
        if (P$x < -1e-08) adjs <- 1
        lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y, col=labcol[i])
        text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
             adj = adjs, cex=cex, col=labcol[i], ...)
      }
    }
    invisible(NULL)
  }

#'@title col2HSV: converts a color to HSV in hexadecimal notation
#'
#'@description
#'col2HSV converts an R color (or a set of colors) into an HSV color model, and
#'then returns the color names in hexadeciaml notation
#'
#'@param color an R color name or a color in hexadecimal notation
#'@return A character vector with the color(s) name(s) in hexadecimal notation
#'@author Gaston Sanchez
#'@seealso \code{\link{wheel}}
#'@export
#'@examples
#'
#' # convert 'tomato'
#' col2HSV("tomato")
#'
col2HSV <-
  function(color)
  {
    # convert to RGB
    rgb_col = col2rgb(color)
    # convert to HSV
    hsv_col = rgb2hsv(rgb_col)
    if (length(color) == 1)
    {
      # get degree
      hue = hsv_col[1]
      sat = hsv_col[2]
      val = hsv_col[3]
      # get colors with hsv
      hex_col = hsv(hue, sat, val)
    }
    if (length(color) > 1)
    {
      hex_col = rep("", length(color))
      for (j in 1:length(color))
      {
        hex_col[j] = hsv(hsv_col[1,j], hsv_col[2,j], hsv_col[3,j])
      }
    }
    hex_col
  }

#'@title Set Colors for a color wheel
#'
#'@description
#'This function set a given number of colors to create a color wheel
#'
#'
#'@param color an R color name or a color in hexadecimal notation
#'@param num integer value indicating how many colors to be added to the wheel
#'@return A character vector with the given color and the set of colors to
#'create a wheel color
#'@author Gaston Sanchez
#'@seealso \code{\link{col2HSV}}
#'@export
#'@examples
#'
#' # create a color wheel based on 'tomato'
#' setColors("tomato", 12)
#' 
#' # set 7 colors for '#3D6DCC'
#' setColors("#3D6DCC", 7)
#'
setColors <-
  function(color, num)
  {
    # convert to RGB
    rgb_col = col2rgb(color)
    # convert to HSV
    hsv_col = rgb2hsv(rgb_col)[,1]
    # get degree
    hue = hsv_col[1]
    sat = hsv_col[2]
    val = hsv_col[3]
    cols = seq(hue, hue + 1, by=1/num)
    cols = cols[1:num]
    cols[cols > 1] <- cols[cols > 1] - 1
    # get colors with hsv
    colors = hsv(cols, sat, val)
    # transparency
    if (substr(color, 1, 1) == "#" && nchar(color) == 9)
    {
      alpha = substr(color, 8, 9)
      colors = paste(colors, alpha, sep="")
    }
    colors
  }

#' @title Adjacent or analogous colors
#' 
#' @description
#' Adjacent color schemes use colors that are next to each other on the color
#' wheel. These colors usually match well and create comfortable designs.
#' 
#' @details
#' The analogous colors are obtained following a color wheel with 12 colors,
#' each one spaced at 30 degrees from each other.
#' 
#' @aliases adjacent analogous
#' @param color an R color name or a color in hexadecimal notation
#' @param plot logical value indicating whether to plot a color wheel with the
#' generated scheme
#' @param bg background color of the plot. Used only when \code{plot=TRUE}
#' @param labcol color for the labels (i.e. names of the colors). Used only when
#' \code{plot=TRUE}
#' @param cex numeric value indicating the character expansion of the labels
#' @param title logical value indicating whether to display a title in the plot.
#' Used only when \code{plot=TRUE}
#' @return A character vector with the given color and the analogous colors in
#' hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{complementary}}, \code{\link{splitComp}},
#' \code{\link{triadic}}, \code{\link{tetradic}}, \code{\link{square}}
#' @export
#' @examples
#' # analogous colors of 'red'
#' adjacent("red", plot = FALSE)
#' 
#' # analogous colors of 'tomato' with default color wheel
#' analogous("tomato")
#' 
#' # analogous colors of '#606FEF' with darker background
#' adjacent("#606FEF", bg = "gray20")
#'
adjacent <-
  function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE)
  {	
    tmp_cols = setColors(color, 12)
    adja_colors <- tmp_cols[c(1,2,12)]
    
    # plot
    if (plot)
    {
      # labels color
      if (is.null(labcol)) 
      {
        lab_col = rep("", 12)
        if (mean(col2rgb(bg)) > 127)
        {
          lab_col[c(1, 2, 12)] <- "black"
          lab_col[c(3:11)] <- col2HSV(bg)
        } else {
          lab_col[c(1, 2, 12)] <- "white"
          lab_col[c(3:11)] <- col2HSV(bg)
        }
      } else {
        lab_col = rep(labcol, 12)
        if (mean(col2rgb(bg)) > 127)
        {
          lab_col[c(1, 2, 12)] <- labcol
          lab_col[c(3:11)] <- col2HSV(bg)
        } else {
          lab_col[c(1, 2, 12)] <- labcol
          lab_col[c(3:11)] <- col2HSV(bg)
        }
      }	
      # hide non-adjacent colors
      tmp_cols[c(3:11)] <- paste(substr(tmp_cols[c(3:11)],1,7), "0D", sep="")
      pizza(tmp_cols, labcol=lab_col, bg=bg, cex=cex)
      # title
      if (title)
        title(paste("Adjacent (analogous) colors of: ", tmp_cols[1]), 
              col.main=lab_col[1], cex.main=0.8)
    }
    # result
    adja_colors
  }

analogous <-
  function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE) 
  {
    adjacent(color, plot=plot, bg=bg, labcol=labcol, cex=cex, title=title)
  }

#'@title sequential HSV colors
#'
#'@description
#'This functions allows to get a sequence of colors in an HSV model with
#'optional pre-especified numbers for saturation, value, and alpha. It is a
#'very flexible function to play with different combinations of saturation,
#'value, and alpha.
#'
#'@details
#'The idea bechind this function is to explore a sequence of colors given some
#'fixed numbers of saturation, valur or alpha for an HSV color model. The
#'argument \code{what} will be taken to generate the sequence in the given
#'\code{percentage} increment steps. In addition, we can specify a number for
#'\code{s, v, alpha}. For example, if \code{what="value"}, we can fix the
#'saturation in \code{s=0.8}, obtaining a sequence of colors with different
#'values but with the same level of saturation.
#'
#'The argument \code{fun} allows to apply a transformation to the generated
#'sequence. By default \code{fun="linear"}, no transformation is applied. If
#'\code{fun="sqrt"}, the square root of the generated sequence will be taken.
#'If \code{fun="log"}, the logarithmic of the generated sequence will be taken.
#'
#'@param color an R color name or a color in hexadeciaml notation
#'@param percentage numeric value indicating the increment steps of the
#'sequence in percentage
#'@param what character string indicating what parameter to taki into account
#'to generate the sequence. Possible values are \code{"saturation"},
#'\code{"value"}, and \code{alpha}
#'@param s optional decimal value (between 0 and 1) to fix the color saturation
#'@param v optional decimal value (between 0 and 1) to fix the color value
#'@param alpha optional decimal value (between 0 and 1) to fix the color alpha
#'transparency
#'@param fun character string indicating the applied transformation to the
#'generated sequence. Possible values are \code{"linear"}, \code{"sqrt"}, and
#'\code{"log"}
#'@param plot logical value indicating whether to plot the sequence
#'@param verbose logical value indicating whether to return the color names of
#'the sequence
#'@author Gaston Sanchez
#'@seealso \code{\link{pizza}}
#'@export
#'@examples
#'
#' # sequence for 'orange'
#' sequential("orange")
#' 
#' # sequence for 'orange' with fun='sqrt' transformation
#' sequential("orange", fun = "sqrt")
#' 
#' # sequence for 'orange' with fun='log' transformation
#' sequential("orange", fun = "log")
#' 
#' # sequential sequence for value with fix saturation s=0.7 and fun='log'
#' sequential("orange", what = "value", s = 0.7, fun = "log")
#' 
#' # sequential sequence for saturation, with fix value s=0.8, alpha=0.5, percentage 10, and fun='log'
#' sequential("orange", 10, what = "value", s = 0.7, alpha = 0.5, fun = "log")
#'
sequential <-
  function(color, percentage=5, what="saturation", 
           s=NULL, v=NULL, alpha=NULL, fun="linear", plot=TRUE, verbose=TRUE)
  {	
    # convert to HSV
    col_hsv = rgb2hsv(col2rgb(color))[,1]
    # transparency
    if (is.null(alpha))
      alpha = 1
    if (substr(color, 1, 1) == "#" && nchar(color) == 9)
      alpha = substr(color, 8, 9)
    # get hue, saturation, and value
    hue = col_hsv[1]
    if (is.null(s)) s = col_hsv[2]
    if (is.null(v)) v = col_hsv[3]
    # sequence function
    getseq = switch(fun, 
                    linear = seq(0, 1, by=percentage/100),
                    sqrt = sqrt(seq(0, 1, by=percentage/100)),
                    log = log1p(seq(0, 1, by=percentage/100)),
                    log10 = log10(seq(0, 1, by=percentage/100))
    )
    # what type of sequence?
    if (what == "saturation") {
      sat = getseq
      fixed = paste("v=", round(v,2), " and alpha=", alpha, sep="")
      if (is.numeric(alpha))
        seq_col = hsv(hue, s=sat, v=v, alpha=alpha)
      if (is.character(alpha)) {
        seq_col = hsv(hue, s=sat, v=v)
        seq_col = paste(seq_col, alpha, sep="")
      }
    }
    if (what == "value") {
      val = getseq
      fixed = paste("s=", round(s,2), " and alpha=", alpha, sep="")
      if (is.numeric(alpha))
        seq_col = hsv(hue, s=s, v=val, alpha=alpha)
      if (is.character(alpha)) {
        seq_col = hsv(hue, s=s, v=val)
        seq_col = paste(seq_col, alpha, sep="")
      }
    }
    if (what == "alpha") {
      alpha = getseq
      fixed = paste("s=", round(s,2), " and v=", round(v,2), sep="")
      seq_col = hsv(hue, s=s, v=v, alpha=alpha)
    }
    # if plot TRUE
    if (plot)
    {
      n = length(seq(0, 1, by=percentage/100))
      fx = unlist(fixed)
      #dev.new()
      plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
      rect(0:(n-1)/n, 0, 1:n/n, 1, col=seq_col, border="lightgray")
      mtext(seq_col, side=1, at=0.5:(n)/n, cex=0.8, las=2)
      title(paste("Sequential colors based on ", what, "\n with fixed ", fx, sep=""),
            cex.main=0.9)
    }
    # result
    if (verbose)
      seq_col
  }

# my own function to generate similar colors from an input vector ####
# Assuming the analogous and sequential functions are already defined

generate_combined_colors_for_vector <- function(colors, sample_sizes) {
  # Check if the length of colors and sample_sizes are the same
  if (length(colors) != length(sample_sizes)) {
    stop("The length of colors vector must match the length of sample_sizes vector.")
  }
  
  # Initialize an empty list to store the final colors
  all_combined_colors <- list()
  
  # Loop through each input color
  for (i in 1:length(colors)) {
    color <- colors[i]
    sample_size <- sample_sizes[i]
    
    # Generate 3 analogous colors
    analogous_colors <- analogous(color, plot = FALSE)
    
    # Initialize an empty list to store the final colors for this color
    final_colors <- list()
    
    # Loop through each analogous color
    for (j in 1:length(analogous_colors)) {
      # Generate a sequence of 21 colors
      seq_colors <- sequential(analogous_colors[j], plot = FALSE)
      
      # Select colors at indices 3, 6, 9, 12, 15, 18, and 21
      #selected_indices <- seq(3, 21, 3)
      #selected_colors <- seq_colors[selected_indices]
      
      # Add the selected colors to the final list for this color
      #final_colors[[j]] <- selected_colors
      final_colors[[j]] <- seq_colors
    }
    
    # Flatten the list into a single vector
    combined_colors <- unlist(final_colors)
    
    # Sample the specified number of colors from the combined colors
    sampled_colors <- sample(combined_colors, sample_size)
    
    # Add the sampled colors to the main list
    all_combined_colors[[color]] <- sampled_colors
  }
  
  return(all_combined_colors)
}

# Example usage
colors <- c("#CC6677", "#CC6677", "#3357FF")
sample_sizes <- c(5, 7, 6)
combined_colors <- generate_combined_colors_for_vector(colors, sample_sizes)

print(combined_colors)

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
custom_theme <-  list(theme_bw(),
                      theme(text = element_text(size = 6), line = element_line(linewidth=0.176), rect = element_rect(linewidth=0.176), strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
                            , axis.text=element_text(size = 5)
                            , panel.grid.major = element_line(linewidth=0.0353, color="black"), panel.grid.minor = element_line(linewidth=0.0353), axis.ticks=element_line(linewidth=0.0353))
)

scatterplot_regression_theme <- list(geom_point(size=0.6, alpha=0.8, stroke=0), scale_color_manual(values=colors_h2aw_scatterplot)
                                     , coord_fixed(), geom_abline(intercept = 0, slope = 1, linewidth=0.1, linetype="dashed"), geom_smooth(alpha=0.8, linewidth=0, fill="grey80", method=lm), geom_line(stat="smooth", method=lm, size=0.2, alpha=0.9)
                                     , theme(legend.position = "none", axis.title = element_blank()))

scatterplot_regression_theme_lower_alpha <- list(geom_point(size=0.6, stroke=0, alpha=0.2, color="#AA4499"), scale_color_manual(values=colors_h2aw_scatterplot)
                                                 , coord_fixed(), geom_abline(intercept = 0, slope = 1, linewidth=0.1, linetype="dashed"), geom_smooth(alpha=0.8, linewidth=0, fill="grey80", method=lm), geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9, color="grey11")
                                                 , theme(legend.position = "none", axis.title = element_blank()))
theme_boxplot_TPM <- list(theme(axis.title.x=element_blank(), panel.grid.major.x = element_blank(), axis.ticks.x=element_blank()) )

theme_violinplot_chromatin_marks <-  list(
  theme_bw(),  scale_fill_manual(values=col_muted), geom_violin(aes(fill=cluster), lwd = 0.24), stat_summary(fun=median, geom="point", size=0.1, color="grey26", position = position_dodge(0.9)),
  theme(text = element_text(size = 6), line = element_line(linewidth=0.176), rect = element_rect(linewidth=0.176), legend.key.size=unit(0.4, "cm"), strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        , axis.text=element_text(size = 5), strip.text=element_text(size = 6), axis.title.x=element_blank(), axis.ticks.x=element_blank()
        , panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth=0.0353, color="black"), panel.grid.minor.y = element_line(linewidth=0.0353), panel.spacing.y=unit(1,"mm"), axis.ticks=element_line(linewidth=0.0353))
)

theme_boxplot_cluster_TPM <-
  list(theme_bw(),
       theme(text = element_text(size = 6),
             line = element_line(linewidth=0.176),
             rect = element_rect(linewidth=0.176),
             strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
             panel.grid.major = element_line(linewidth=0.0353, color="black"),
             panel.grid.minor = element_line(linewidth=0.0353),
             panel.grid.major.x = element_blank(),
             axis.text=element_text(size = 5),
             axis.ticks.y=element_line(linewidth=0.0353),
             axis.ticks.x=element_blank(),
             axis.title.x=element_blank())
  )

theme_horizontal_nature <- list(
  theme(
    panel.grid.minor = element_blank()
    , panel.grid.major.y = element_blank()
    , panel.grid.major.x = element_line(linewidth = ggplot_line_width_1pt/2, color="#EBEBEB")
    , panel.border = element_blank()
    , line = element_line(linewidth = ggplot_line_width_1pt/2)
    , rect = element_rect(linewidth = ggplot_line_width_1pt/2)
    , axis.text.x.top = element_blank()
    , axis.ticks.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , axis.ticks = element_blank()
    , text = element_text(size = 5, family = "Arial") # change font size of all text
    , axis.text = element_text(size = 5, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    , legend.title = element_text(size = 6, family = "Arial") # change font size of legend title
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent")
  )
)


theme_vertical_nature <- list(
  theme(
    panel.grid.minor = element_blank()
    , panel.grid.major.y = element_line(linewidth = ggplot_line_width_1pt/2)
    , panel.border = element_blank()
    , panel.grid.major.x = element_blank()
    , line = element_line(linewidth = ggplot_line_width_1pt/2)
    , rect = element_rect(linewidth = ggplot_line_width_1pt/2)
    , axis.text.x.top = element_blank()
    , axis.ticks.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , axis.ticks = element_blank()
    , text = element_text(size = 5, family = "Arial") # change font size of all text
    , axis.text = element_text(size = 5, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    , legend.title = element_text(size = 6, family = "Arial") # change font size of legend title
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent", color = NA)
  )
)

theme_scatterplot_nature <- list(
  theme(
    panel.grid.minor = element_blank()
    , panel.grid.major.y = element_line(linewidth = ggplot_line_width_1pt/2)
    , panel.grid.major.x = element_line(linewidth = ggplot_line_width_1pt/2)
    , line = element_line(linewidth = ggplot_line_width_1pt/2)
    , rect = element_rect(linewidth = ggplot_line_width_1pt/2)
    , axis.text.x.top = element_blank()
    , axis.ticks.x.top = element_blank()
    , axis.line.x.top = element_blank()
    , text = element_text(size = 5, family = "Arial") # change font size of all text
    , axis.text = element_text(size = 5, family = "Arial") # change font size of axis text
    , axis.title = element_text(size = 6, family = "Arial") # change font size of axis titles
    , plot.title = element_text(size = 7, family = "Arial") # change font size of plot title
    , legend.text = element_text(size = 6, family = "Arial") # change font size of legend text
    , legend.title = element_text(size = 6, family = "Arial") # change font size of legend title
    , panel.background = element_rect(fill="transparent")
    , plot.background = element_rect(fill="transparent", color = NA)
  )
)
#
