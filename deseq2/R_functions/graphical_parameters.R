# colors ####

suppressPackageStartupMessages(library("scico"))
suppressPackageStartupMessages(library("circlize"))
col_muted <- c("#CC6677", "#88CCEE", "#DDCC77", "#117733", "#332288", "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD") # source is https://personal.sron.nl/~pault/
col_vibrant <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB')
col_high_contrast <- c("#FFFFFF", '#004488', '#DDAA33', '#BB5566', '#000000')
col_bright <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
col_rlog <- colorRamp2(seq(from=1, to=10, by=1), scico(n=10, direction=-1, palette="lajolla"))
many_colors <- c(col_vibrant, col_high_contrast[-1], col_bright[-7], col_muted)


#