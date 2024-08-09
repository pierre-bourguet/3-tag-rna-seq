# import DESeq2 environment
load("../06_DESeq2/DESeq2_environment.RData")
setwd("../07_analysis/")

# create output directory
dir.create("figures", showWarnings = F, recursive = T)

#
# import libraries ####
library(ggbreak)
library(patchwork)
library("pheatmap")
library(svglite)
library(ComplexHeatmap)
library(agricolae)
library(ggbeeswarm)

#
# number of up TEs in mutants ####

# find TEs upregulated in each mutant
up_TEs <- lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features")

# extract nb of upregulated TEs for each genotype and convert to a dataframe
nb_up_TEs <- lapply(FUN = length, X = up_TEs)
nb_up_TEs_df <- data.frame(
  condition = gsub("res_", "", names(nb_up_TEs)),
  Count = unlist(nb_up_TEs)
)

#
# histogram: number of up TEs in cdca7 mutants ####

# look at all the data
ggplot(nb_up_TEs_df, aes(y = condition, x = Count)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  theme_minimal()

# filter out ddm1 samples
nb_up_TEs_df_no_ddm1 <- nb_up_TEs_df %>% filter(!str_detect(condition, "ddm1_[^2]|F2")) %>%
  filter(!str_detect(condition, "ddm1_2_G5"))

# reorder levels
nb_up_TEs_df_no_ddm1$condition <- factor(nb_up_TEs_df_no_ddm1$condition, levels = rev(c("a_1", "a_2", "b_1", "b_2", "a_long_1", "a_long_2", "a_long_3", "ab_1", "ab_2", "a_long_b", "ddm1_2_G2")))

# add a WT condition
fig1_nb_up_TEs_df <- nb_up_TEs_df_no_ddm1 %>%
  filter(!str_detect(condition, "long")) %>%
  bind_rows(tibble(condition = "WT", Count = NA))

# reorder levels
fig1_nb_up_TEs_df$condition <- factor(fig1_nb_up_TEs_df$condition, levels = rev(c("WT", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))

fig1_up_TEs_hist_WT <- ggplot(fig1_nb_up_TEs_df, aes(y = condition, x = Count, fill = condition)) +
  geom_col() +
  labs(x = "Count", y = "Genotype") +
  scale_y_discrete(position = "left") +
  scale_fill_manual(values = rev(col_muted_2_replicates[c(19,1:4,7,8,5)])) +
  scale_x_continuous(breaks = c(0,100,200,300), limits=c(0,1130), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(350, 990), ticklabels=c(1000,1100), scales=0.3, space = 0.1, expand = F) +
  theme_horizontal_nature +
  theme(
    panel.border = element_rect(fill=NA, colour = "black", linewidth = pt_0.5_to_mm),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(margin = margin(t = -5)), # reduce space between axis titles and axis labels
    axis.title.y = element_text(margin = margin(r = -5), angle = 90),
    plot.margin = margin(t = -5, r = -5)
    , legend.position = "none"
  ) ; fig1_up_TEs_hist_WT

# export in svg, works fine but fonts are not exported properly (open the svg file in a text editor to see).
svglite::svglite(filename = "figures/up_TEs_cdca7_histogram_with_WT.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
fig1_up_TEs_hist_WT
dev.off()

#
# histogram: number of up TEs in ddm1 cdca7 mutants ####

# filter samples
nb_up_TEs_df_w_ddm1 <- nb_up_TEs_df %>%
  # filter(str_detect(condition, "ddm1|ab")) %>% # use this one include cdca7ab controls
  filter(str_detect(condition, "ddm1")) %>%
  filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
  filter(!condition %in% "F2_ddm1") %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_ddm1 <- ggplot(nb_up_TEs_df_w_ddm1, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,4,3,5)], rep(col_muted[5], 5)))) +
  theme(legend.position = "none") ; barplot_nb_up_TEs_ddm1

#
# superplot: rlog of TEs up in ddm1 cdca7 mutants ####

################## isolate TEs up in samples of interest

# Filter samples of interest for the plot
filtered_names <- names(up_TEs)[
  grepl("ddm1|Col_0|ab", names(up_TEs)) &                # Include if matches ddm1, Col_0, or ab
    !grepl("long|NBD|ddm1_2_G2", names(up_TEs)) &               # Exclude if matches long, NBD, or ddm1_2_G2
    !(names(up_TEs) %in% "F2_ddm1")                  # Exclude specific name F2_ddm1
]

# Select elements from the list using the filtered names
up_TEs_ddm1_cdca7 <- up_TEs[filtered_names]
up_TEs_ddm1_cdca7_ids <- unique(unlist(up_TEs_ddm1_cdca7))

################## filter normalized counts and format for superplot

# Specify the normalization method
normalization_method <- "rlog"  # "RPM" or "rlog"

# Prepare the matrices and other variables based on the normalization method
if (normalization_method == "RPM") {
  z <- RPM_merged %>%
    tidyr::pivot_longer(-Geneid, names_to = "sample", values_to = "expression") %>%
    dplyr::mutate(condition = gsub("_R[0-9]", "", sample)) %>%
    dplyr::filter(str_detect(condition, "ddm1|Col_0|ab")) %>%
    dplyr::filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
    dplyr::filter(!condition %in% "F2_ddm1")
  # Select up TEs & log2 transform
  m <- z %>%
    dplyr::filter(Geneid %in% up_TEs_ddm1_cdca7_ids) %>%
    dplyr::mutate(expression = log2(expression + 1))
} else if (normalization_method == "rlog") {
  z <- rlog %>%
    tidyr::pivot_longer(-Geneid, names_to = "sample", values_to = "expression") %>%
    dplyr::mutate(condition = gsub("_R[0-9]", "", sample)) %>%
    #dplyr::filter(str_detect(condition, "ddm1|Col_0|ab")) %>% # use this one to include Col-0 & ab controls
    dplyr::filter(str_detect(condition, "ddm1|Col_0")) %>% # use this one to include Col-0 only
    dplyr::filter(!str_detect(condition, "long|NBD|ddm1_2_G2")) %>%
    dplyr::filter(!condition %in% "F2_ddm1")
  # Select up TEs
  m <- z[z$Geneid %in% up_TEs_ddm1_cdca7_ids,]
}

# import the superplot function
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/superplot_w_boxplot.R")

# Define the specific order vector
order_vector <- c("Col_0", "ab_1", "ab_2" , "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")
order_vector <- c("Col_0", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")

# Define the color palette
my_colors <- c("grey55", col_muted[c(4,4,3)], rep(col_muted[5], 6)) # to include cdca7ab controls
my_colors <- c("grey55", col_muted[c(3)], rep(col_muted[5], 6)) # without
# rename the dataframe for compatibility with the function
names(m)[3] <- "value"

# run the function and customize the plot
superplot_upTEs_ddm1 <- superplot_w_boxplot(m, rev(order_vector), rev(my_colors)) + 
  coord_flip(ylim=c(0,10)) +
  xlab("Genotype") + ylab("Transcript levels (rlog)") + labs(title="upregulated TEs")
superplot_upTEs_ddm1

# write svg output
svglite::svglite(filename = "figures/fig_rlog_TEs_ddm1_cdca7.svg", width = 80*mm_to_inches, height = 50*mm_to_inches)
set.seed(55)
barplot_nb_up_TEs_ddm1 + superplot_upTEs_ddm1 + patchwork::plot_layout(axes = 'collect') & theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- m %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
# cdca7 long mutants: number of up TEs, heatmap and superplots ####

up_TEs <- lapply(FUN = length, X = (lapply(FUN = up, X = all_res, y="TEs")))
nb_up_TEs_df <- data.frame(
  condition = gsub("res_", "", names(up_TEs)),
  Count = unlist(up_TEs)
)

# filter samples
nb_up_TEs_df_a_long_no_ddm1 <- nb_up_TEs_df %>% filter(str_detect(condition, "long|a_1|a_2|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

nb_up_TEs_df_a_long_only_ab <- nb_up_TEs_df %>% filter(str_detect(condition, "a_long_b|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0)) %>%
  mutate(condition = factor(condition, levels = rev(c("Col_0", "ab_1", "ab_2", "a_long_b")))) %>%
  arrange(condition)
  
nb_up_TEs_df_a_long_w_ddm1 <- nb_up_TEs_df %>% filter(str_detect(condition, "ddm1") & str_detect(condition, "G5|long|ab")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
# nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_a_long <- ggplot(rbind(nb_up_TEs_df_a_long_no_ddm1, nb_up_TEs_df_a_long_w_ddm1), aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  #scale_fill_manual(values = rev(c(col_muted[c(10,4,4,3,5)], rep(col_muted[5], 5)))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long

barplot_nb_up_TEs_a_long_b <- ggplot(nb_up_TEs_df_a_long_only_ab, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,4,6)]))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long_b

## heatmap at TEs up in cdca7a/b
cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7a_long_mutants"
            , z = rld_df %>% dplyr::select(matches("long_b|ab|Geneid") & -matches("ddm1"))
            , n = "rlog")
wide_rld_avg

# superplot
#cdca7_ab_up_TEs_union <- union(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))

a_long_b_rld <- rld_df %>%
  dplyr::filter(Geneid %in% cdca7_ab_up_TEs) %>%
  dplyr::select(matches("long_b|ab|Geneid|Col_0") & -matches("ddm1")) %>%
  tidyr::pivot_longer(cols = -Geneid, names_to = "sample", values_to = "value") %>%
  dplyr::mutate(condition=stringr::str_replace(sample, "_R[123]", ""))

superplot_a_long_b <- superplot_w_boxplot(data = a_long_b_rld, condition_order = rev(c("Col_0", "ab_1", "ab_2", "a_long_b"))
                                          , colors = rev(c("grey55", col_muted[c(4,4,6)]))) + 
  coord_flip(ylim=c(1.2,9)) +
  xlab("Genotype") + ylab("Transcript levels (rlog)") + labs(title="upregulated TEs")

# write svg output
svglite::svglite(filename = "figures/superplot_a_long_b.svg", width = 3, height = 1.5)
barplot_nb_up_TEs_a_long_b + superplot_a_long_b + plot_layout(axes = "collect", widths = c(2,3)) #& theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- a_long_b_rld %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
# heatmap of upregulated TEs ####
upTEs_RPM <- as_tibble(read.delim("DEGs_batch/normalized_counts/batch_upTEs_RPM.tsv", header=T, sep="\t"))
upTEs_rlog <- as_tibble(read.delim("DEGs_batch/normalized_counts/batch_upTEs_rlog.tsv", header=T, sep="\t"))

# filter out samples with ddm1 mutations, mom1 mutants, F2 segregants
upTEs_RPM_no_ddm1 <- upTEs_RPM[, grep("ddm1|F2|mom1", colnames(upTEs_RPM), invert = T)]
upTEs_rlog_no_ddm1 <- upTEs_rlog[, grep("ddm1|F2|mom1", colnames(upTEs_rlog), invert = T)]

# heatmap of upregulated TEs
normalization <- "RPM"
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=RPM_merged, n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=rld_df, n = "rlog") ) # DOESN'T WORK FOR SOME REASON

### TEs up regulated in cdca7a/b or ddm1

cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
ddm1_G2_up_TEs <- up(all_res$res_ddm1_2_G2, "TEs")

## using log2 (RPM + 1) & rlog at cdca7 upTEs

# with all samples
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs", z = RPM_merged, n = normalization)

# only cdca7 mutants
# log2(RPM+1)
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants", z = RPM_merged %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = normalization)
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants_avg", z = RPM_merged_avg %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = normalization)
# rlog
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants", z = rld_df %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = 'rlog')
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7_mutants_avg", z = wide_rld_avg %>% dplyr::select(-matches("ddm1|F2|mom1|long")), n = 'rlog')

## heatmap of rlog at cdca7-upTEs for fig1
# Specify the normalization method
normalization_method <- "rlog"  # "RPM" or "rlog"

# Prepare the matrices and other variables based on the normalization method
if (normalization_method == "RPM") {
  z <- RPM_merged_avg %>%
    dplyr::select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
  # Select up TEs & log2 transform & remove Geneid
  m <- as.matrix(log2(z[z$Geneid %in% cdca7_ab_up_TEs, -1] + 1))
  color_scale <- colorRamp2(seq(from = 0, to = 9, by = 1), scico(n = 10, direction = -1, palette = "lajolla"))
  heatmap_name <- "log2\n(RPM+1)"
} else if (normalization_method == "rlog") {
  z <- wide_rld_avg %>%
    dplyr::select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
  # Select up TEs & remove Geneid
  m <- as.matrix(z[z$Geneid %in% cdca7_ab_up_TEs, -1])
  color_scale <- colorRamp2(seq(from = 3, to = 10, by = 1), scico(n = 8, direction = -1, palette = "lajolla"))
  heatmap_name <- "rlog"
}

# See data spread to define the color scale
summary(m)

# Prepare labels
italicized_labels_heatmap <- c(
  expression(plain('WT')),
  expression(italic('cdca7α-1')),
  expression(italic('cdca7α-2')),
  expression(italic('cdca7β-1')),
  expression(italic('cdca7β-2')),
  expression(italic('cdca7α/β-1')),
  expression(italic('cdca7α/β-2')),
  expression(italic('ddm1'))
)

# Plot the heatmap
Heatmap(m,
        name = heatmap_name,
        col = color_scale,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width = 15, height = 15,
        column_title = paste0("cdca7α/β upTEs\nn=", length(cdca7_ab_up_TEs)),
        column_labels = italicized_labels_heatmap
)

# remove clustering & sort the heatmap by row sums (TEs with highest expression across genotypes at the top)
m <- m[order(rowSums(m), decreasing = T),] 
fig1_up_TEs_heatmap <- Heatmap(t(m),
                                name = heatmap_name,
                                col = color_scale,
                                cluster_rows = F,
                                cluster_columns = F,
                                width=15, height=15,
                                column_title=NULL,
                                row_labels = italicized_labels_heatmap,
                                row_names_side = "left",
                                row_names_gp = grid::gpar(fontsize = 6),
                                use_raster = T,
                                border= T,
                                heatmap_legend_param = list(
                                  title = "log2\n(RPM+1)", at = c(0, 5, 10), 
                                  labels = c("0", "5", "10"),
                                  legend_height = unit(1, "cm"),
                                  legend_width = unit(1, "cm"),
                                  labels_gp = gpar(fontsize = 6),
                                  title_gp = gpar(fontsize = 6)
                                )
) ; fig1_up_TEs_heatmap

svglite::svglite(filename = "figures/up_TEs_cdca7_heatmap.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
svglite::svglite(filename = "figures/up_TEs_cdca7_heatmap_rlog.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
draw(fig1_up_TEs_heatmap)
dev.off()

## using log2 (RPM + 1) at ddm1 upTEs

# select samples of interest
z <- RPM_merged_avg %>%
  select(any_of(c("Geneid", "Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))
# log2 transform
m <- as.matrix( log2(z[z$Geneid %in% ddm1_G2_up_TEs,-1] + 1))
# see data spread to scale the heatmap
summary(m)
col_log2_RPM <- colorRamp2(seq(from=0, to=9, by=1), scico(n=10, direction=-1, palette="lajolla"))
Heatmap(m,
        name = "log2\n(RPM+1)",
        col = col_log2_RPM,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15,
        column_title=paste0("ddm1-2 upTEs\nn=", length(ddm1_G2_up_TEs))
)

## Z-score heatmap, scaling on RPM
z <- RPM_merged_avg %>%
  select(-matches("ddm1|F2|mom1|long"))
col_order <- colnames(m)[c(1:3,6,7,4,5)]

m <- as.matrix(z[z$Geneid %in% cdca7_ab_up_TEs,-1])
m <- t(scale(t(m)))
summary(m)
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
Heatmap(m,
        name = "z-score",
        col = col_zscore,
        column_order = col_order, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15
)

# z-score scaling after log2 + 1 transformation
m <- as.matrix( log2(z[z$Geneid %in% cdca7_ab_up_TEs,-1] + 1))
m <- t(scale(t(m)))
summary(m)
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
Heatmap(m,
        name = "z-score",
        col = col_zscore,
        column_order = col_order, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15
)

# Save the heatmap to a PDF file
ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
file_name <- paste0(output_dir, n, "/heatmap_", y, ".pdf")
pdf(file_name, width = 15, height = 7)
draw(ht, heatmap_legend_side = "right")
dev.off()

# heatmap at ddm1 up TEs
ddm1_G5_up_TEs <- up(all_res$res_ddm1_2_G5, "TEs")
# with all samples
DEG_heatmap(ddm1_G5_up_TEs, "ddm1_G5", z = RPM_merged, n = normalization)

#

# looking at TEs ####

# all upregulated TEs at all samples

d <- RPM %>%
  filter(Geneid %in% DEGs$upTEs) %>%
  pivot_longer(cols=-Geneid) %>%
  mutate(condition=substr(name, 1, nchar(name)-3))  %>%
  mutate(log2_RPM = log2(value + 1)) %>% # distinguish ddm1 conditions from others
  mutate(ddm1 = if_else(str_detect(condition, "ddm1"), TRUE, FALSE))

ggplot(d %>% filter(ddm1 == T), aes(x = name, y = log2_RPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = paste0("quantification of TE transcripts\nn = ", length(unique(d$Geneid))), y = "Log2 (RPM + 1)", x = "Genotype") +
  coord_flip() +
  theme(legend.position = "none")

ggplot(d %>% filter(ddm1 == F), aes(x = name, y = log2_RPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = paste0("quantification of TE transcripts\nn = ", length(unique(d$Geneid))), y = "Log2 (RPM + 1)", x = "Genotype") +
  coord_flip() +
  theme(legend.position = "none") +
  ylim(0,8)

# Calculating the median of log2_RPM for each name and merge with the transgene expression levels

medians <- d %>%
  group_by(name, condition) %>%
  summarise(median_log2_RPM = median(log2_RPM, na.rm = TRUE), .groups = 'drop')

# merge with transgene expression levels
mTurq_data <- RPM %>%
  filter(Geneid %in% "mTurq_3xcMyc") %>%
  pivot_longer(cols=-Geneid) %>%
  mutate(condition=substr(name, 1, nchar(name)-3)) %>%
  mutate(log2_RPM = log2(value + 1)) %>%
  select(name, mTurq_log2_RPM = log2_RPM)

medians <- left_join(medians, mTurq_data, by = c("name"))

# plot all data
ggplot(medians, aes(y = condition, x = median_log2_RPM, fill = condition, size = mTurq_log2_RPM)) +
  geom_jitter(width = 0.01, height = 0.2, shape = 21, show.legend = F) +
  theme_minimal(base_family = "Arial") +  
  theme(text = element_text(family = "Arial")) +
  labs(title = paste0("upregulated TE transcripts\nn = ", length(DEGs$upTEs)),
       y = "", x = "Median log2 (RPM + 1)", fill="genotype")# +
scale_fill_manual(values=col_muted[c(4,2,1,10)]) +
  scale_y_discrete(labels= rev(italicized_labels_tagseq)) +
  theme(axis.title.y = element_blank())


# all upregulated TEs in samples of interest ####

samples <- c("WT", "cdca7_a", "cdca7_b", "cdca7_ab")

dd <- d %>%
  filter(condition %in% samples)

ggplot(dd, aes(x = name, y = log2_RPM, fill = condition)) +
  geom_boxplot(outlier.shape=NA) +
  theme_minimal() +
  labs(title = "quantification of TE transcripts", y = "log2 (RPM + 1)", x = "genotype") +
  coord_flip()

# now plotting the median of replicates
# Filtering and ordering data based on specific samples
dd <- d %>%
  filter(condition %in% samples) %>%
  mutate(condition = factor(condition, levels = rev(samples)))

# Calculating the median of log2_RPM for each name
medians <- dd %>%
  group_by(name, condition) %>%
  summarise(median_log2_RPM = median(log2_RPM, na.rm = TRUE), .groups = 'drop')

# Creating the plot with ordered conditions
italicized_labels_tagseq <- c(
  expression(plain('WT')),
  expression(italic('cdca7-α')),
  expression(italic('cdca7-β')),
  expression(italic('cdca7-α/β'))
)
grDevices::cairo_pdf("figures/CDCA7-ab_log2_RPM_at_TEs.pdf", width = 4, height = 2)
set.seed(130)
p2 <- ggplot(medians, aes(y = condition, x = median_log2_RPM, fill = condition)) +
  geom_jitter(width = 0.01, height = 0.2, shape = 21, size = 3, show.legend = F) +
  theme_minimal(base_family = "Arial") +  
  theme(text = element_text(family = "Arial")) +
  labs(title = paste0("upregulated TE transcripts\nn = ", length(DEGs$upTEs)),
       y = "", x = "Median log2 (RPM + 1)", fill="genotype") +
  scale_fill_manual(values=col_muted[c(4,2,1,10)]) +
  scale_y_discrete(labels= rev(italicized_labels_tagseq)) +
  theme(axis.title.y = element_blank()) ; p2
dev.off()

library(patchwork)
setwd("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/tagseq_03_cdca7/")
set.seed(130)
magnifier <- 1.5 ; grDevices::cairo_pdf("figures/poster_fig2A_2B.pdf", width = 5*magnifier, height = 4*magnifier)
p1 + p2 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = "bold"),
                                                    text = element_text(size=14, family="Arial"),
                                                    axis.text.y = element_text(size=16)
)
dev.off()



# make linear correlations of rlogs between samples ####

### between cdca7-ab mutants
cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
cdca7_ab_up_PCGs <- intersect(up(all_res$res_ab_1, "PCGs"), up(all_res$res_ab_2, "PCGs"))
ddm1_up_PCGs <- up(all_res$res_ddm1_2_G2, "PCGs")
cdca7_ab_rld <- wide_rld_avg %>%
  dplyr::filter(Geneid %in% cdca7_ab_up_TEs)

# linear correlation between samples
cor(cdca7_ab_rld %>% dplyr::select(-Geneid), method = "pearson")
# scatterplot of cdca7-ab mutants
ggplot(cdca7_ab_rld, aes(x = ab_1, y = a_long_b)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  geom_abline() +
  labs(title = "cdca7-ab mutants", x = "rlog(ab_1)", y = "rlog(ab_2)") +
  coord_fixed(ylim=c(2,11), xlim=c(2,11))

data <- cdca7_ab_rld %>% dplyr::select(c("Col_0", "ab_1", "a_long_b", "ab_2", "ddm1_2_G2", "ddm1_2_G5"))
data <- cdca7_ab_rld %>% dplyr::select(c("Col_0", "a_1", "a_2", "a_long_1", "a_long_2", "a_long_3"))
pairs(data)
library(GGally)
ggpairs(data, title = "Scatter Plot Matrix for mtcars Dataset", axisLabels = "show") +
  coord_fixed(ylim=c(1,11), xlim=c(1,11)) +
  geom_abline(intercept = 0, slope =1)


#
### NOT UPDATED from there ####

grDevices::cairo_pdf("figures/CDCA7-ab_log2_RPM_at_TEs.pdf", width = 4, height = 2)
set.seed(130)

# plot heatmaps ####
suppressPackageStartupMessages(library("RColorBrewer")) ; suppressPackageStartupMessages(library("pheatmap"))
DEG_heatmap <- function(x, y, z, n) {
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      if (n=="rlog") {
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\n", n)
      }
      else {
        sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x)+1))
        title <- paste0(y, "\nn=", length(x),"\nlog2(", n, "+1)")
      }
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(output_dir, n, "/heatmap_", y, ".pdf"), main=title, cluster_cols = F)
    }
  }
}

mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,9:(9+length(conditions))], n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_norm[,9:(8+nrow(samples))], n = normalization) )


# write file with RPM / RPM values for all samples, averaged over replicates ####
write.table(x=cts_summary_norm_mean, file=paste0(args[1], normalization, strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# function that adds annotations to a df using its row.names as Geneid
annotate_df <- function(df) { 
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}
# median of ratios (from DESeq2, median of ratios to geometric mean). Write tables and heatmaps ####
cts_MoR <- annotate_df( as.data.frame(counts(dds, normalized=T)) )
cts_MoR_mean <- average_replicates(cts_MoR)
write.table(x=cts_MoR, file=paste0(args[1], "ESF", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=cts_MoR_mean, file=paste0(args[1], "ESF", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_MoR[,9:(8+nrow(samples))], n = "MoR") )
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_MoR_mean[,9:(9+length(conditions))], n = "MoR") )
# exploring differences between median of ratios normalization by DESeq and RPM / RPM ####
if (FALSE==TRUE) { # just protecting this code so it's not executed when i run the script
  test <- cts_summary[,sample_columns]
  par(pty = "s")
  plot(x=colSums(test) / colSums(test)[6], y=sizeFactors(dds), xlim=c(0.5,1.7), ylim=c(0.5,1.7), xlab=normalization)
  abline(a = 0, b=1)
  text(x=colSums(test) / colSums(test)[6], sizeFactors(dds), labels=names(sizeFactors(dds)), cex= 0.5, pos=1)
}
