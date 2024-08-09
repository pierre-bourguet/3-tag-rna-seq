## I want to create a scatterplot to compare counts from salmon pseudomapping with counts from STAR + salmon counts

# import STAR and salmon sense counts
STAR_counts_sense <- as_tibble(
  read.delim(
    "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/STAR_test_02/02_counts/cdca7_ab_R1/STAR_mapping_salmon_quant_cdca7_ab_R1/quant.sf",
    header = TRUE,
    sep = "\t"
  )
)
salmon_counts_sense <- as_tibble(
  read.delim(
    "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/STAR_test_02/02_counts/cdca7_ab_R1/quant_cdca7_ab_R1/quant.sf",
    header = TRUE,
    sep = "\t"
  )
)
  
# merge the two tables
counts_comparison <- inner_join(STAR_counts_sense, salmon_counts_sense, by = "Name", suffix = c("_STAR", "_salmon"))

# create a scatterplot of all isoforms
ggplot(counts_comparison, aes(x = log2(NumReads_STAR+1), y = log2(NumReads_salmon+1))) +
  geom_point(alpha=0.25, size=0.75) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = "STAR counts (log2+1)", y = "Salmon counts (log2+1)", title = "Comparison of counts from STAR and Salmon\nall isoforms\nsample:cdca7_ab_R1 from tagseq_03") +
  coord_fixed() +
  theme_minimal()

# Merge the counts at the gene level
counts_comparison_merged <- counts_comparison %>%
  dplyr::select(Name, NumReads_STAR, NumReads_salmon) %>%
  dplyr::mutate(Name = gsub("\\.\\d+", "", Name)) %>%  # Remove the .1, .2, etc. at the end of Geneid
  dplyr::group_by(Name) %>%
  dplyr::summarise(across(everything(), sum)) # Sum the counts of all isoforms for each gene

# create a scatterplot of merged isoforms
ggplot(counts_comparison_merged, aes(x = log2(NumReads_STAR+1), y = log2(NumReads_salmon+1))) +
  geom_point(alpha=0.5, size=1) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = "STAR counts (log2+1)", y = "Salmon counts (log2+1)", title = "Comparison of counts from STAR and Salmon\nmerged isoforms\nsample:cdca7_ab_R1 from tagseq_03") +
  coord_fixed() +
  theme_minimal()

# compare how many loci have 0 counts in STAR & salmon
counts_comparison_merged %>%
  dplyr::filter(NumReads_STAR == 0) %>%
  nrow() # 0
counts_comparison_merged %>%
  dplyr::filter(NumReads_salmon == 0) %>%
  nrow() # 0
