library(tidyverse)
library(ggplot2)
library(reshape2)

# Check if command-line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  stop("Usage: Rscript script.R <filename>")
}

# Get the filename from the command-line arguments
raw_counts <- commandArgs(trailingOnly = TRUE)[1]
trimmed_counts <- commandArgs(trailingOnly = TRUE)[2]
umi_removed_counts <- commandArgs(trailingOnly = TRUE)[3]
xaxis_order <- commandArgs(trailingOnly = TRUE)[4]

if (xaxis_order == 1){
	raw_rc <- read_tsv(raw_counts, col_names = c("Sample","Reads")) %>% 
		mutate(Sample = str_sub(Sample, end=-13))
	fastp_trimmed_rc <- read_tsv(trimmed_counts, col_names = c("Sample","Reads")) %>% 
  		mutate(Sample = str_sub(Sample, end=-13))
	umi_collapsed_rc <- read_tsv(umi_removed_counts, col_names = c("Sample","Reads")) %>% 
		mutate(Sample = str_sub(Sample, end=-13))
} else {
	raw_rc <- read_tsv(raw_counts, col_names = c("Sample","Reads"))
        fastp_trimmed_rc <- read_tsv(trimmed_counts, col_names = c("Sample","Reads")) 
        umi_collapsed_rc <- read_tsv(umi_removed_counts, col_names = c("Sample","Reads"))
}

all_rc <- raw_rc %>% left_join(fastp_trimmed_rc, by="Sample") %>% 
  left_join(umi_collapsed_rc, by="Sample") %>% 
  rename(Raw=Reads.x, Trimmed=Reads.y, Umi_Collapsed=Reads)

write_tsv(all_rc, file = "all_read_counts_summarized.txt")

g3 <- melt(all_rc)


if (xaxis_order == 1){
  pdf(file = "all_read_counts_summarized.pdf")
  print(
    my_plot <- g3 %>% ggplot(aes(x=factor(Sample, levels=Sample[c(7,8,9,10,17,18,1,2,3,4,15,16,19,20,21,22,11,12,13,14,5,6)]), y = value, fill = variable))+
      geom_bar(stat = "identity", position=position_dodge(width=0.8))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      xlab("Genotype")+
      ylab("# of Reads")+
      scale_y_continuous(labels = scales::comma_format())+
      scale_fill_manual(values = c("black", "grey", "green4"))+
      labs(fill="")
  )
  dev.off()
} else {
  pdf(file = "all_read_counts_summarized.pdf")
  print(
    my_plot <- g3 %>% ggplot(aes(x=factor(Sample), y = value, fill = variable))+
        geom_bar(stat = "identity", position=position_dodge(width=0.8))+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        xlab("Genotype")+
        ylab("# of Reads")+
        scale_y_continuous(labels = scales::comma_format())+
        scale_fill_manual(values = c("black", "grey", "green4"))+
        labs(fill="")
        )
  dev.off()
}
