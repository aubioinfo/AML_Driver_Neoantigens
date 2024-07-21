#### CIBERSORTx Analysis ####

# Load necessary libraries
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggrastr)
library(reshape2)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\06_clinicalData\\CIBERSORTx")

# Load and preprocess data
dat <- read.csv("plot_cibersortx.csv", header = TRUE)
dat_long <- melt(dat)
dat_filtered <- subset(dat_long, variable %in% c('Tn_MAL', 'Tm_CD52', 'Tex'))

# Create violin plot with additional layers
p <- ggplot(dat_filtered, aes(x = group, y = value, color = group)) + 
  geom_violin() +
  geom_quasirandom_rast(width = 0.25, alpha = 0.5) +
  geom_boxplot_jitter(width = 0.4, alpha = 0.75) + 
  theme_pubr() + 
  facet_wrap(~ variable, scales = "free") +
  labs(x = "", y = "Predicted cellular abundance") +
  theme(strip.text.x = element_text(size = 13),
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") +
  scale_color_manual(values = c("#7570B3", "#1B9E77")) +
  stat_compare_means(comparisons = list(c('Yes', 'No')), 
                     method = 'wilcox.test', 
                     label = "p.signif")

# Save plot
ggsave(p, filename = "01_AMLHighLow_cibersortx.pdf", width = 5, height = 3.3)