###### Neoantigen presentation features
###### Jin Peng 
###### 2021.12.15

# Set global options and load libraries
options(scipen = 10)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(ggpubr)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\04_PresentationFeature")

# Load data
t1 <- read_xlsx("01.Unique_pMHC_DriverMutation_Fusions.xlsx")
t1 <- as.data.frame(t1)

# Define theme for plots
custom_theme <- theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")
  )

# Plot histograms
p1 <- ggplot(t1, aes(x = MT_Score, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 20) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Binding Affinity (nM)") +
  custom_theme

p2 <- ggplot(t1, aes(x = Stability, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_log10() +
  labs(title = "Binding Stability (Hours)") +
  custom_theme

p3 <- ggplot(t1, aes(x = Hydrophobicity_fraction, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 10) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Hydrophobicity Fraction") +
  custom_theme

# Split data for further analysis
t_mutation <- subset(t1, Variant_type != "Fusion")
t_fusion <- subset(t1, Variant_type == "Fusion")

p4 <- ggplot(t_mutation, aes(x = expression_tpm, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 25) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_log10() +
  labs(title = "Tumor Abundance") +
  custom_theme

p5 <- ggplot(t_fusion, aes(x = expression_tpm, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 20) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_log10() +
  labs(title = "Tumor Abundance") +
  custom_theme

# Combine plots and save to file
combined_plots <- p1 + p2 + p3 + p4 + p5
ggsave(combined_plots, filename = "01.presentation_features.pdf", width = 6, height = 4.4)

# Mutation position analysis
m1 <- unique(t_mutation[, c("MT_Epitope_Seq", "HLA_Allele", "Peptide_Length", "Mutation_Position")])
m2 <- data.frame(table(m1[, c("Peptide_Length", "Mutation_Position")]))
m3 <- subset(m2, Mutation_Position != 0)

p_cor <- ggplot(m3, aes(x = Mutation_Position, y = Peptide_Length)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradientn(colors = c("#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61")) +
  labs(title = "Mutation Position") +
  custom_theme

ggsave(p_cor, filename = "01.mutation_position.pdf", width = 3.8, height = 3)

# Variant type presentation features
t1 <- read_xlsx("figure2_plotdata.xlsx")
t1 <- as.data.frame(t1)
t2 <- subset(t1, Variant_type != "Fusion")

# Plot Binding Affinity
p1 <- ggplot(t1, aes(x = Variant_type, y = MT_Score, fill = Variant_type)) +
  geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Binding Affinity (nM)") +
  custom_theme +
  stat_compare_means(aes(x = Variant_type, y = MT_Score), ref.group = "Indel")

# Plot Binding Stability
p2 <- ggplot(t1, aes(x = Variant_type, y = log10(Stability), fill = Variant_type)) +
  geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Binding Stability (Hours)") +
  custom_theme +
  stat_compare_means(aes(x = Variant_type, y = log10(Stability)), ref.group = "Indel")

# Plot Tumor Abundance
p3 <- ggplot(t2, aes(x = Variant_type, y = expression_tpm, fill = Variant_type)) +
  geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Tumor Abundance (TPM)") +
  scale_y_log10() +
  custom_theme +
  stat_compare_means(aes(x = Variant_type, y = expression_tpm), ref.group = "Indel")

# Plot Hydrophobicity Fraction
p4 <- ggplot(t1, aes(x = Variant_type, y = Hydrophobicity_fraction, fill = Variant_type)) +
  geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Fraction Hydrophobicity") +
  custom_theme +
  stat_compare_means(aes(x = Variant_type, y = Hydrophobicity_fraction), ref.group = "Indel")

# Combine variant type plots and save to file
variant_type_plots <- p1 + p2 + p3 + p4
ggsave(variant_type_plots, filename = "02.presentation_features_snv_indel_fusion.pdf", width = 5, height = 4.8)

## immunogenicity score for fusion, indel and snp

# Load required libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\04_PresentationFeature")

# Load data
t1 <- read_xlsx("figure2_plotdata.xlsx")
t1 <- as.data.frame(t1)

# Define factor levels and comparisons
t1$Variant_type <- factor(t1$Variant_type, levels = c("Indel", "Fusion", "SNP"))
my_comparisons <- list(c("Indel", "SNP"), c("Indel", "Fusion"), c("SNP", "Fusion"))

# Define a custom color palette
palette <- c("#FC8D62", "#66C2A5", "#8DA0CB")

# Plot PRIME immunogenicity
p0 <- ggviolin(
  t1, x = "Variant_type", y = "PRIME", 
  fill = "Variant_type", xlab = "",
  palette = palette,
  add = "boxplot", 
  add.params = list(fill = "white")
) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

# Plot DeepImmuno immunogenicity
p1 <- ggviolin(
  t1, x = "Variant_type", y = "DeepImmuno_immunogenicity", 
  fill = "Variant_type", xlab = "",
  palette = palette,
  add = "boxplot", 
  add.params = list(fill = "white")
) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

# Plot DeepHLApan immunogenicity
p2 <- ggviolin(
  t1, x = "Variant_type", y = "DeepHLApan_immunogenic_score", 
  fill = "Variant_type", xlab = "",
  palette = palette,
  add = "boxplot", 
  add.params = list(fill = "white")
) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

# Plot IEDB immunogenicity
p3 <- ggviolin(
  t1, x = "Variant_type", y = "IEDB_Immunogenicity", 
  fill = "Variant_type", xlab = "",
  palette = palette,
  add = "boxplot", 
  add.params = list(fill = "white")
) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

# Combine plots
plots <- p0 + p2 + p3

# Save combined plots to file
ggsave(plots, filename = "04.Immunogenicity_DeepImmuno_DeepHLApan_IEDB_snv_indel_fusion.pdf", width = 9.6, height = 4.6)


## Using the TESLA .etc data to train classifier

# Load required libraries
library(pROC)
library(readxl)
library(ggplot2)
library(ggsci)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\04_PresentationFeature\\NeoantigenDB")

# Load data
dat <- read_xlsx("PRIME_tableS1_cell2020_formated.xlsx")
dat <- as.data.frame(dat)

# Define a function for ROC plot
plot_roc <- function(data, method_col, method_name, filename) {
  data_clean <- na.omit(data[, c(1, 7, method_col)])
  roc_obj <- roc(Immunogenicity ~ get(names(data_clean)[3]), data_clean)
  
  pdf(filename, width = 4, height = 4)
  plot(roc_obj, legacy.axes = TRUE,
       main = paste("Optimal cut-off point\n(", method_name, ")", sep = ""),
       thresholds = "best",
       print.thres = "best")
  dev.off()
}

# Plot ROC curves
plot_roc(dat, 10, "PRIME", "Optimal_curoff_PRIME.pdf")
plot_roc(dat, 4, "DeepHLApan", "Optimal_curoff_DeepHLAPan.pdf")
plot_roc(dat, 11, "DeepImmuno", "Optimal_curoff_DeepImmuno.pdf")
plot_roc(dat, 12, "IEDB", "Optimal_curoff_IEDB.pdf")

######
###### Immunogenic yes or no number distribution
######

# Set working directory for plotting
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\04_PresentationFeature")

# Load data for plotting
t <- read_xlsx("01.Unique_pMHC_DriverMutation_Fusions.xlsx")
t1 <- as.data.frame(t)
t2 <- t1[, c("Variant_type", "Immunogenic2")]

# Prepare data for plotting
vartype <- as.data.frame(table(t2))
vartype$x <- "b"
vartype$y <- "a"
vartype$cat1 <- paste0("a_", 1:6)
vartype$cat2 <- vartype$Immunogenic2

# Calculate percentages and labels
dat1 <- aggregate(Freq ~ cat1, data = vartype, sum)
dat1$per1 <- dat1$Freq / sum(dat1$Freq)
dat1$per.y1 <- cumsum(dat1$per1) - dat1$per1 / 2
dat1$label1 <- paste('(', round(dat1$per1 * 100, 2), '%)', sep = '')

dat2 <- aggregate(Freq ~ cat2, data = vartype, sum)
dat2$per2 <- dat2$Freq / sum(dat2$Freq)
dat2$per.y2 <- cumsum(dat2$per2) - dat2$per2 / 2
dat2$label2 <- paste(dat2$cat2, '(', round(dat2$per2 * 100, 2), '%)', sep = '')

vartype <- merge(vartype, dat1[, c("cat1", "per1", "per.y1", "label1")], by.x = "cat1", by.y = "cat1")
vartype <- merge(vartype, dat2[, c("cat2", "per2", "per.y2", "label2")], by.x = "cat2", by.y = "cat2")

# Plot the data
p <- ggplot(vartype) +
  geom_bar(aes(x, per1, fill = cat1), stat = 'identity', width = 0.8, color = 'white') +
  geom_text(aes(x = 2, y = per.y1, label = label1), size = 2.5, color = 'black') +
  geom_bar(aes(y, per2 / 3, fill = cat2), stat = 'identity', width = 1.3) +
  geom_text(aes(x = 1.25, y = per.y2, label = label2), size = 2.5, color = 'black') +
  scale_y_continuous(labels = scales::percent) +
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_igv() + 
  theme(legend.position = 'top')

# Save plot
ggsave(p, filename = "04.1_Immun_Yes_No_Numbers.pdf", width = 4, height = 4.5)


###### Presentation features for Binders (Immunogenic or not)

# Load required libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\04_PresentationFeature")

# Load and prepare data
t1 <- read_xlsx("figure2_plotdata.xlsx")
t1 <- as.data.frame(t1)
t1$Immunogenic2 <- factor(t1$Immunogenic2, levels = c("Yes", "No"))

# Define comparison list
my_comparisons <- list(c("Yes", "No"))

# Subset data
t2 <- subset(t1, Variant_type != "Fusion")
t3 <- subset(t1, Variant_type == "Fusion")

# Define a function for creating violin plots
create_violin_plot <- function(data, y, ylab, palette, scale_y_log = FALSE) {
  p <- ggviolin(data, x = "Immunogenic2", y = y, fill = "Immunogenic2", xlab = "", ylab = ylab, 
                palette = palette, add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons)
  
  if (scale_y_log) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}

# Create plots
p1 <- create_violin_plot(t1, "MT_Score", "Binding Affinity (nM)", c("#7FC97F", "#BEAED4"))
p2 <- create_violin_plot(t1, "Stability", "Binding Stability (Hours)", c("#7FC97F", "#BEAED4"), TRUE)
p3 <- create_violin_plot(t2, "expression_tpm", "Tumor Abundance (TPM)", c("#7FC97F", "#BEAED4"), TRUE)
p3.1 <- create_violin_plot(t3, "expression_tpm", "Tumor Abundance (FFPM)", c("#7FC97F", "#BEAED4"), TRUE)
p4 <- create_violin_plot(t1, "Hydrophobicity_fraction", "Fraction Hydrophobicity", c("#7FC97F", "#BEAED4"))

# Combine and save plots
combined_plots1 <- p1 + p2 + p3 + p3.1 + p4
ggsave(combined_plots1, filename = "05.presentation_features_Immunogenic_Not.pdf", width = 5.6, height = 6)

###### Immunogenicity score for (Immunogenic or not)

# Load and prepare data again
t1 <- read_xlsx("figure2_plotdata.xlsx")
t1 <- as.data.frame(t1)
t1$Immunogenic2 <- factor(t1$Immunogenic2, levels = c("Yes", "No"))

# Create immunogenicity plots
p1 <- create_violin_plot(t1, "DeepImmuno_immunogenicity", "DeepImmuno Immunogenicity", c("#7FC97F", "#BEAED4"))
p2 <- create_violin_plot(t1, "DeepHLApan_immunogenic_score", "DeepHLApan Immunogenicity", c("#7FC97F", "#BEAED4"))
p3 <- create_violin_plot(t1, "IEDB_Immunogenicity", "IEDB Immunogenicity", c("#7FC97F", "#BEAED4"))
p4 <- create_violin_plot(t1, "PRIME", "PRIME Immunogenicity", c("#7FC97F", "#BEAED4"))

# Combine and save plots
combined_plots2 <- p1 + p2 + p3 + p4
ggsave(combined_plots2, filename = "06.Immunogenicity_yes_not.pdf", width = 6.4, height = 9.2)


# mutation position vs. Immunogenicity
# Load required libraries
library(readxl)
library(ggplot2)
library(RColorBrewer)

# Define file path
file_path <- "01.Unique_pMHC_DriverMutation_Fusions.xlsx"

# Load and prepare data
t1 <- read_xlsx(file_path) %>% as.data.frame()
t1 <- subset(t1, Mutation_Position != 0)

# Mutation Position vs. Immunogenicity
t2 <- as.data.frame(table(na.omit(t1[, c("Mutation_Position", "Immunogenic2")])))
t2_F <- subset(t2, Immunogenic2 == "Yes")
t2_T <- subset(t2, Immunogenic2 == "No")
t2_F$all <- sum(t2_F$Freq)
t2_T$all <- sum(t2_T$Freq)
t3 <- rbind(t2_F, t2_T)
t3$ratio <- t3$Freq / t3$all

# Plot Mutation Position vs. Immunogenicity
p <- ggplot(data = t3, mapping = aes(x = Mutation_Position, y = ratio, fill = Immunogenic2)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), alpha = 0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "Mutation Position", y = "Fraction of All") +
  scale_y_continuous(limits = c(0, 0.13)) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

ggsave(p, filename = "07.Mutation_position_Immunogenic_Yes_No.pdf", width = 3, height = 3)

# Fisher exact test for mutation position
options(scipen = 200)
n_yes <- as.data.frame(table(t1$Immunogenic2))[2, 2]
n_no <- as.data.frame(table(t1$Immunogenic2))[1, 2]

stat <- data.frame(matrix(NA, nrow = 11, ncol = 4))
row.names(stat) <- paste0("Position_", 1:11)
colnames(stat) <- c("P.value", "OR", "95%CI_lower", "95%CI_higher")

for (p in 1:11) {
  m1 <- subset(t2, Mutation_Position == p)
  m2 <- matrix(c(m1[2, 3], n_yes - m1[2, 3], m1[1, 3], n_no - m1[1, 3]), ncol = 2, nrow = 2, byrow = TRUE)
  s1 <- fisher.test(m2)
  stat[p, ] <- c(round(s1$p.value, 3), round(s1$estimate, 3), round(s1$conf.int[1], 3), round(s1$conf.int[2], 3))
}

write.csv(stat, "Fisher_exact_test_position_Immunogenic_Yes_No.csv")

# Length & Mutation Position Analysis
t1 <- read_xlsx(file_path) %>% as.data.frame()
split_data <- split(t1, t1$Peptide_Length)
list_stat <- list()

for (l in 1:4) {
  t1 <- split_data[[l]]
  t2 <- as.data.frame(table(t1[, c("Mutation_Position", "Immunogenic2")]))
  n_yes <- as.data.frame(table(t1$Immunogenic2))[2, 2]
  n_no <- as.data.frame(table(t1$Immunogenic2))[1, 2]
  
  stat <- matrix(NA, nrow = 11, ncol = 4)
  
  for (p in 1:(l + 7)) {
    m1 <- subset(t2, Mutation_Position == p)
    m2 <- matrix(c(m1[2, 3], n_yes - m1[2, 3], m1[1, 3], n_no - m1[1, 3]), ncol = 2, nrow = 2, byrow = TRUE)
    s1 <- fisher.test(m2)
    stat[p, ] <- c(s1$p.value, s1$estimate, s1$conf.int[1], s1$conf.int[2])
  }
  
  row.names(stat) <- paste0("Position_", 1:11)
  colnames(stat) <- c("P.value", "OR", "95%CI_lower", "95%CI_higher")
  list_stat[[l]] <- stat
}

names(list_stat) <- paste0("Length", 8:11)
Stats <- do.call(rbind, list_stat)
Mut_Position <- rep(1:11, 4)
Length <- rep(8:11, each = 11)
Stats <- cbind(Stats, Mutation_Position = Mut_Position, Length = Length)
Stats$OR_100 <- round(as.numeric(Stats$OR), 2) * 100

write.csv(Stats, "07.Fisher_exact_test_mutation_position_length.csv")

# Plot Mutation Position vs. Peptide Length
p1 <- ggplot(Stats, aes(factor(Mutation_Position), factor(Length))) + 
  geom_tile(aes(fill = OR_100), colour = "white") +
  scale_fill_gradient2(low = "#4292C6", high = "#A50F15", midpoint = 100, limit = c(0, 300)) +
  labs(x = "Mutation Position", y = "Peptide Length") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")
  )

ggsave(p1, filename = "07.mutation_position_length.pdf", width = 3.8, height = 3)

# Peptide Length for Strong Binder Enrichment
t1 <- read_xlsx(file_path) %>% as.data.frame()
t2 <- as.data.frame(table(t1[, c("Peptide_Length", "Immunogenic2")]))

p <- ggplot(data = t2, mapping = aes(x = Peptide_Length, y = Freq, fill = Immunogenic2)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), alpha = 0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "", y = "Number of p:MHC") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

ggsave(p, filename = "00.PepLength_Yes_No.pdf", width = 2.6, height = 2.5)

# Peptide Length vs. Strong Binder Ratio
t1 <- read_xlsx(file_path) %>% as.data.frame()
t2 <- as.data.frame(table(na.omit(t1[, c("Peptide_Length", "Immunogenic2")])))
t2_F <- subset(t2, Immunogenic2 == "Yes")
t2_T <- subset(t2, Immunogenic2 == "No")
t2_F$all <- sum(t2_F$Freq)
t2_T$all <- sum(t2_T$Freq)
t3 <- rbind(t2_F, t2_T)
t3$ratio <- t3$Freq / t3$all

p <- ggplot(data = t3, mapping = aes(x = Peptide_Length, y = ratio, fill = Immunogenic2)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), alpha = 0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "Peptide Length", y = "Fraction of All") +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

ggsave(p, filename = "00.peplength_Yes_No_ratio.pdf", width = 2.6, height = 2.5)

# Fisher exact test for peptide length

# Load required libraries
library(readxl)
library(ggplot2)
library(RColorBrewer)

# Define file path
file_path <- "01.Unique_pMHC_DriverMutation_Fusions.xlsx"

# Load and prepare data
t1 <- read_xlsx(file_path) %>% as.data.frame()
t2 <- as.data.frame(table(t1[, c("Peptide_Length", "Immunogenic2")]))

# Fisher exact test for peptide length
options(scipen = 200)
n_yes <- as.data.frame(table(t1$Immunogenic2))[2, 2]
n_no <- as.data.frame(table(t1$Immunogenic2))[1, 2]

stat <- data.frame(matrix(NA, nrow = 4, ncol = 4))
row.names(stat) <- paste0("PepLength_", 8:11)
colnames(stat) <- c("P.value", "OR", "95%CI_lower", "95%CI_higher")

for (p in 1:4) {
  m1 <- subset(t2, Peptide_Length == (p + 7))
  m2 <- matrix(c(m1[2, 3], n_yes - m1[2, 3], m1[1, 3], n_no - m1[1, 3]), ncol = 2, byrow = TRUE)
  fisher_test <- fisher.test(m2)
  stat[p, ] <- c(
    round(fisher_test$p.value, 10),
    round(fisher_test$estimate, 10),
    round(fisher_test$conf.int[1], 10),
    round(fisher_test$conf.int[2], 10)
  )
}

write.csv(stat, "Fisher_exact_test_pepLength_High_Low.csv")

# Variant type for Strong binder enrichment
t1 <- read_xlsx(file_path) %>% as.data.frame()
t2 <- as.data.frame(table(t1[, c("Variant_type", "Immunogenic2")]))
t2$Variant_type <- factor(t2$Variant_type, levels = c("Indel", "SNP", "Fusion"))

# Plot Variant type vs. Immunogenicity
p <- ggplot(data = t2, mapping = aes(x = Variant_type, y = Freq, fill = Immunogenic2)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), alpha = 0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "", y = "Number of p:MHC") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

ggsave(p, filename = "00.VariantType_strong_weak.pdf", width = 2.6, height = 2.5)


# Load required libraries
library(readxl)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggcorrplot)
library(ggrastr) # For rasterized ggplot2 plots

# Set working directory
setwd("D:/01.Projects/V4_AML_mutation_Fusion_Neoantigens/04_PresentationFeature")

# Load and preprocess data
t <- read_xlsx("figure2_plotdata.xlsx") %>% as.data.frame()
t2 <- unique(t[, c("expression_tpm", "Variant_type", "MT_Score", "Stability", "Hydrophobicity_fraction")])
t2$expression_tpm <- log10(t2$expression_tpm + 1)
t2$MT_Score <- log10(t2$MT_Score + 1)

# Subset data based on Variant_type
mat_mut <- subset(t2, Variant_type != "Fusion")
mat_fusion <- subset(t2, Variant_type == "Fusion")
mat_indel <- subset(t2, Variant_type == "Indel")
mat_snv <- subset(t2, Variant_type == "SNP")

# Remove Variant_type column and compute correlation matrices
mat1 <- unique(mat_indel[, -2])
mat2 <- unique(mat_snv[, -2])
mat3 <- unique(mat_fusion[, -2])

corr1 <- round(cor(mat1, method = "pearson"), 3)
corr2 <- round(cor(mat2, method = "pearson"), 3)
corr3 <- round(cor(mat3, method = "pearson"), 3)

# Define function to get lower triangle of correlation matrix
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  cormat[cormat == 1] <- NA
  return(cormat)
}

# Apply function and melt correlation matrices
corr1 <- melt(get_lower_tri(corr1), na.rm = TRUE)
corr2 <- melt(get_lower_tri(corr2), na.rm = TRUE)
corr3 <- melt(get_lower_tri(corr3), na.rm = TRUE)

# Plot correlation matrices
p1 <- ggplot(corr1, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4292C6", high = "#A50F15", mid = "white",
                       midpoint = 0, limits = c(-0.5, 0.5)) +
  theme_minimal() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  labs(x = "", y = "")

p2 <- ggplot(corr2, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4292C6", high = "#A50F15", mid = "white",
                       midpoint = 0, limits = c(-0.5, 0.5)) +
  theme_minimal() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  labs(x = "", y = "")

p3 <- ggplot(corr3, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4292C6", high = "#A50F15", mid = "white",
                       midpoint = 0, limits = c(-0.5, 0.5)) +
  theme_minimal() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  labs(x = "", y = "")

# Combine and save plots
plots <- p1 + p2 + p3
ggsave(plots, filename = "08.Mutation_fusion_corr.pdf", width = 9.6, height = 3.6)

# Binding affinity vs. stability
t <- read_xlsx("figure2_plotdata.xlsx") %>% as.data.frame()

col <- c("#5AAE61", "#9970AB", "black")
p1 <- ggplot(data = t, aes(MT_Score, Stability, color = Immunogenic)) +
  geom_point_rast(size = 3, alpha = 0.5, shape = 16, stroke = 1) +
  labs(x = "Binding Affinity", y = "Binding Stability") +
  scale_color_manual(values = col) +
  theme_bw() +
  geom_smooth(method = "lm", se = TRUE, color = col[3]) +
  scale_x_log10() +
  scale_y_log10()

stat <- cor.test(log10(t$MT_Score), log10(t$Stability), method = "pearson")
anno <- c(stat$estimate, stat$p.value)
names(anno) <- c("R", "P")
lab <- paste(names(anno), "=", round(anno, 6), collapse = "\n")

p2 <- p1 + annotate("text", 10, 15, label = lab)
ggsave(p2, filename = "09.corr_Binding_Affinity_Stability.pdf", width = 7, height = 5)

# Bar plot for three variant types
setwd("D:/01.Projects/V4_AML_mutation_Fusion_Neoantigens/04_PresentationFeature")
t <- read_xlsx("01.Unique_pMHC_DriverMutation_Fusions.xlsx") %>% as.data.frame()
t1 <- t[, c("Variant_type", "Presented", "Immunogenic2")]

t1_pre <- subset(t1, Presented == "Yes")
t1_pre_imm <- subset(t1_pre, Immunogenic2 == "Yes")

tab1 <- data.frame(table(t1_pre$Variant_type))
tab2 <- data.frame(table(t1_pre_imm$Variant_type))

tab1$type <- "Presented"
tab2$type <- "Presented_Immu"

tab3 <- rbind(tab1, tab2)

source("D:/03.R_functions/JPComUse/R/PlotBar.R")
p1 <- PlotBar(data = tab3, x.column = "type", y.column = "Freq",
              col.palette = "Set2", fill.column = "Var1", coord_flip = FALSE)

# Save the plot
ggsave(p1, filename = "10.VariantType_barplot.pdf", width = 3.4, height = 3)

