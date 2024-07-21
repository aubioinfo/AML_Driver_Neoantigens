##########
########## immunologic milieu within the AML micro-environment 
##########

# Load required libraries
library(readxl)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggstatsplot)
library(patchwork)
library(ggrastr)

# Set working directory and load data
setwd("D:/01.Projects/V4_AML_mutation_Fusion_Neoantigens/06_clinicalData/ImmuneCells")
data <- read_xlsx("../03_Merged_Cohorts_Neo_Clinical.xlsx")

# Read and preprocess count data
count <- read.table("TCGA179_BeatAML477_RJAML351.count.txt", header = TRUE, row.names = 1)
count <- as.matrix(count)
batch <- c(rep(1, 477), rep(2, 351), rep(3, 179))
adjusted <- ComBat_seq(count, batch = batch, group = NULL)
write.csv(adjusted, "BeatAML477_RJAML351_TCGA179_Adjusted_count.csv")

# Read and preprocess TPM data
tpm <- read.csv("BeatAML477_RJAML351_TCGA179_Adjusted_TPM.csv", header = TRUE, row.names = 1)
markers <- read.csv("TcellGenes.csv", header = TRUE)

# Merge marker data with TPM data
marker_mat <- merge(markers, tpm, by = "Symbol")
marker_mat1 <- marker_mat[, c("Symbol", "Type", as.character(data$sample))]
write.csv(marker_mat1, "Marker_Three_Cohorts.tpm.csv")

# Process and aggregate TPM data
matched_tpm <- subset(tpm, Type == "protein_coding")
matched_tpm_agg <- aggregate(matched_tpm[, -c(1:2)], by = list(matched_tpm$Symbol), FUN = mean)
write.csv(matched_tpm_agg, "Three_matched_Adjusted_TPM.csv")

# Plotting
# Load additional libraries
library(ggExtra)

# Read marker scores data
t1 <- read.csv("Top25%_bottom25%_Groups_MarkerScores_Three_Cohorts.csv", header = TRUE)
t1$group <- factor(t1$group, levels = c("Yes", "No"))

# Define plot function
create_violin_plot <- function(data, y_var, y_label, label_y) {
  ggplot(data, aes(x = group, y = log2(get(y_var)), fill = factor(group))) +
    scale_fill_brewer(palette = "Pastel1") +
    geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA) +
    geom_point_rast(aes(fill = factor(group), alpha = 0.4), size = 0.2, position = position_jitterdodge(0.75)) +
    theme_bw() +
    ylab(y_label) +
    xlab("") +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
          panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
          panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")) +
    stat_compare_means(mapping = aes(x = group, y = get(y_var)), label.y = label_y)
}

# Create and save plots
p1 <- create_violin_plot(t1, "Cytolytic_score", "Cytolytic score", 7.5)
p2 <- create_violin_plot(t1, "CTL", "Cytotoxic T lymphocyte score (TIDE)", 3)
p3 <- create_violin_plot(t1, "Dysfunction", "Dysfunction score (TIDE)", 2.5)

plot1 <- p1 + p2 + p3
ggsave("01_Genes.pdf", plot = plot1, width = 5, height = 3)

# Define correlation plot function
create_correlation_plot <- function(data, group_col, color_val) {
  ggplot(data, aes(CTL, Dysfunction, color = group)) + 
    geom_point_rast(size = 1, alpha = 1, shape = 16, stroke = 1) + 
    labs(x = "Cytotoxic T lymphocyte score", y = "Dysfunction score") +
    scale_color_manual(values = color_val) + theme_bw() + 
    geom_smooth(method = "lm", se = TRUE, color = "black")
}

# Correlation analysis and plotting
correlation_and_plot <- function(data, group_name, color_val) {
  stat <- cor.test(data$CTL, data$Dysfunction, method = "pearson")
  anno <- c(stat$estimate, stat$p.value)
  names(anno) <- c("R", "P")
  lab <- paste(names(anno), "=", round(anno, 6), collapse = "\n")
  
  p <- create_correlation_plot(data, group_name, color_val) +
    annotate("text", 0, 2, label = lab)
  return(p)
}

p1 <- correlation_and_plot(t1, "Yes", brewer.pal(3, "Pastel1")[1])
p2 <- correlation_and_plot(subset(t1, group == "Yes"), "Yes", brewer.pal(3, "Pastel1")[1])
p3 <- correlation_and_plot(subset(t1, group == "No"), "No", brewer.pal(3, "Pastel1")[2])

plot2 <- p1 + p2 + p3
ggsave("02_Cor.pdf", plot = plot2, width = 9.6, height = 3)


####################################
###### DE expression analysis ######
####################################
# Load necessary libraries
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\06_clinicalData\\ImmuneCells")

# Source custom functions
source("D:\\03.R_functions\\JPComUse\\R\\RunDESeq2.R")
source("D:\\03.R_functions\\JPComUse\\R\\plotVolcanoV3.R")
source("D:\\03.R_functions\\JPComUse\\R\\RunGSEA.R")

# Read TPM and count matrix
tpm <- read.csv("Top25%_bottom25%_TPM_matrix.csv", header = TRUE)
mat <- read.csv("Top25%_bottom25%_Count_matrix.csv", header = TRUE, row.names = 1)

# Filter and reorder count matrix
mat <- mat[rowMeans(mat) > 0, ]
mat <- mat[, c(172:342, 1:171)]

# Run DESeq2 analysis
RunDESeq2(count_mat = mat, n.cont = 171, n.treat = 171, 
          prefix = "DESeq2_", sort.p = FALSE, 
          merge.normalized = TRUE, normalized_mat = tpm)

###########################
### Volcano plot #########
###########################
# Read DESeq2 results
m1 <- read.table("DESeq2__with_normalized_mat.txt", header = TRUE)

# Define genes of interest
gene_selecte <- c("STAT1", "IFIT1", "CXCL10", "MX1", "IRF1", "CD8A", "CD8B", "PRF1", "GZMB", "CD274", "LAG3", "CTLA4", "IDO1")

# Create volcano plot
p <- plotVolcano(mat = m1, gene.col = "Symbol", x.col = "log2FoldChange", y.col = "pvalue", 
                 x_cut1 = 0.3, x_cut2 = 0.5, y_cut1 = 0.05, y_cut2 = 0.01, x.lim = 5, y.lim = 15, 
                 label = TRUE, selected_genes = gene_selecte,
                 labx = bquote(~Log[2]~"(fold change)"), 
                 laby = bquote(~-Log[10]~italic("FDR")),
                 title = "High vs. Low")
ggsave(p, filename = "03_Valcanoplot.pdf", width = 5, height = 5)

###########################
### GSEA plot ############
###########################
# Prepare gene list for GSEA
DEgene_mat <- read.table("DESeq2__with_normalized_mat.txt", header = TRUE, sep = "\t")
DEgene_mat1 <- DEgene_mat[rowMeans(DEgene_mat[, 10:351]) > 0.1, ]
DE_list <- DEgene_mat1[, c(9, 3)]
colnames(DE_list)[1] <- "ID"
DE_list <- distinct(DE_list, ID, .keep_all = TRUE) %>% na.omit()
geneList <- sort(setNames(DE_list$log2FoldChange, DE_list$ID), decreasing = TRUE)

# Load gene sets and run GSEA
geneSet <- read.gmt("h.all.v7.5.1.symbols.gmt")
set.seed(123456)
gsea.enrich <- GSEA(geneList, TERM2GENE = geneSet, pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE)
a <- gsea.enrich@result

# Plot GSEA results
p <- gseaplot2(gsea.enrich, c(3:4, 9), color = brewer.pal(8, "Set2")[1:3],
                title = "",
                subplots = c(1:3),
                rel_heights = c(0.75, 0.2, 0.3),
                pvalue_table = FALSE)
ggsave(p, filename = "04_GSEA.pdf", width = 6, height = 5)

# Save GSEA results
write.csv(a, "04.GSEA.out.csv")
