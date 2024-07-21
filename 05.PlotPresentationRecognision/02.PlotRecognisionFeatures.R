###### Neoantigen presentation features
###### Jin Peng 
###### 2022.1.12

## Recognition features
options(scipen = 10)

# Load required libraries
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggcorrplot)
library(reshape2)

# Set working directory
setwd("D:/01.Projects/V4_AML_mutation_Fusion_Neoantigens/05_RecognitionFeatures")

# Load and preprocess data
t1 <- read_xlsx("Presented_pMHC.xlsx") %>% as.data.frame()

# Plot Agretopicity histogram
p1 <- ggplot(t1, aes(x = Agretopicity, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 10) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_log10() +
  labs(title = "Agretopicity") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")
  )

# Plot Foreignness histogram
p2 <- ggplot(t1, aes(x = Foreignness, fill = Variant_type)) +
  geom_histogram(alpha = 0.5, bins = 15) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_log10() +
  labs(title = "Foreignness") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")
  )

# Combine and save histograms
plots <- p1 + p2
ggsave(plots, filename = "00.RecognitionFeatures.pdf", width = 5.5, height = 3)

# Load data for correlation plot
t1 <- read_xlsx("Mutation_Presented_pMHC.xlsx") %>% as.data.frame()
t2 <- t1[, c("MT_Score", "Stability", "expression_tpm", "Foreignness", "Agretopicity")]
t2$expression_tpm <- log10(t2$expression_tpm + 1)

# Compute correlation matrix
corr <- round(cor(t2, method = "pearson"), 5)

# Plot correlation matrix
p <- ggcorrplot(
  corr, hc.order = FALSE,
  outline.col = "white",
  ggtheme = ggplot2::theme_void,
  colors = c("#6D9EC1", "white", "#E46726")
)
ggsave(p, filename = "01.corrplot.pdf", width = 6, height = 4.6)

# Melt and plot detailed correlation matrix
corr1 <- melt(corr, na.rm = TRUE)
p1 <- ggplot(corr1, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#6D9EC1", high = "#E46726", mid = "white",
    midpoint = 0, limits = c(-0.3, 0.3)
  ) +
  theme_minimal() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  labs(x = "", y = "")
ggsave(p1, filename = "01.corrplot_detailed.pdf", width = 6, height = 4.6)

# Plot Foreignness vs. Agretopicity
col <- c("#5AAE61", "#9970AB", "black")
p2 <- ggplot(data = t1, aes(Foreignness, Agretopicity, color = Immunogenic2)) +
  geom_point(size = 3, alpha = 0.5, shape = 16) +
  scale_color_manual(values = col) +
  labs(x = "Foreignness", y = "Agretopicity") +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()
ggsave(p2, filename = "02.Agreto_Fore.pdf", width = 5, height = 3.6)

# Fisher test for recognition
t_recog <- subset(t1, Recognized == "Yes")
t_not <- subset(t1, Recognized == "No")

recog_mat <- data.frame(table(t_recog$Immunogenic2))
notrecog_mat <- data.frame(table(t_not$Immunogenic2))

fishermat <- matrix(c(
  recog_mat[2, 2], recog_mat[1, 2],
  notrecog_mat[2, 2], notrecog_mat[1, 2]
), ncol = 2, nrow = 2, byrow = FALSE)
row.names(fishermat) <- c("Yes", "No")
colnames(fishermat) <- c("Recognized", "Not_recognized")
s1 <- fisher.test(fishermat)
print(s1)

# Bar plot for variant types
t <- read.csv("Recognized_VariantType.csv", header = TRUE)
t1_rec_imm <- subset(t, Immunogenic2 == "Yes")

tab1 <- data.frame(table(t$Variant_type))
tab2 <- data.frame(table(t1_rec_imm$Variant_type))

tab1$type <- "Recognized"
tab2$type <- "Recognized_Immu"

tab3 <- rbind(tab1, tab2)

source("D:/03.R_functions/JPComUse/R/PlotBar.R")
p3 <- PlotBar(data = tab3, x.column = "type", y.column = "Freq",
              col.palette = "Set2", fill.column = "Var1", coord_flip = FALSE)
ggsave(p3, filename = "03.VariantType_barplot.pdf", width = 3.1, height = 3)

# Fisher test for presented and recognized (Mutation)
fishermat_mut <- matrix(c(
  89, 792,
  36, 2906
), ncol = 2, nrow = 2, byrow = TRUE)
row.names(fishermat_mut) <- c("Yes", "No")
colnames(fishermat_mut) <- c("Presented_Recognized", "Not_Presented_Recognized")
s1_mut <- fisher.test(fishermat_mut)
print(s1_mut)

# Fisher test for presented and recognized (Fusion)
fishermat_fusion <- matrix(c(
  5, 29,
  2, 163
), ncol = 2, nrow = 2, byrow = TRUE)
row.names(fishermat_fusion) <- c("Yes", "No")
colnames(fishermat_fusion) <- c("Presented_Recognized", "Not_Presented_Recognized")
s1_fusion <- fisher.test(fishermat_fusion)
print(s1_fusion)


###### PRROC Analysis

# Load required libraries
library(PRROC)
library(dplyr)

# Set working directory
setwd("D:/01.Projects/V4_AML_mutation_Fusion_Neoantigens/05_RecognitionFeatures")

# Function to calculate AUPRC
#' @param ranked.list Data frame with columns:
#'   - hla: HLA allele
#'   - rank: Rank of the pMHC
#'   - sequence: Protein sequence of the putative neoepitope
#' @param tested.list Data frame with columns:
#'   - hla: HLA allele
#'   - sequence: Protein sequence of the tested peptide
#'   - validated: Boolean indicating if the peptide is validated (TRUE) or not (FALSE)
#' @return AUPRC score
auprc_calculation <- function(ranked.list, tested.list) {
  # Merge the lists
  merged.list <- inner_join(ranked.list, tested.list, by = c('hla', 'sequence'))
  
  # Check if the merged list is not empty
  if (nrow(merged.list) > 0) {
    validated.ranks <- merged.list %>% filter(validated == TRUE) %>% pull(rank)
    not.validated.ranks <- merged.list %>% filter(validated == FALSE) %>% pull(rank)
    
    # Calculate AUPRC
    auprc <- PRROC::pr.curve(-validated.ranks, -not.validated.ranks, curve = TRUE)
    auprc.score <- auprc$auc.integral
    
    return(auprc.score)
  } else {
    return(NA)
  }
}

# Define a helper function to process rank lists and save plots
process_ranklist <- function(ranked_file, tested_file, output_pdf) {
  ranked.list <- read.table(ranked_file, header = TRUE)
  tested.list <- read.table(tested_file, header = TRUE)
  
  # Calculate AUPRC
  auprc_score <- auprc_calculation(ranked.list, tested.list)
  
  # Print AUPRC score
  print(paste("AUPRC Score for", ranked_file, ":", auprc_score))
  
  # Plot AUPRC
  merged.list <- inner_join(ranked.list, tested.list, by = c('hla', 'sequence'))
  validated.ranks <- merged.list %>% filter(validated == TRUE) %>% pull(rank)
  not.validated.ranks <- merged.list %>% filter(validated == FALSE) %>% pull(rank)
  auprc <- PRROC::pr.curve(-validated.ranks, -not.validated.ranks, curve = TRUE)
  
  pdf(output_pdf, width = 5, height = 4.5)
  plot(auprc)
  dev.off()
}

# Process Mutation-derived rank lists
process_ranklist("BindingAffinity_ranklist.txt", "testlist.txt", "BindingAffinity_prc.pdf")
process_ranklist("Prsented_ranklist.txt", "testlist.txt", "Prsented_ranklist_prc.pdf")
process_ranklist("Prsented_Recognized_ranklist.txt", "testlist.txt", "Prsented_Recognized_ranklist_prc.pdf")

# Process Fusion-derived rank lists
process_ranklist("BindingAffinity_ranklist_fusion.txt", "testlist_fusion.txt", "BindingAffinity_prc_fusion.pdf")
process_ranklist("Prsented_ranklist_fusion.txt", "testlist_fusion.txt", "Prsented_ranklist_prc_fusion.pdf")
process_ranklist("Prsented_Recognized_ranklist_fusion.txt", "testlist_fusion.txt", "Prsented_Recognized_ranklist_prc_fusion.pdf")
