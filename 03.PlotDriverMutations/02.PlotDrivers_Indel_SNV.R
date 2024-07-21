### 2021/8/3
### Jin Peng
### Mutation landscape of AML (Five Cohorts)

# Load required packages
library(ggplot2)
library(ggrepel)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(readxl)
library(reshape2)
library(ggsignif)
library(TCGAmutations)
library(ggpubr)
library(tidyr)

# Set working directory
setwd("D:\\Documents\\Desktop\\V3_AML_DriverNeo\\03.DriverGeneMutations")

# Load mutation data
driver_genes <- read.table("49_driverMut.txt", header = TRUE)
mutations <- read.table("FiveCohortsCombined_mutations.txt", header = TRUE)
mutations <- merge(mutations, driver_genes, by = "Hugo_Symbol")

# Create unique mutation matrix and calculate mutation frequencies
unique_mutations <- unique(mutations[, c(1, 3)])
mutation_frequencies <- data.frame(table(unique_mutations$Hugo_Symbol))
mutation_frequencies <- mutation_frequencies[order(mutation_frequencies$Freq, decreasing = TRUE), ]

# Load custom functions
source("C:/JP_R.packages/JPComUse/R/Create.matrix.R")

# Create wide-format mutation matrix
mutation_matrix <- create.matrix(mutations, 
Sample_ID = "Tumor_Sample_Barcode", GeneSymbol = "Hugo_Symbol", column = "Variant_Classification")
mutation_matrix[mutation_matrix == 0] <- ''

# Define colors for variant types
variant_colors <- c(pal_nejm("default")(8), "#8DD3C7", "#CCEBC5", "#FC8D62", "#8DA0CB")

# Define functions for drawing mutation types
alter_functions <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = '#ebebe1ff', col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[2], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[1], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[3], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[4], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[5], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[6], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[7], col = NA))
  },
  Splice_Region = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[8], col = NA))
  },
  Intron = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[9], col = NA))
  },
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[11], col = NA))
  },
  Silent = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[10], col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = variant_colors[12], col = NA))
  }
)

# Define variant type colors
variant_colormap <- c(
  "Missense_Mutation" = variant_colors[2], "Frame_Shift_Ins" = variant_colors[1], 
  "Nonsense_Mutation" = variant_colors[3], "Frame_Shift_Del" = variant_colors[4], 
  "In_Frame_Ins" = variant_colors[5], "Splice_Site" = variant_colors[6], 
  "In_Frame_Del" = variant_colors[7], "Splice_Region" = variant_colors[8], 
  "Intron" = variant_colors[9], "Translation_Start_Site" = variant_colors[11], 
  "Silent" = variant_colors[10], "Nonstop_Mutation" = variant_colors[12]
)

# Load clinical data
clinical_data <- read.csv("CombinedAML_clinicalData.csv", header = TRUE)
pdata <- subset(clinical_data, clinical_data$Tumor_Sample_Barcode %in% colnames(mutation_matrix))
mutation_matrix <- mutation_matrix[, pdata$Tumor_Sample_Barcode]

# Define clinical annotations
clinical_annotations <- HeatmapAnnotation(
  Platform = pdata$Platform,
  Center = pdata$Center,
  show_annotation_name = TRUE,
  col = list(
    Platform = c("WES" = "#A6D854", "Panel" = "#FC8D62"),
    Center = c("RJ_APL" = "#D95F02", "RJ_AML" = "#7570B3", "TCGA_LAML" = "#E7298A", "BeatAML" = "#66A61E", "AMLSG" = "#E6AB02")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

# Plot the heatmap
heatmap <- oncoPrint(
  mutation_matrix, col = variant_colormap, alter_fun = alter_functions,
  remove_empty_columns = FALSE, alter_fun_is_vectorized = FALSE,
  column_title = "Driver mutation landscape of AML",
  row_names_side = "left", pct_side = "right", bottom_annotation = clinical_annotations
)

# Save the heatmap plot
pdf("01.AML.DriverMutLand.pdf", width = 18, height = 9.5)
draw(heatmap)
dev.off()

# Load driver mutation data
mutation_data <- read_xlsx("02.All_Driver_mut_type.xlsx")[, 2:3]
mutation_table <- table(mutation_data)
mutation_data_melted <- melt(mutation_table)
mutation_data_melted$Variant_Type <- factor(mutation_data_melted$Variant_Type, levels = c("All", "Driver"))
mutation_data_melted$Type <- factor(mutation_data_melted$Type, levels = c("SNV", "Indel"))

# Load custom plotting functions
source("C:\\JP_R.packages\\JPComUse\\R\\PlotBar.R")

# Plot driver mutation types
driver_mutation_plot <- PlotBar(
  data = mutation_data_melted, x.column = "Variant_Type", y.column = "value",
  col.palette = "Set2", fill.column = "Type", coord_flip = FALSE
)

# Save the driver mutation plot
ggsave(driver_mutation_plot, filename = "02.All_Driver_mut_type.pdf", width = 2, height = 3)

# Fisher's exact test
fisher_table <- matrix(c(mutation_table[2, 1], mutation_table[1, 1], mutation_table[2, 2], mutation_table[1, 2]), ncol = 2, nrow = 2)
rownames(fisher_table) <- c("Driver", "All")
colnames(fisher_table) <- c("Indel", "SNV")
fisher_test_result <- fisher.test(fisher_table)


# Calculate Indel/(Indel + SNV) for pan-cancer (Variant level)

# Load required libraries
library(TCGAmutations)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# Set working directory
setwd("D:/Documents/Desktop/V3_AML_DriverNeo/03.DriverGeneMutations")

# Load TCGA mutation data for all available studies except LAML
mat_ava <- tcga_available()
tumor_type <- mat_ava$Study_Abbreviation[-34]
tumor_IndelPro_varian_out <- list()

for (tuty in tumor_type) {
  cancer <- tcga_load(study = tuty, source = "MC3")
  maf <- cancer@data
  maf <- maf %>% filter(Variant_Type %in% c("SNP", "DEL", "INS"))
  maf <- maf %>% mutate(Var_Type = ifelse(Variant_Type == "SNP", "SNV", "Indel"))
  maf <- maf %>% unite("Mut_ID", colnames(maf)[1:6], sep = "|", remove = TRUE) %>% select(Mut_ID, Var_Type) %>% distinct()
  maf <- maf %>% mutate(Tumor_type = tuty)
  tumor_IndelPro_varian_out[[tuty]] <- maf
}

pancan_variants <- bind_rows(tumor_IndelPro_varian_out)

# Read the AML WES data and combine with pan-cancer data
aml_wes <- read.table("03.AML_WES_all_driver_unique_variants.txt", header = TRUE)
pancan_variants <- bind_rows(pancan_variants, aml_wes)

# Prepare data for plotting
pancan_variants <- pancan_variants %>% select(Tumor_type, Var_Type)
data1 <- table(pancan_variants)
data2 <- melt(data1)
data2 <- data2 %>% filter(Tumor_type != "LAML") %>% mutate(Var_Type = factor(Var_Type, levels = c("SNV", "Indel")))

# Plot Indel proportions across tumor types
source("C:/JP_R.packages/JPComUse/R/PlotBar.R")
p <- PlotBar(data = data2, x.column = "Tumor_type", y.column = "value", col.palette = "Set2", fill.column = "Var_Type", coord_flip = FALSE)

# Save the plot
ggsave(p, filename = "03.Pancancer_mut_type.pdf", width = 6, height = 3.5)

# Calculate Indel proportion for all mutations
vrt_all <- read.table("04_1.All_variants.txt", header = TRUE)
va_2 <- vrt_all %>% select(Tumor_Sample_Barcode, Type)
va_3 <- table(va_2)
va_4 <- subset(va_3, Type == "Indel")
colnames(va_4)[3] <- "Indel_count"
va_5 <- aggregate(Freq ~ Tumor_Sample_Barcode, data = va_3, FUN = sum)
va_6 <- cbind(va_4, va_5)
va_6 <- va_6 %>% mutate(P_indels = Indel_count / Freq)

# Calculate Indel proportion for driver mutations
vrt_driver <- read.table("04_2.Driver_variants.txt", header = TRUE)
vd_2 <- vrt_driver %>% select(Tumor_Sample_Barcode, Type)
vd_3 <- table(vd_2)
vd_4 <- subset(vd_3, Type == "Indel")
colnames(vd_4)[3] <- "Indel_count"
vd_5 <- aggregate(Freq ~ Tumor_Sample_Barcode, data = vd_3, FUN = sum)
vd_6 <- cbind(vd_4, vd_5)
vd_6 <- vd_6 %>% mutate(P_indels = Indel_count / Freq)

# Add mutation type and prepare data for plotting
va_6 <- va_6 %>% mutate(Tumor_type = "All")
vd_6 <- vd_6 %>% mutate(Tumor_type = "Driver")
vad <- bind_rows(va_6, vd_6) %>% select(-Type, -Freq)

# Plot the proportion of indels
p1 <- ggplot(data = vad) + 
  geom_boxplot(aes(x = Tumor_type, y = P_indels, fill = Tumor_type), outlier.colour = NA, alpha = 0.6) +
  theme_classic() + 
  labs(x = "", y = "Proportion of indels") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1, color = "black"),
        panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey"))

# Save the plot
ggsave(p1, filename = "04.AML_Proportion_indels_wes_panel.pdf", width = 2, height = 3)


# Pancancer "Proportion of indels" analysis
# Load required libraries
library(TCGAmutations)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(RColorBrewer)

# Fetch available tumor types
mat_ava <- tcga_available()
tumor_types <- as.vector(mat_ava$Study_Abbreviation)
tumor_types <- tumor_types[-34]  # Remove LAML
tumor_PI_out <- list()

# Process each tumor type
for (tumor in tumor_types) {
  cancer <- tcga_load(study = tumor, source = "MC3")
  maf <- cancer@data
  
  # Count variants by sample and type
  count <- data.frame(table(maf[, c(8, 16)]))
  count <- subset(count, Variant_Type %in% c("SNP", "INS", "DEL"))
  
  all_count <- aggregate(Freq ~ Tumor_Sample_Barcode, data = count, FUN = sum)
  indel <- subset(count, Variant_Type %in% c("INS", "DEL"))
  indel_agg <- aggregate(Freq ~ Tumor_Sample_Barcode, data = indel, FUN = sum)
  colnames(indel_agg)[2] <- "Indel_count"
  
  # Merge counts and calculate indel proportion
  merged_mat <- merge(indel_agg, all_count, by = "Tumor_Sample_Barcode")
  merged_mat$P_indels <- merged_mat$Indel_count / merged_mat$Freq
  merged_mat$Tumor_type <- tumor
  
  tumor_PI_out[[tumor]] <- merged_mat
}

# Combine results
merged <- do.call(rbind, tumor_PI_out)
merged <- merged[, -3]
merged <- subset(merged, Tumor_type != "LAML")

# Read AML WES sample IDs and merge with data
wes_aml_id <- read.csv("WES_AML_sampleID.csv", header = TRUE)
merge_aml <- rbind(merged, vad)

# Plot settings
mycolors <- unique(c(
  pal_nejm("default")(8), pal_npg("nrc")(8), brewer.pal(8, "Set1"), 
  brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), "#8DD3C7", "#CCEBC5", 
  "#FC8D62", "#8DA0CB"
))

# Order tumor types by median indel proportion
order_type <- aggregate(P_indels ~ Tumor_type, data = merge_aml, FUN = median)
order_type <- order_type[order(order_type$P_indels, decreasing = FALSE),]
merge_aml$Tumor_type <- factor(merge_aml$Tumor_type, levels = order_type$Tumor_type)

# Create boxplot
p2 <- ggplot(data = merge_aml) +
  geom_boxplot(aes(x = Tumor_type, y = P_indels, fill = factor(Tumor_type)), 
               outlier.colour = NA, alpha = 0.6) +
  theme_classic() + 
  labs(x = "", y = "Proportion of Indels") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  )

# Adjust y-axis limit and add color palette and statistical comparison
p3 <- p2 + 
  scale_y_continuous(limits = c(0, 0.6)) + 
  scale_fill_manual(values = mycolors) + 
  stat_compare_means(mapping = aes(x = Tumor_type, y = P_indels))

# Save the plot
ggsave(p3, filename = "05.Prop_indels_Pancancer.pdf", width = 7, height = 3.5)


### Plot the Fusion landscape of rjAML cohort
# Load required libraries
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(readxl)

# Set working directory
setwd("D:\\Documents\\Desktop\\Driver_Neos\\03_plot")

# Read and process data for the heatmap
mat1 <- read.table("02.driver_Fusion_plot.txt", header = TRUE)
source("C:\\JP_R.packages\\JPComUse\\R\\Create.matrix.R")
mut_matrix1 <- create.matrix(mat1, Sample_ID = "Tumor_Sample_Barcode", GeneSymbol = "Hugo_Symbol", column = "Variant_Classification")
mut_matrix1[mut_matrix1 == 0] <- ''

# Define colors and alteration functions
mycol <- c("#2D7597FF", "#C37469FF", pal_jco("default")(10))
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#ebebe1ff', col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycol[1], col = NA))
  }
)

# Read clinical data and create annotation
aml_clin <- read_xlsx("02.driver_Fusion_clinical.xlsx")
clini <- HeatmapAnnotation(Center = aml_clin$Center, show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 7))

# Create and save heatmap
p <- oncoPrint(
  t(mut_matrix1), col = mycol, alter_fun = alter_fun,
  top_annotation = NULL, alter_fun_is_vectorized = FALSE, 
  bottom_annotation = clini
)

pdf("02_Fusion_landscape.pdf", width = 10, height = 4)
p
dev.off()

# Read and process data for the bar plot
t1 <- read.table("03_driver_variants_plot.txt", header = TRUE)
t2 <- data.frame(table(t1$Variant_Type))

# Create and save bar plot
p <- ggplot(t2, aes(Var1, Freq, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.65) +
  labs(x = '', y = 'No. driver variants') +
  theme_classic() +
  theme(
    axis.text.x = element_text(hjust = 1, color = "black"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave(p, filename = "03_Num_driver_variants.pdf", width = 2.5, height = 2.8)
