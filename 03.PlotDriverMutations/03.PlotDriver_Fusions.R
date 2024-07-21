## 2021.10.27

# Set working directory
setwd("D:\\Documents\\Desktop\\V3_AML_DriverNeo\\04.Neoantigen.Pred\\PvacFuse")

# Load required libraries
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ggrastr)

# Read and process data for the box plot
fusion_data <- read.csv("Plot_fusions_hg38.csv", header = TRUE)
fusion_counts <- data.frame(table(fusion_data$Fusion_ID))

# Set factor levels for Fusion
fusion_levels <- c("KMT2A--MLLT1", "KMT2A--MLLT10", "KMT2A--ELL", "KMT2A--MLLT3", "PML--RARA", "KMT2A--MLLT4", "BCR--ABL1", "CBFB--MYH11", "RUNX1--RUNX1T1")
fusion_data$Fusion <- factor(fusion_data$Fusion, levels = fusion_levels)

# Create and save the box plot
p1 <- ggplot(data = fusion_data, aes(Fusion, FFPM)) +
  stat_boxplot(geom = "errorbar", width = 0.25, aes(color = "black")) + 
  geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") + 
  geom_jitter(aes(fill = Fusion), width = 0.2, shape = 21, size = 2, alpha = 0.6) + 
  scale_fill_manual(values = brewer.pal(9, "Set1")) +  
  scale_color_manual(values = c("black", "black")) + 
  theme_classic() + 
  ylab("Fusion fragments per million total reads (FFPM)") + 
  xlab("") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed", color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = "dashed", color = "grey")
  ) +
  scale_y_log10()

ggsave(p1, filename = "Fusion_expression.pdf", width = 5, height = 4)

# Load additional libraries for the heatmap
library(circlize)
library(ComplexHeatmap)

# Read and process data for the heatmap
MutMat <- read.table("Plot_fusions_Landscape.txt", header = TRUE)
source("C:/JP_R.packages/JPComUse/R/Create.matrix.R")
Mutmat_wide <- create.matrix(MutMat, Sample_ID = "Tumor_Sample_Barcode", GeneSymbol = "Hugo_Symbol", column = "Variant_Classification")
Mutmat_wide[Mutmat_wide == 0] <- ''

# Define colors and alteration functions for the heatmap
mycol <- c(pal_nejm("default")(8), "#8DD3C7", "#CCEBC5")

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.8, "mm"), gp = gpar(fill = '#ebebe1ff', col = NA))
  },
  inversion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.8, "mm"), gp = gpar(fill = mycol[2], col = NA))
  },
  translocation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0, "mm"), h - unit(0.8, "mm"), gp = gpar(fill = mycol[1], col = NA))
  }
)

col <- c("inversion" = mycol[2], "translocation" = mycol[1])

# Import and process clinical data
clinic <- read.csv("Plot_fusions_dataset.csv", header = TRUE)
pdata <- subset(clinic, clinic$Tumor_Sample_Barcode %in% colnames(Mutmat_wide))
Mutmat_wide <- Mutmat_wide[, pdata$Tumor_Sample_Barcode]

# Define the clinical annotations for the heatmap
ClinicAnno <- HeatmapAnnotation(
  DataSet = pdata$DataSet,
  show_annotation_name = TRUE,
  col = list(DataSet = c("Leucegene" = "#A6D854", "RJAPL" = "#FC8D62", "BeatAML" = "#8DA0CB", "RJAML" = "#E78AC3", "TCGA" = "#66C2A5")),
  annotation_name_gp = gpar(fontsize = 7)
)

# Create and save the heatmap
p1 <- oncoPrint(
  Mutmat_wide, col = col, 
  top_annotation = NULL,
  right_annotation = NULL,
  alter_fun = alter_fun,
  remove_empty_columns = TRUE,
  alter_fun_is_vectorized = FALSE,
  row_names_side = "left", 
  pct_side = "",
  bottom_annotation = ClinicAnno
)

pdf("AML_fusionMutLand.pdf", width = 9, height = 4)
p1
dev.off()

