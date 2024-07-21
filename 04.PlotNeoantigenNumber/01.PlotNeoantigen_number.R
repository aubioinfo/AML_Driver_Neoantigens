###### Count the frequency of HLA in 1028 AML patients
###### 2021.12.14
###### Jin Peng

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\03_HLA")

# Load necessary libraries
library(reshape2)
library(gg.gap)
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(ggpmisc)
library(ggrepel)
library(RColorBrewer)
library(ggsci)
library(ggsignif)
library(ggpubr)

# Read and process HLA data
hla_data <- read.table("1028_Matched_DNA_RNA_HLAs.txt", header = TRUE)
hla_data_long <- melt(hla_data, id.vars = "Sample")
hla_data_unique <- unique(hla_data_long[, -2])
hla_count <- data.frame(table(hla_data_unique$value))
hla_count$Ratio <- hla_count$Freq / 1028
write.csv(hla_count, "01.192_HLAs_1028Patients.csv")

# Plot top 15 HLA frequency
top15_hla <- read.table("top15_hla_freq_plot.txt", header = TRUE)
top15_hla$HLA <- factor(top15_hla$HLA, levels = top15_hla$HLA)

p1 <- ggplot(top15_hla, aes(HLA, Ratio, fill = Type)) + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) + 
  labs(x = '', y = 'Allele Frequency') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")) +
  scale_fill_brewer(palette = "Accent")

ggsave(p1, filename = "01.Top15_HLA_allele.pdf", width = 8, height = 3)

# Process pMHC numbers per mutation/fusion
mutation_data <- read_xlsx("Matched_DNA_RNA_mutation_filtered_epitopes.xlsx")
mutation_data <- subset(mutation_data, Mut_type == "Driver" & Variant_type != "DNP")
unique_pmhc <- unique(mutation_data[, c("Gene_Name", "Variant_type", "ID.1")])
unique_pmhc <- unique_pmhc[, c(1, 2)]
gene_pmhc <- data.frame(table(unique_pmhc))
gene_freq <- aggregate(gene_pmhc$Freq, by = list(Gene_Name = gene_pmhc$Gene_Name), sum)
gene_freq <- gene_freq[order(gene_freq$x, decreasing = TRUE), ]
gene_pmhc$Gene_Name <- factor(gene_pmhc$Gene_Name, levels = gene_freq$Gene_Name)

gene_freq1 <- subset(gene_freq, x >= 30)
gene_freq2 <- subset(gene_freq, x < 30)

gene_pmhc1 <- merge(gene_freq1, gene_pmhc, by = "Gene_Name")
gene_pmhc1$Gene_Name <- factor(gene_pmhc1$Gene_Name, levels = gene_freq1$Gene_Name)
gene_pmhc2 <- merge(gene_freq2, gene_pmhc, by = "Gene_Name")
gene_pmhc2$Gene_Name <- factor(gene_pmhc2$Gene_Name, levels = gene_freq2$Gene_Name)

p1 <- ggplot(gene_pmhc1, aes(Gene_Name, Freq, fill = Variant_type)) + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) + 
  labs(x = '', y = 'Number of peptides predicted') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        legend.position = "top") +
  scale_fill_brewer(palette = "Set1")

ggsave(p1, filename = "02.Top_Num_neoantigens_DriversMutation.pdf", width = 6, height = 3.7)

p2 <- ggplot(gene_pmhc2, aes(Gene_Name, Freq, fill = Variant_type)) + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) + 
  labs(x = '', y = 'Number of peptides predicted') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        legend.position = "top") +
  scale_fill_brewer(palette = "Set1")

ggsave(p2, filename = "02.Bottom_Num_neoantigens_DriversMutation.pdf", width = 4, height = 3)

# Pie plot for mutation types
source("D:\\03.R_functions\\JPComUse\\R\\PlotPie.R")
t4 <- gene_pmhc[, 2:3]
t5 <- aggregate(t4$Freq, by = list(group = t4$Variant_type), sum)
colnames(t5) <- c("Group", "value")
p3 <- PlotPie(data = t5, rev.group.order = c("Indel", "SNP"), col.sef = FALSE, palette = "Set1")
ggsave(p3, filename = "02.PiePlot_Indel_SNP.pdf", width = 4.5, height = 4)

# Process fusion-derived p:MHCs
fusion_data <- read_xlsx("Matched_DNA_RNA_fusion_filtered_epitopes.xlsx")
fusion_data <- unique(fusion_data[, c("Gene_Name", "pMHC")])
fusion_data$Gene_Name <- gsub("RARA--PML", "PML--RARA", fusion_data$Gene_Name)
fusion_data$Gene_Name <- gsub("ELL--KMT2A", "KMT2A--ELL", fusion_data$Gene_Name)
fusion_pmhc <- data.frame(table(fusion_data$Gene_Name))
fusion_pmhc <- fusion_pmhc[order(fusion_pmhc$Freq, decreasing = TRUE), ]
fusion_pmhc$Var1 <- factor(fusion_pmhc$Var1, levels = fusion_pmhc$Var1)

p1 <- ggplot(fusion_pmhc, aes(Var1, Freq)) + 
  geom_bar(stat = "identity", alpha = 0.6, width = 0.6) + 
  labs(x = '', y = 'Number of p:MHC predicted') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "top") +
  scale_fill_brewer(palette = "Set3")

ggsave(p1, filename = "02.Num_neoantigens_DriversFusions.pdf", width = 3.8, height = 3.5)

# Correlation neoantigen vs. variants for mutations
mut_neo_data <- read_xlsx("Matched_DNA_RNA_mutation_filtered_epitopes.xlsx")
mut_neo_data <- subset(mut_neo_data, Mut_type == "Driver" & Variant_type != "DNP")
unique_mut_neo <- unique(mut_neo_data[, c("Gene_Name", "HLA_Allele", "MT_Epitope_Seq")])
num_mut_neo <- data.frame(table(unique_mut_neo$Gene))
colnames(num_mut_neo) <- c("Gene", "Mut_neo")

num_mut_variant <- unique(mut_neo_data[, c("Gene_Name", "Index")])
num_mut_variant <- data.frame(table(num_mut_variant$Gene_Name))
colnames(num_mut_variant) <- c("Gene", "Mut_variant")
mat3 <- merge(num_mut_neo, num_mut_variant, by = "Gene")

p1 <- ggplot(mat3, aes(x = log2(Mut_variant), y = log2(Mut_neo), fill = Gene)) + 
  geom_point(alpha = 0.5, shape = 21, stroke = 1, size = 4) + 
  labs(x = "Number of mutation variants", y = "Number of peptides predicted") + 
  theme_classic() +
  geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 25,
                  segment.color = "black", segment.size = 0.15) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey"))

p2 <- p1 + stat_smooth(color = "grey", fill = "skyblue", method = "lm") +
  stat_poly_eq(formula = y ~ x, parse = TRUE)
cor.test(mat3$Mut_neo, mat3$Mut_variant, method = "pearson")
ggsave(p1, filename = "03.corr_num_MutationVariants_neoantigens.pdf", width = 4.5, height = 4)


### Fusions
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Read and process data
mat1 <- read_xlsx("Matched_DNA_RNA_fusion_filtered_epitopes.xlsx")
mat1 <- data.frame(mat1)
mat2 <- unique(mat1[, c("Gene_Name", "HLA_Allele", "MT_Epitope_Seq")])

# Create data frames for fusions and fusion variants
num_Fusions <- data.frame(table(mat2$Gene_Name))
colnames(num_Fusions) <- c("Fusion", "Fusion_neo")
num_fuVa <- unique(mat1[, c("Gene_Name", "hg38_id")])
num_fuVa <- data.frame(table(num_fuVa$Gene_Name))
colnames(num_fuVa) <- c("Fusion", "Fusion_variant")

# Merge data frames
mat3 <- merge(num_Fusions, num_fuVa, by = "Fusion")

# Plot fusion data
p1 <- ggplot(mat3, aes(x = Fusion_variant, y = Fusion_neo, fill = Fusion)) + 
  geom_point(alpha = 0.5, shape = 21, stroke = 1, size = 4) + 
  labs(x = "Number of Fusion variants", y = "Number of peptides predicted") + 
  theme_classic() +
  geom_text_repel(aes(label = Fusion), size = 3, max.overlaps = 25,
                  segment.color = "black", segment.size = 0.15) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey"))

p1
p2 <- p1 + stat_smooth(color = "grey", fill = "skyblue", method = "lm") +
  stat_poly_eq(formula = y ~ x, parse = TRUE)
cor.test(mat3$Fusion_neo, mat3$Fusion_variant, method = "pearson")

ggsave(p1, filename = "03.corr_num_FusionVariants_neoantigens.pdf", width = 4.5, height = 4)

###### Driver gene-specific Length distribution
setwd("D:\\Documents\\Desktop\\V4_AML_mutation_Fusion_Neoantigens\\03_HLA")
t1 <- read_xlsx("DriverMutation_Fusions.xlsx")
t1 <- data.frame(t1)
t2 <- unique(t1[, c("Gene_Name", "Variant_type", "HLA_Allele", "Peptide_Length", "MT_Epitope_Seq")])
t2$HLA_Allele <- gsub("HLA-", "", t2$HLA_Allele)
t2 <- t2[, c("Gene_Name", "Peptide_Length")]

# Frequency tables and ordering
t3 <- data.frame(table(t2))
gene_freq <- data.frame(table(t2$Gene_Name))
gene_freq <- gene_freq[order(gene_freq$Freq, decreasing = TRUE), ]
gene_order1 <- as.character(gene_freq$Var1[c(1:10, 12:18, 20:28, 30:34, 36:39, 42:44, 46:52, 54:55)])
gene_order2 <- as.character(gene_freq$Var1[c(11, 19, 29, 35, 40:41, 45, 53, 56)])
gene_order3 <- c(gene_order1, gene_order2)
t3$Gene_Name <- factor(t3$Gene_Name, levels = gene_order3)

# Plot driver gene-specific length distribution
source("C:\\JP_R.packages\\JPComUse\\R\\PlotBar.R")
p <- PlotBar(data = t3, x.column = "Gene_Name", y.column = "Freq", 
             col.palette = "Set2", fill.column = "Peptide_Length", coord_flip = FALSE)

ggsave(p, filename = "04.DriverslengthDistribution.pdf", width = 7, height = 4)

# Subset and merge data
t4 <- subset(t3, Peptide_Length %in% c(9, 10))
colnames(gene_freq) <- c("Gene_Name", "all_pep")
gene_freq$Type <- c(rep("Mutation", 10), rep("Fusion", 1), rep("Mutation", 7), rep("Fusion", 1),
                    rep("Mutation", 9), rep("Fusion", 1), rep("Mutation", 5), rep("Fusion", 1),
                    rep("Mutation", 4), rep("Fusion", 2), rep("Mutation", 3), rep("Fusion", 1),
                    rep("Mutation", 7), rep("Fusion", 1), rep("Mutation", 2), rep("Fusion", 1))

t5 <- merge(t4, gene_freq, by = "Gene_Name")
t5$ratio <- t5$Freq / t5$all_pep
t5$Type <- factor(t5$Type, levels = c("Mutation", "Fusion"))

# Plot length distribution
p <- ggplot(data = t5, aes(x = Peptide_Length, y = ratio, fill = Peptide_Length)) + 
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.75), size = 0.8, color = "black") +
  geom_boxplot(width = 0.25, notch = FALSE, outlier.size = -1, color = "black", lwd = 0.8, alpha = 0.4) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(), color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_classic() + facet_grid(. ~ Type) +
  stat_compare_means() + 
  theme(legend.position = "none",
        panel.grid.major = element_line(size = 0.1, linetype = "longdash", color = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = "longdash", color = "grey")) +
  ylab("Relative ratio") + xlab("")  

ggsave(p, filename = "04.PepLength_Mutation_Fusions.pdf", width = 3.6, height = 4)

###### Tumor specific neoantigen numbers
# Count number of variants and peptides
t1 <- read_xlsx("DriverMutation_Fusions.xlsx")
t1 <- data.frame(t1)
t2 <- unique(t1[, c(1, 4)])
t3 <- data.frame(table(t2$Variant_type))
t3[1, 2] <- 41  # Adjust this value if necessary

t4 <- unique(t1[, c(3, 4, 2)])
t5 <- data.frame(table(t4$Variant_type))

# Count number of Mutant-specific peptides
t6 <- subset(t1, Fold_Change > 1)
t7 <- unique(t6[, c(1, 4)])
t8 <- data.frame(table(t7$Variant_type))
t9 <- unique(t6[, c(3, 4, 2)])
t10 <- data.frame(table(t9$Variant_type))

# Combine and reshape data
Num_dat <- data.frame(t3, Non_specific_pMHC = t5$Freq)
colnames(Num_dat) <- c("Mut_type", "Non_specific_variant", "Non_specific_pMHC")
Num_dat$Non_specific_ratio <- Num_dat$Non_specific_pMHC / Num_dat$Non_specific_variant
Num_dat$Mut_type <- factor(Num_dat$Mut_type, levels = c("Indel", "Fusion", "SNP"))

# Plot neoantigen numbers
p <- ggplot(Num_dat, aes(Mut_type, weight = Non_specific_ratio, fill = Mut_type)) +
  geom_bar(color = "black", width = 0.7, position = 'dodge') +
  labs(x = "", y = 'p:MHC per variant') +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "top") +
  theme_bw()

ggsave(p, filename = "04.MutType_NumPmhc.pdf", width = 3.7, height = 4)

# Combine and process data
Num_dat <- data.frame(
  Mut_type = t3$Mut_type,
  Non_specific_variant = t5$Freq,
  Non_specific_pMHC = t5$Freq,
  Specific_variant = t8$Freq,
  Specific_pMHC = t10$Freq
)

# Compute ratios
Num_dat$Non_specific_ratio <- Num_dat$Non_specific_pMHC / Num_dat$Non_specific_variant
Num_dat$specific_ratio <- Num_dat$Specific_pMHC / Num_dat$Specific_variant

# Prepare data for plotting
pMHC_dat <- Num_dat[, c("Mut_type", "Non_specific_ratio", "specific_ratio")]
pMHC_dat_long <- melt(pMHC_dat, id.vars = "Mut_type")

# Plot
p <- ggplot(pMHC_dat_long, aes(x = variable, y = value, fill = Mut_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, position = "dodge") +
  labs(x = "", y = "p:MHC per variant") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# Save plot
ggsave(p, filename = "04.Tumor_specific_NumPmhc.pdf", width = 3.7, height = 4)
