#### Single cell RNA-seq analysis
#### Jin Peng
#### 2024

## GSE241989 Analysis - 18 AML Patients with Matched RNA-seq and scRNA-seq

# Load necessary libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)

# Define color palette for clusters
cluster_cols <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F",
                   "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6",
                   "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
                   "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B", 
                   "#66A61E", "#F1788D", "#E6550D", "#E7969C")

# Load and preprocess data
scaml_gse241989 <- readRDS("~/projects/AML_plasticity/jp_2024/scAML_GSE241989_103690Cells.rds")
scaml_gse241989 <- subset(scaml_gse241989, orig.ident %in% c("11H103", "11H097", "06H088", "12H138", "12H106", "05H034", "14H007", "05H066", "09H070", "09H032", "08H087", "05H193", "07H134", "12H010", "03H112", "09H010", "04H096", "02H017"))

# Add mitochondrial percentage
scaml_gse241989[["percent.mt"]] <- PercentageFeatureSet(scaml_gse241989, pattern = "^MT-")

# Plot quality control metrics
VlnPlot(scaml_gse241989, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize, find variable features, and scale data
scaml_gse241989 <- NormalizeData(scaml_gse241989)
scaml_gse241989 <- FindVariableFeatures(scaml_gse241989, selection.method = "vst", nfeatures = 2000)
scaml_gse241989 <- ScaleData(scaml_gse241989)

# Perform PCA and clustering
scaml_gse241989 <- RunPCA(scaml_gse241989, features = VariableFeatures(scaml_gse241989))
scaml_gse241989 <- FindNeighbors(scaml_gse241989, dims = 1:15)
scaml_gse241989 <- FindClusters(scaml_gse241989, resolution = 0.5)
scaml_gse241989 <- RunUMAP(scaml_gse241989, dims = 1:15)

# Generate plots
p0 <- DimPlot(scaml_gse241989, reduction = "umap", label = TRUE)
p1 <- DimPlot(scaml_gse241989, reduction = "umap", group.by = "orig.ident", raster = TRUE, cols = cluster_cols, label = TRUE, pt.size = 2)
p2 <- FeaturePlot(scaml_gse241989, features = c("CD3D", "CD3G", "CD3E"), raster = TRUE, pt.size = 2)

# Save plots
ggsave(p1, filename = "gse241989_dimplot_patients.pdf", width = 6, height = 5)
ggsave(p2, filename = "gse241989_dimplot_CD3.pdf", width = 6, height = 5)

# Extract T cells
t_cells <- subset(scaml_gse241989, idents = "14")
# Save T cells subset
# saveRDS(t_cells, "scAML_GSE241989_TCells.rds")


## GSE116256 Analysis - 4 AML Patients with Matched RNA-seq and scRNA-seq

# Load and preprocess data
scaml_gse116256 <- readRDS("~/projects/AML_plasticity/jp_2024/scAML_GSE116256_30712Cells.rds")
scaml_gse116256 <- subset(scaml_gse116256, orig.ident %in% c("AML556-D0", "AML707B-D0", "AML916-D0", "AML921A-D0"))

# Add mitochondrial percentage
scaml_gse116256[["percent.mt"]] <- PercentageFeatureSet(scaml_gse116256, pattern = "^MT-")

# Plot quality control metrics
VlnPlot(scaml_gse116256, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize, find variable features, and scale data
scaml_gse116256 <- NormalizeData(scaml_gse116256)
scaml_gse116256 <- FindVariableFeatures(scaml_gse116256, selection.method = "vst", nfeatures = 2000)
scaml_gse116256 <- ScaleData(scaml_gse116256)

# Perform PCA and clustering
scaml_gse116256 <- RunPCA(scaml_gse116256, features = VariableFeatures(scaml_gse116256))
scaml_gse116256 <- FindNeighbors(scaml_gse116256, dims = 1:10)
scaml_gse116256 <- FindClusters(scaml_gse116256, resolution = 0.5)
scaml_gse116256 <- RunUMAP(scaml_gse116256, dims = 1:10)

# Generate plots
p0 <- DimPlot(scaml_gse116256, reduction = "umap", label = TRUE, raster = TRUE)
p1 <- DimPlot(scaml_gse116256, reduction = "umap", group.by = "orig.ident", raster = TRUE, cols = cluster_cols, label = TRUE, pt.size = 2)
p2 <- FeaturePlot(scaml_gse116256, features = c("CD3D", "CD3G", "CD3E"), raster = TRUE, pt.size = 2)

# Save plots
ggsave(p1, filename = "gse116256_dimplot_patients.pdf", width = 6.5, height = 5)
ggsave(p2, filename = "gse116256_dimplot_CD3.pdf", width = 6, height = 5)

# Extract T cells
t_cells <- subset(scaml_gse116256, idents = "8")
# Save T cells subset
saveRDS(scaml_gse116256, "~/projects/AML_plasticity/jp_2024/scAML_GSE116256_30712Cells.rds")
saveRDS(t_cells, "scAML_GSE116256_TCells.rds")


## This study 4 patients
## Analysis of scRNA-seq Data for 4 AML Patients
# Load necessary libraries
library(DropletUtils)
library(Seurat)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(BiocParallel)
library(harmony)
library(DoubletFinder)
library(clustree)
library(patchwork)
library(tidyverse)
library(paletteer)

#######################################
# 1. Load data and create Seurat objects #
#######################################

# Define file paths
file_paths <- list(
  XL = "./scRNA_rawdata/MS24051049_result/outs/filtered_feature_bc_matrix.h5",
  ZZH = "./scRNA_rawdata/MS24051050_result/outs/filtered_feature_bc_matrix.h5",
  YLL = "./scRNA_rawdata/MS24051051_result/outs/filtered_feature_bc_matrix.h5",
  LY = "./scRNA_rawdata/MS24051141_result/outs/filtered_feature_bc_matrix.h5"
)

# Load data and create Seurat objects
seurat_list <- lapply(names(file_paths), function(name) {
  data <- Read10X_h5(file_paths[[name]])
  seurat_obj <- CreateSeuratObject(counts = data, project = name)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
})

# Assign names to the list of Seurat objects
names(seurat_list) <- names(file_paths)

####################################################################################################################################
# 2. Identify and remove empty droplets (this step is important for the library with low-quality cells)                             #
####################################################################################################################################

# Process each Seurat object
for (i in seq_along(seurat_list)) {
  seurat_obj <- seurat_list[[i]]
  n_scaml <- seurat_obj$nCount_RNA  # Identify total RNA count for each cell
  l_scaml <- sort(n_scaml)[2]  # Find the second least RNA count as the lower bound
  e_scaml <- emptyDrops(seurat_obj@assays$RNA@counts, lower = l_scaml, BPPARAM = MulticoreParam())  # Identify potential empty droplets
  
  # Filter out empty droplets
  seurat_no_empty <- seurat_obj[, which(e_scaml$FDR <= 0.01)]
  cat(sprintf("Removed %d droplets from sample %s\n", 
              ncol(seurat_obj) - ncol(seurat_no_empty), names(seurat_list)[i]))
  cat(sprintf("Percentage of retained cells: %.2f%%\n", 
              100 * ncol(seurat_no_empty) / ncol(seurat_obj)))
  
  # Identify and remove low-RNA content cells
  seurat_no_empty <- subset(seurat_no_empty, 
                            subset = nFeature_RNA > 500 & 
                            nFeature_RNA < 7000 & 
                            percent.mt < 15)
  
  # Store the cleaned Seurat object
  seurat_list[[i]] <- seurat_no_empty
}

##################################
# 3. Identify and remove doublets #
##################################

# Define a function for processing Seurat objects
process_seurat <- function(seurat_obj) {
  # Standard preprocessing
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  
  # pK Identification
  sweep_res <- paramSweep_v3(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Run DoubletFinder
  doublet_rate <- ncol(seurat_obj) * 8e-6  # ~1000 cells, 0.8% multiplet rate
  n_exp_poi <- round(doublet_rate * ncol(seurat_obj))
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = mpK,
                                  nExp = n_exp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Remove doublets
  seurat_obj$doubletfinder <- seurat_obj@meta.data[, paste0("DF.classifications_0.25_", mpK, "_", n_exp_poi)]
  seurat_singlet <- subset(seurat_obj, doubletfinder == "Singlet")
  
  return(seurat_singlet)
}

# Apply processing to each Seurat object
seurat_list_processed <- lapply(seurat_list, process_seurat)

# Combine Seurat objects
scaml_combined <- merge(seurat_list_processed[[1]],
                        y = seurat_list_processed[-1],
                        add.cell.ids = names(file_paths),
                        project = "scAML")

# Plot quality control metrics
VlnPlot(scaml_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Save the combined Seurat object
saveRDS(scaml_combined, "combined_4patients_qc.rds")

# Load necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load and subset Seurat object
scaml_4pts <- readRDS("combined_4patients_qc.rds")

# Quality control: Filter cells based on feature counts and mitochondrial percentage
scaml_4pts <- subset(scaml_4pts, 
                     subset = nFeature_RNA > 1000 & 
                              nFeature_RNA < 6000 & 
                              percent.mt < 15)

# Recompute mitochondrial percentage feature
scaml_4pts[["percent.mt"]] <- PercentageFeatureSet(scaml_4pts, pattern = "^MT-")

# Plot quality control metrics
VlnPlot(scaml_4pts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Perform standard preprocessing steps
scaml_4pts <- NormalizeData(scaml_4pts)
scaml_4pts <- FindVariableFeatures(scaml_4pts, selection.method = "vst", nfeatures = 2000)
scaml_4pts <- ScaleData(scaml_4pts)
scaml_4pts <- RunPCA(scaml_4pts, features = VariableFeatures(scaml_4pts))
scaml_4pts <- FindNeighbors(scaml_4pts, dims = 1:10)
scaml_4pts <- FindClusters(scaml_4pts, resolution = 0.5)
scaml_4pts <- RunUMAP(scaml_4pts, dims = 1:10)

# Define color palette for clusters
clusterCols <- c("#FF7F00", "#1B9E77", "#FDAE6B", 
                 "#66A61E", "#F1788D", "#E6550D", "#E7969C",
                 "#843C39", "#8C6D31", "#3182BD", "#54990F",
                 "#BD9E39", "#E7BA52", "#31A354", "#E41A1C",
                 "#6BAED6", "#9ECAE1", "#AD494A", "#A1D99B",
                 "#C7E9C0", "#99600F", "#C3BC3F", "#D6616B")

# Create UMAP plots
p0 <- DimPlot(scaml_4pts, reduction = "umap", label = TRUE, raster = TRUE)
p1 <- DimPlot(scaml_4pts, reduction = "umap", group.by = "orig.ident", raster = TRUE, cols = clusterCols, label = TRUE, pt.size = 2)
p2 <- FeaturePlot(scaml_4pts, features = c("CD3D", "CD3G", "CD3E"), raster = TRUE, pt.size = 2)

# Subset T cells based on cluster IDs
t_cells <- subset(scaml_4pts, seurat_clusters %in% c(9, 13, 16, 20, 23))

# Save results
saveRDS(scaml_4pts, "combined_4patients_qc_42329cells.rds")
saveRDS(t_cells, "combined_4patients_TCells.rds")

# Save plots
ggsave(p1, filename = "combined_4patients_dimplot_patients.pdf", width = 6, height = 5.2)
ggsave(p2, filename = "combined_4patients_dimplot_CD3.pdf", width = 6, height = 5)


## Projectils Analysis
####################################
###### Neo_high vs. Neo_low ######
####################################

# Load required libraries
library(ProjecTILs)
library(Seurat)
library(ggplot2)
library(reshape2)

# Set working directory
setwd("D:\\01.Projects\\V4_AML_mutation_Fusion_Neoantigens\\07_scRNA_Bulk_Neo_CD8T")

# Define number of cores for parallel processing
ncores <- 4

# Load reference map
ref <- readRDS("ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")

# Plot reference map
p0 <- DimPlot(ref, label = FALSE, cols = refCols) + NoLegend()
ggsave(p0, filename = "00_ref.map.pdf", width = 4, height = 4)

# Load datasets
dat1 <- readRDS("scAML_GSE116256_TCells.rds")
dat2 <- readRDS("scAML_GSE241989_TCells.rds")
dat3 <- readRDS("combined_4patients_TCells.rds")

# Merge datasets
scaml.list <- list(GSE116256 = dat1, GSE241989 = dat2, `This study` = dat3)
scaml <- merge(scaml.list[[1]], y = scaml.list[-1], project = "scAML")

# Assign neoantigen groups
scaml$neo_group <- ifelse(scaml$orig.ident %in% c("07H134", "08H087", "09H070", "14H007", "AML556-D0",
                                                  "AML916-D0", "AML921A-D0", "09H010", "YLL", "ZZH"), 
                          "Neo_high", "Neo_low")

# Split data by neoantigen group
data.split <- SplitObject(scaml, split.by = "neo_group")

# Make projections using reference map
query.projected <- make.projection(data.split, ref, filter.cells = FALSE, ncores = ncores)

# Predict cell states
query.projected <- lapply(query.projected, function(x) {
  cellstate.predict(ref = ref, query = x, reduction = "umap", ndim = 2)
})

# Merge projections
query.projected.merged <- suppressMessages(Reduce(ProjecTILs:::merge.Seurat.embeddings, query.projected))

# Split merged projections by neoantigen group
query.sub.byExp <- SplitObject(query.projected.merged, split.by = "neo_group")
n_neohigh <- length(table(query.sub.byExp$Neo_high$orig.ident))
n_neolow <- length(table(query.sub.byExp$Neo_low$orig.ident))

# Plot projections
a <- plot.projection(ref, query = query.sub.byExp$Neo_high, linesize = 0.5, pointsize = 0.25) +
  ggtitle(paste0("Neo-High; n=", n_neohigh))
b <- plot.projection(ref, query = query.sub.byExp$Neo_low, linesize = 0.5, pointsize = 0.25) +
  ggtitle(paste0("Neo-Low; n=", n_neolow))

# Plot state composition
c <- plot.statepred.composition(ref, query.sub.byExp$Neo_high, metric = "percent") + ylim(0, 65) + 
  ggtitle(paste0("Neo-High; n=", n_neohigh)) + theme_classic() + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
d <- plot.statepred.composition(ref, query.sub.byExp$Neo_low, metric = "percent") + ylim(0, 65) + 
  ggtitle(paste0("Neo-Low; n=", n_neolow)) + theme_classic() + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))

# Combine and save plots
plot1 <- a + b 
plot2 <- c + d
# ggsave(plot1, filename = "01_AML_T_Cells_projections.pdf", width = 12, height = 4.5)
# ggsave(plot2, filename = "01_AML_T_Cells_barplots.pdf", width = 12, height = 4.5)

# Compute and plot fold changes
which.types <- table(query.projected.merged$functional.cluster) > 30
stateColors_func <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc",
                      "#FF0000", "#87f6a5", "#e812dd")
states_all <- levels(ref$functional.cluster)
names(stateColors_func) <- states_all
cols_use <- stateColors_func[names(which.types)][which.types]

query.projected.merged$functional.cluster <- factor(query.projected.merged$functional.cluster, levels = states_all)
query.list <- SplitObject(query.projected.merged, split.by = "neo_group")

norm.neohigh <- table(query.list[["Neo_high"]]$functional.cluster) / sum(table(query.list[["Neo_high"]]$functional.cluster))
norm.neolow <- table(query.list[["Neo_low"]]$functional.cluster) / sum(table(query.list[["Neo_low"]]$functional.cluster))

foldchange <- norm.neohigh[which.types] / norm.neolow[which.types]
foldchange <- sort(foldchange, decreasing = TRUE)

tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cell_state", "Fold_change")

p <- ggplot(tb.m, aes(x = Cell_state, y = Fold_change, fill = Cell_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols_use) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = "log2") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "right") +
  ggtitle("AML vs. Normal")

# Save fold change plot
ggsave(p, filename = "02_NeoHigh_NeoLow_T_subsets_FC.pdf", width = 6, height = 4)

