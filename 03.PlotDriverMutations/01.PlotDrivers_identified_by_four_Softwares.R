### 2021/7/26
### Jin Peng
### Driver mutation plots using 4 different softwares

# Set working directory
setwd("D:/01.Projects/V3_AML_DriverNeo/02.DriverMutationAnalysis")

# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

# Function to create ggplot for MutSigCV
plot_mutsigcv <- function(data) {
  data$q <- ifelse(data$q == 0, 1e-15, data$q)
  data$color <- brewer.pal(8, "Dark2")[2]
  
  ggplot(data, aes(x=log10(Freq / 2840), y=-log10(q), size=1.3, fill=color)) +
    geom_point(alpha=.5, shape=21, stroke=1) +
    labs(x=bquote(~Log[10]~"(Mutation frequency)"), y=bquote(~-Log[10]~italic("Q-value"))) +
    scale_fill_manual(values=brewer.pal(8, "Dark2")[2]) +
    scale_color_manual(values="black") +
    geom_text_repel(aes(label=gene), size=2, point.padding=unit(0.15, "lines"), max.overlaps=30, segment.color="black", segment.size=0.15) +
    theme_bw(base_size=12, base_family="Times") +
    theme(
      legend.position="none",
      legend.title=element_blank(),
      panel.grid.major=element_line(size=0.1, linetype="longdash", color="grey"),
      panel.grid.minor=element_line(size=0.1, linetype="longdash", color="grey")
    )
}

# Function to create ggplot for 2020plus
plot_2020plus <- function(data) {
  data$q <- ifelse(data$q == 0, 1e-6, data$q)
  
  ggplot(data, aes(x=log10(Freq / 2840), y=-log10(q), size=1.3, fill=type)) +
    geom_point(alpha=.5, shape=21, stroke=1) +
    labs(x=bquote(~Log[10]~"(Mutation frequency)"), y=bquote(~-Log[10]~italic("Q-value"))) +
    scale_fill_manual(values=brewer.pal(8, "Dark2")[3:4]) +
    scale_color_manual(values=c("black", "black")) +
    geom_text_repel(aes(label=gene), size=2, point.padding=unit(0.15, "lines"), max.overlaps=30, segment.color="black", segment.size=0.15) +
    theme_bw(base_size=12, base_family="Times") +
    theme(
      legend.position="none",
      legend.title=element_blank(),
      panel.grid.major=element_line(size=0.1, linetype="longdash", color="grey"),
      panel.grid.minor=element_line(size=0.1, linetype="longdash", color="grey")
    )
}

# Function to create ggplot for OncodriveCLUST
plot_oncodriveclust <- function(data) {
  data$color <- brewer.pal(8, "Dark2")[5]
  
  ggplot(data, aes(x=log10(Freq / 2840), y=-log10(P), size=1.3, fill=color)) +
    geom_point(alpha=.5, shape=21, stroke=1) +
    labs(x=bquote(~Log[10]~"(Mutation frequency)"), y=bquote(~-Log[10]~italic("P-value"))) +
    scale_fill_manual(values=brewer.pal(8, "Dark2")[5]) +
    scale_color_manual(values=c("black", "black")) +
    geom_text_repel(aes(label=gene), size=2, point.padding=unit(0.15, "lines"), max.overlaps=30, segment.color="black", segment.size=0.15) +
    theme_bw(base_size=12, base_family="Times") +
    theme(
      legend.position="none",
      legend.title=element_blank(),
      panel.grid.major=element_line(size=0.1, linetype="longdash", color="grey"),
      panel.grid.minor=element_line(size=0.1, linetype="longdash", color="grey")
    )
}

# Function to create ggplot for OncodriveFML
plot_oncodrivefml <- function(data) {
  data$color <- brewer.pal(8, "Dark2")[6]
  
  ggplot(data, aes(x=log10(Freq / 2840), y=-log10(q), size=1.3, fill=color)) +
    geom_point(alpha=.5, shape=21, stroke=1) +
    labs(x=bquote(~Log[10]~"(Mutation frequency)"), y=bquote(~-Log[10]~italic("P-value"))) +
    scale_fill_manual(values=brewer.pal(8, "Dark2")[6]) +
    scale_color_manual(values=c("black", "black")) +
    geom_text_repel(aes(label=gene), size=2, point.padding=unit(0.15, "lines"), max.overlaps=30, segment.color="black", segment.size=0.15) +
    theme_bw(base_size=12, base_family="Times") +
    theme(
      legend.position="none",
      legend.title=element_blank(),
      panel.grid.major=element_line(size=0.1, linetype="longdash", color="grey"),
      panel.grid.minor=element_line(size=0.1, linetype="longdash", color="grey")
    )
}

# Load data and create plots
mat1 <- read.csv("01.MutSigCV.OUT.csv", header=TRUE)
mat2 <- read.csv("02.2020plus.OUT.csv", header=TRUE)
mat3 <- read.csv("03.OncodriveCLUST.OUT.csv", header=TRUE)
mat4 <- read.csv("04.OncodriveFML.OUT.csv", header=TRUE)

p1 <- plot_mutsigcv(mat1)
p2 <- plot_2020plus(mat2)
p3 <- plot_oncodriveclust(mat3)
p4 <- plot_oncodrivefml(mat4)

# Combine plots
library(patchwork)
combined_plot <- p1 + p2 + p3 + p4
ggsave("01.Driver_mutations_discoveried_4_softwares.pdf", plot=combined_plot, width=9, height=9)

# Plot the frequency identified by 4 softwares
mat <- read.csv("soft_mut.csv", header=TRUE, row.names=1)
order <- read.table("01.mut_freq.txt", header=TRUE)
mat <- as.matrix(mat)
mycol <- pal_npg("nrc")(8)

alter_fun <- list(
  background=function(x, y, w, h) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp=gpar(fill='#ebebe1ff', col=NA)),
  MutSigCV=function(x, y, w, h) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp=gpar(fill=mycol[2], col=NA)),
  '2020Plus'=function(x, y, w, h) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp=gpar(fill=mycol[3], col=NA)),
  OncodriveCLUST=function(x, y, w, h) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp=gpar(fill=mycol[5], col=NA)),
  OncodriveFML=function(x, y, w, h) grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp=gpar(fill=mycol[6], col=NA))
)

col <- c(MutSigCV=mycol[2], '2020Plus'=mycol[3], OncodriveCLUST=mycol[5], OncodriveFML=mycol[7])

p5 <- oncoPrint(mat, col=col, show_pct=FALSE, top_annotation=NULL, right_annotation=NULL, alter_fun=alter_fun, row_order=order$Var1)
pdf("02.Driver_mutations_heatmap.pdf", width=3, height=5)
draw(p5)
dev.off()

# Plot the frequency of mutations
order <- order[order(order$Freq), ]
order$Var1 <- factor(order$Var1, levels=order$Var1)

p6 <- ggplot(order, aes(Var1, Freq)) +
  geom_point(size=1.75, color="#0072B5FF", fill=alpha("lightblue", 0.5), alpha=0.75, shape=21, stroke=0.75) +
  geom_segment(aes(x=Var1, xend=Var1, y=0, yend=Freq), size=0.5, linetype="dashed", color="#838B8B") +
  theme_bw(base_size=12, base_family="Times") +
  theme(legend.position="none", legend.title=element_blank()) +
  coord_flip()
ggsave("03.Mut_freq.pdf", plot=p6, width=2, height = 5.35)
