### Jin Peng
### Clinical data analysis
### 2022.6
# Load required libraries
library(survival)
library(survminer)
library(readxl)

# Load the data
data <- read_xlsx("03_Merged_Cohorts_Neo_Clinical.xlsx")

# Function to create survival plot
#' @param df Data frame with survival data
#' @param cohort Cohort name for the plot title
#' @param neo Column name for neoantigen data
#' @return A ggplot2 object of the survival plot
survival_plot <- function(df, cohort, neo) {
  df$group <- ifelse(df[[neo]] >= 1, 'Yes', 'No')
  df$group <- factor(df$group, levels = c("Yes", "No"))

  # Create survival object and fit the survival model
  surv_obj <- Surv(df$Time, df$Status)
  fit <- survfit(surv_obj ~ group, data = df)
  
  # Perform log-rank test
  surv_diff <- survdiff(surv_obj ~ group, data = df)
  p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  hr <- (surv_diff$obs[1] / surv_diff$exp[1]) / (surv_diff$obs[2] / surv_diff$exp[2])
  up95 <- exp(log(hr) + qnorm(0.975) * sqrt(1 / surv_diff$exp[2] + 1 / surv_diff$exp[1]))
  low95 <- exp(log(hr) - qnorm(0.975) * sqrt(1 / surv_diff$exp[2] + 1 / surv_diff$exp[1]))

  # Format HR and CI strings
  hr_text <- paste("Hazard Ratio = ", round(hr, 2), sep = "")
  ci_text <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

  # Create the plot
  plot <- ggsurvplot(
    fit,
    palette = "npg",
    pval = paste(
      ifelse(p_val < 0.001, "p < 0.001", paste("p = ", round(p_val, 3), sep = "")),
      hr_text, ci_text, sep = "\n"
    ),
    pval.size = 5,
    xlab = "Survival time in days",
    conf.int = FALSE,
    combine = TRUE,
    risk.table = FALSE,
    title = cohort,
    surv.median.line = 'hv',
    legend.title = neo,
    legend.labs = c("Yes (>=1)", "No (=0)"),
    legend = c(0.65, 0.85)
  )
  
  return(plot)
}

# Create a list of plots
plots <- list()

# Define cohorts and neos
cohorts <- c("TCGA", "BeatAML", "RJAML")
neonames <- c("Indel_neo", "SNV_neo", "FusionNeo", "Mutation_Neo", "DriverNeo")

# Function to generate plots for each cohort
generate_plots <- function(cohort_name, data) {
  subset_data <- subset(data, Center == cohort_name)
  plots <- lapply(neonames, function(neo) {
    survival_plot(subset_data, cohort = paste(cohort_name), neo = neo)
  })
  return(plots)
}

# Generate and save plots for each cohort
for (cohort in cohorts) {
  cohort_plots <- generate_plots(cohort, data)
  pdf_filename <- paste0(cohort, "_Survival_plots.pdf")
  pdf(pdf_filename, width = 9, height = 15)
  arrange_ggsurvplots(cohort_plots, print = FALSE, nrow = 5, ncol = 3)
  dev.off()
}

# Combined plots
combined_plots <- lapply(neonames, function(neo) {
  survival_plot(data, cohort = "Combined", neo = neo)
})
pdf("Combined_Survival_plots.pdf", width = 15, height = 3)
arrange_ggsurvplots(combined_plots, print = FALSE, nrow = 1, ncol = 5)
dev.off()

# TCGA DFS plots
tcga_dfs_data <- subset(data, Center == "TCGA")
tcga_dfs_data <- tcga_dfs_data[, -c(3, 4)]
colnames(tcga_dfs_data)[37:38] <- c("Time", "Status")
tcga_dfs_plots <- lapply(neonames, function(neo) {
  survival_plot(tcga_dfs_data, cohort = "TCGA (DFS)", neo = neo)
})
pdf("TCGA_DFS_Survival_plots.pdf", width = 15, height = 3)
arrange_ggsurvplots(tcga_dfs_plots, print = FALSE, nrow = 1, ncol = 5)
dev.off()

# RJAML EFS plots
rjam_efs_data <- subset(data, Center == "RJAML")
rjam_efs_data <- rjam_efs_data[, -c(3, 4)]
colnames(rjam_efs_data)[33:34] <- c("Time", "Status")
rjam_efs_data$Time <- as.numeric(rjam_efs_data$Time)
rjam_efs_data$Status <- as.numeric(rjam_efs_data$Status)
rjam_efs_plots <- lapply(neonames, function(neo) {
  survival_plot(rjam_efs_data, cohort = "RJAML (EFS)", neo = neo)
})
pdf("RJAML_EFS_Survival_plots.pdf", width = 15, height = 3)
arrange_ggsurvplots(rjam_efs_plots, print = FALSE, nrow = 1, ncol = 5)
dev.off()
