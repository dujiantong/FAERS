#' FAERS Data Visualization
#' @description Publication-ready figures

library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

source("R/utils.R")

OUTPUT_DIR <- "output/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load results
signals <- fread("output/tables/significant_signals.csv")

#' Forest plot for top signals
plot_forest <- function(data, drug_filter = NULL, top_n = 20) {
  if (!is.null(drug_filter)) {
    plot_data <- data[drug == drug_filter]
  } else {
    plot_data <- head(data[order(-ror)], top_n)
  }
  
  plot_data <- plot_data %>%
    .[, event := factor(event, levels = rev(event))] %>%
    .[ror_ci_lower > 0]  # Ensure positive for log scale
  
  p <- ggplot(plot_data, aes(x = event, y = ror)) +
    geom_point(aes(size = case_count), color = "#2E86AB") +
    geom_errorbar(aes(ymin = ror_ci_lower, ymax = ror_ci_upper), 
                  width = 0.2, color = "#A23B72") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_y_log10() +
    labs(
      title = sprintf("Disproportionality Analysis: %s", 
                     ifelse(is.null(drug_filter), "Top Signals", drug_filter)),
      x = "MedDRA Preferred Term",
      y = "Reporting Odds Ratio (ROR, log scale)",
      size = "Case Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

# Generate plots
p1 <- plot_forest(signals, top_n = 15)
ggsave(file.path(OUTPUT_DIR, "forest_plot_top_signals.pdf"), p1, 
       width = 10, height = 8)

# Drug-specific plots
for (drug in unique(signals$drug)) {
  p <- plot_forest(signals, drug_filter = drug, top_n = 10)
  ggsave(file.path(OUTPUT_DIR, sprintf("forest_%s.pdf", drug)), p,
         width = 10, height = 6)
}

#' Volcano plot
plot_volcano <- function(data) {
  plot_data <- data %>%
    .[, `:=`(
      log_ror = log2(ror),
      log_p = -log10(0.05 / nrow(data)),  # Bonferroni corrected
      highlight = ror > 5 & case_count > 10
    )]
  
  p <- ggplot(plot_data, aes(x = log_ror, y = case_count)) +
    geom_point(aes(color = highlight, size = ror), alpha = 0.6) +
    geom_text_repel(data = plot_data[highlight == TRUE],
                    aes(label = event), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("grey60", "#F18F01")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      title = "Signal Detection Volcano Plot",
      x = "Log2(ROR)",
      y = "Number of Cases",
      color = "Significant Signal"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

p2 <- plot_volcano(signals)
ggsave(file.path(OUTPUT_DIR, "volcano_plot.pdf"), p2, width = 10, height = 8)

message("Visualizations saved to ", OUTPUT_DIR)
