# ==============================================================================
# 4_Publication_Figures_CORRECTED_3B.R
# FA versus source-annotated FVPTC
#
# Purpose
# -------
# Generate publication figures from the methodologically corrected ML outputs
# produced by 3B_METHOD_CORRECTIONS_PROFESSOR.R.
#
# This script:
# - does NOT clear the R workspace;
# - does NOT retrain any machine-learning model;
# - keeps the initial seven-model screen from the original 3.R outputs;
# - uses nested repeated-CV, cross-fitted PI, normalized Borda, rank-stability,
#   panel-size, and all-8-gene correlation outputs from stage2_methodologically_corrected;
# - writes figures to publication-figures/clean-final-corrected so that the old
#   clean-final folder is preserved.
# ============================================================================== 

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales",
  "grid",
  "pROC"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)
library(pROC)

# ------------------------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------------------------

project_dir <- getwd()
stage1_dir <- file.path(
  project_dir,
  "ml-output",
  "stage1_screening_reproducible"
)

stage2_old_dir <- file.path(
  project_dir,
  "ml-output",
  "stage2_top3_final_reproducible"
)

stage2_corrected_dir <- file.path(
  project_dir,
  "ml-output",
  "stage2_methodologically_corrected"
)

figure_dir <- file.path(
  project_dir,
  "publication-figures",
  "clean-final-corrected"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(figure_dir)) {
  stop("The figure output directory could not be created:\n", figure_dir)
}

# ------------------------------------------------------------------------------
# 3. General settings and helper functions
# ------------------------------------------------------------------------------

font_family <- "Arial"
export_dpi <- 600

# Use these exact eight genes throughout the corrected analysis.
gene_list <- c(
  "FN1",
  "KRT19",
  "TIMP1",
  "CHI3L1",
  "CLU",
  "SPP1",
  "APOE",
  "EFEMP1"
)

group_colors <- c(
  "FA" = "#3B78A8",
  "FVPTC" = "#D9534F"
)

model_colors <- c(
  "XGBoost" = "#0072B2",
  "Naive Bayes" = "#D55E00",
  "KNN" = "#009E73"
)

# Shared default model colors used by multiple figure sections.
default_model_colors <- model_colors

panel_colors <- c(
  "Primary" = "#1B9E77",
  "Nested core" = "#66A61E",
  "Secondary" = "#7570B3",
  "Weak/model-limited" = "#E6AB02",
  "Not supported" = "#999999"
)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

wrap_gene_names <- function(x, width = 35) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return("None")
  }
  paste(strwrap(paste(x, collapse = ", "), width = width), collapse = "\n")
}

read_required_csv <- function(path, required_columns = NULL) {
  if (!file.exists(path)) {
    stop("Required file was not found:\n", path)
  }
  output <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(required_columns)) {
    missing_columns <- setdiff(required_columns, colnames(output))
    if (length(missing_columns) > 0) {
      stop(
        "Missing columns in:\n", path,
        "\nMissing: ", paste(missing_columns, collapse = ", ")
      )
    }
  }
  output
}

save_publication_plot <- function(plot, filename_base, width, height, dpi = export_dpi) {
  output_base <- file.path(figure_dir, filename_base)
  
  ggsave(
    filename = paste0(output_base, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = paste0(output_base, ".tiff"),
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  
  ggsave(
    filename = paste0(output_base, ".pdf"),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

format_p_value <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
  )
}

parse_panel_genes <- function(panel_table, panel_name) {
  row <- panel_table %>% filter(.data$Panel == panel_name)
  if (nrow(row) == 0) {
    return(character(0))
  }
  genes_string <- as.character(row$Genes[1])
  if (is.na(genes_string) || !nzchar(genes_string)) {
    return(character(0))
  }
  trimws(unlist(strsplit(genes_string, ",")))
}

base_theme <- theme_classic(base_size = 11, base_family = font_family) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 12.5, hjust = 0),
    plot.subtitle = element_text(size = 9.3, color = "gray30", margin = margin(b = 7)),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.2),
    plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
  )

# ------------------------------------------------------------------------------
# 4. Required inputs
# ------------------------------------------------------------------------------

required_input_files <- c(
  file.path(stage1_dir, "Table_S0_Initial_7_Models_Run_Level_ROC.csv"),
  file.path(stage1_dir, "Table_S0_Initial_7_Models_Numeric_Summary.csv"),
  file.path(stage2_old_dir, "Table_Final_Selected_Models_From_Stage1.csv"),
  file.path(stage2_old_dir, "Report_Sample_Label_Order_Used_By_ML.csv"),
  file.path(stage2_corrected_dir, "Table_NestedCV_Final_Model_Performance_Aggregated.csv"),
  file.path(stage2_corrected_dir, "Table_NestedCV_Nadeau_Bengio_Corrected_Tests.csv"),
  file.path(stage2_corrected_dir, "Table_NestedCV_Aggregated_Predictions.csv"),
  file.path(stage2_corrected_dir, "Table_Corrected_Consensus_Normalized_Borda.csv"),
  file.path(stage2_corrected_dir, "Table_Selected_Panels_Corrected.csv"),
  file.path(stage2_corrected_dir, "Table_CrossFitted_PI_Summary_All_Models.csv"),
  file.path(stage2_corrected_dir, "Table_CrossFitted_PI_Rank_Stability.csv"),
  file.path(stage2_corrected_dir, "Table_All8_Genes_Spearman_Correlation_Matrix.csv"),
  file.path(stage2_corrected_dir, "Table_Panel_Size_Sensitivity_Performance.csv"),
  file.path(stage2_corrected_dir, "Table_Corrected_Univariate_Expression_Wilcoxon.csv"),
  file.path(project_dir, "final_data_clean.rds")
)

missing_input_files <- required_input_files[!file.exists(required_input_files)]

if (length(missing_input_files) > 0) {
  stop(
    "These required corrected-figure source files are missing:\n",
    paste(missing_input_files, collapse = "\n")
  )
}

# ------------------------------------------------------------------------------
# 5. Load shared data
# ------------------------------------------------------------------------------

selected_models_table <- read_required_csv(
  file.path(stage2_old_dir, "Table_Final_Selected_Models_From_Stage1.csv"),
  required_columns = c("Final_model_rank", "Model")
) %>%
  arrange(.data$Final_model_rank)

final_model_names <- as.character(selected_models_table$Model)

if (!identical(final_model_names, c("XGBoost", "Naive Bayes", "KNN"))) {
  stop(
    "Expected final models to be XGBoost, Naive Bayes, KNN. Found:\n",
    paste(final_model_names, collapse = ", ")
  )
}

model_colors <- model_colors[final_model_names]

sample_labels <- read_required_csv(
  file.path(stage2_old_dir, "Report_Sample_Label_Order_Used_By_ML.csv"),
  required_columns = c("Sample_order", "Sample_name", "Label")
) %>%
  arrange(.data$Sample_order) %>%
  mutate(
    Label = factor(.data$Label, levels = c("FA", "FVPTC")),
    Batch = factor(
      c(rep("GSE29315", 30), rep("GSE54958", 13)),
      levels = c("GSE29315", "GSE54958")
    )
  )

if (nrow(sample_labels) != 43) {
  stop("This script expects exactly 43 samples for the 30/13 dataset assignment.")
}

performance_table <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_NestedCV_Final_Model_Performance_Aggregated.csv"),
  required_columns = c("Model", "AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "MCC", "TP", "TN", "FP", "FN")
) %>%
  mutate(Model = factor(.data$Model, levels = final_model_names)) %>%
  arrange(.data$Model)

pairwise_statistics <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_NestedCV_Nadeau_Bengio_Corrected_Tests.csv"),
  required_columns = c("Comparison", "Mean_AUC_difference", "Corrected_p_value", "N_paired_resamples")
)

corrected_predictions <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_NestedCV_Aggregated_Predictions.csv"),
  required_columns = c("Model", "rowIndex", "Sample", "obs", "Mean_Probability_FVPTC")
) %>%
  mutate(
    Model = factor(.data$Model, levels = final_model_names),
    obs = factor(.data$obs, levels = c("FA", "FVPTC")),
    pred = factor(
      ifelse(.data$Mean_Probability_FVPTC >= 0.5, "FVPTC", "FA"),
      levels = c("FA", "FVPTC")
    )
  )

consensus_table <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_Corrected_Consensus_Normalized_Borda.csv"),
  required_columns = c("Consensus_rank", "Gene", "Normalized_consensus_score", "Score_drop", "Positive_models", "Wilcoxon_BH_FDR")
) %>%
  arrange(.data$Consensus_rank)

selected_panels_table <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_Selected_Panels_Corrected.csv"),
  required_columns = c("Panel", "Genes")
)

primary_panel <- parse_panel_genes(
  selected_panels_table,
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- parse_panel_genes(
  selected_panels_table,
  "Corrected_secondary_candidates"
)

zero_consensus_genes <- parse_panel_genes(
  selected_panels_table,
  "Zero_normalized_consensus_genes"
)

nested_core_panel <- primary_panel[primary_panel %in% c("CHI3L1", "FN1", "EFEMP1", "CLU")]

if (length(primary_panel) == 0) {
  stop("Could not read corrected primary panel from Table_Selected_Panels_Corrected.csv")
}

pi_summary <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_CrossFitted_PI_Summary_All_Models.csv"),
  required_columns = c("Model", "Gene", "Mean", "SD", "Positive_fraction", "Positive_permutation_fraction")
) %>%
  mutate(
    Model = factor(.data$Model, levels = final_model_names),
    Gene = factor(.data$Gene, levels = rev(consensus_table$Gene))
  )

rank_stability <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_CrossFitted_PI_Rank_Stability.csv"),
  required_columns = c("Model", "Gene", "Positive_resample_frequency", "Top4_resample_frequency", "Median_positive_rank", "Rank_IQR")
) %>%
  mutate(
    Model = factor(.data$Model, levels = final_model_names),
    Gene = factor(.data$Gene, levels = rev(consensus_table$Gene))
  )

correlation_matrix <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_All8_Genes_Spearman_Correlation_Matrix.csv")
)

panel_size_table <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_Panel_Size_Sensitivity_Performance.csv"),
  required_columns = c("Number_genes", "Genes", "Model", "AUC", "Within_0_01_of_model_maximum")
) %>%
  mutate(Model = factor(.data$Model, levels = final_model_names))

wilcoxon_table <- read_required_csv(
  file.path(stage2_corrected_dir, "Table_Corrected_Univariate_Expression_Wilcoxon.csv"),
  required_columns = c("Gene", "Median_FA", "Median_FVPTC", "Expression_direction", "Wilcoxon_BH_FDR")
)

final_data <- readRDS(file.path(project_dir, "final_data_clean.rds"))

missing_genes <- setdiff(gene_list, rownames(final_data))
if (length(missing_genes) > 0) {
  stop("These genes are missing from final_data_clean.rds:\n", paste(missing_genes, collapse = ", "))
}

if (!identical(as.character(sample_labels$Sample_name), colnames(final_data))) {
  stop("Sample order mismatch between final_data_clean.rds and Report_Sample_Label_Order_Used_By_ML.csv")
}

expression_data <- as.data.frame(
  t(final_data[gene_list, , drop = FALSE]),
  stringsAsFactors = FALSE
)

for (gene_name in gene_list) {
  expression_data[[gene_name]] <- as.numeric(expression_data[[gene_name]])
}

expression_data$Sample <- sample_labels$Sample_name
expression_data$rowIndex <- sample_labels$Sample_order
expression_data$Group <- sample_labels$Label
expression_data$Batch <- sample_labels$Batch











# ==============================================================================
# FIGURE 1
# CORRECTED WORKFLOW OF BIOMARKER DISCOVERY AND INTERNAL ML PRIORITIZATION
# ==============================================================================

cat("\n--- Generating FIGURE 1: Corrected workflow ---\n")

# ------------------------------------------------------------------------------
# 4. Load workflow information from corrected 3B outputs
# ------------------------------------------------------------------------------

selected_models_file <- file.path(
  stage2_old_dir,
  "Table_Final_Selected_Models_From_Stage1.csv"
)

selected_panels_corrected_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

workflow_required_files <- c(
  selected_models_file,
  selected_panels_corrected_file
)

workflow_missing_files <- workflow_required_files[
  !file.exists(workflow_required_files)
]

if (length(workflow_missing_files) > 0) {
  stop(
    "These workflow source files are missing:\n",
    paste(workflow_missing_files, collapse = "\n")
  )
}

selected_models_table <- read.csv(
  selected_models_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(
  c("Final_model_rank", "Model") %in%
  colnames(selected_models_table)
)) {
  stop(
    "The selected-model table must contain ",
    "'Final_model_rank' and 'Model'."
  )
}

selected_models_table <- selected_models_table %>%
  arrange(Final_model_rank)

final_model_names <- as.character(
  selected_models_table$Model
)

selected_panels_corrected <- read.csv(
  selected_panels_corrected_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("Panel", "Genes") %in% colnames(selected_panels_corrected))) {
  stop(
    "The corrected selected-panel table must contain ",
    "'Panel' and 'Genes'."
  )
}

extract_corrected_genes <- function(panel_name) {
  
  selected_row <- selected_panels_corrected %>%
    filter(Panel == panel_name)
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

primary_panel <- extract_corrected_genes(
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- extract_corrected_genes(
  "Corrected_secondary_candidates"
)

expression_only_genes <- extract_corrected_genes(
  "Expression_confirmed_not_ML_supported"
)

zero_consensus_genes <- extract_corrected_genes(
  "Zero_normalized_consensus_genes"
)

final_models_text <- wrap_gene_names(
  final_model_names,
  width = 33
)

primary_panel_text <- wrap_gene_names(
  primary_panel,
  width = 33
)

secondary_text <- wrap_gene_names(
  secondary_candidates,
  width = 33
)

expression_only_text <- wrap_gene_names(
  expression_only_genes,
  width = 33
)

zero_consensus_text <- wrap_gene_names(
  zero_consensus_genes,
  width = 33
)

# ------------------------------------------------------------------------------
# 5. Workflow drawing helpers
# ------------------------------------------------------------------------------

workflow_blue <- list(
  border = "#2874A6",
  face = "#EBF5FB",
  title = "#154360"
)

workflow_orange <- list(
  border = "#D68910",
  face = "#FEF5E7",
  title = "#7E5109"
)

workflow_green <- list(
  border = "#229954",
  face = "#E9F7EF",
  title = "#145A32"
)

workflow_gray <- "#34495E"
workflow_shadow <- "#BDC3C7"
workflow_light_gray <- "#7F8C8D"

add_phase_background <- function(
    x,
    y,
    width,
    height,
    edge_color,
    face_color,
    title,
    title_color
) {
  
  grid.roundrect(
    x = unit(x + width / 2, "native"),
    y = unit(y + height / 2, "native"),
    width = unit(width, "native"),
    height = unit(height, "native"),
    r = unit(0.015, "snpc"),
    gp = gpar(
      col = edge_color,
      fill = grDevices::adjustcolor(
        face_color,
        alpha.f = 0.65
      ),
      lwd = 1.8,
      lty = "dashed"
    )
  )
  
  grid.text(
    label = title,
    x = unit(x + width / 2, "native"),
    y = unit(y + height - 0.15, "native"),
    just = c("center", "top"),
    gp = gpar(
      fontfamily = font_family,
      fontface = "bold",
      fontsize = 15,
      col = title_color
    )
  )
}

add_workflow_box <- function(
    x,
    y,
    width,
    height,
    title,
    body,
    edge_color,
    face_color = "white",
    title_size = 15,  # اندازه تیتر باز هم بزرگ‌تر شد
    body_size = 13,   # اندازه متن باز هم بزرگ‌تر شد
    shadow = TRUE
) {
  
  if (isTRUE(shadow)) {
    grid.roundrect(
      x = unit(x + width / 2 + 0.04, "native"),
      y = unit(y + height / 2 - 0.04, "native"),
      width = unit(width, "native"),
      height = unit(height, "native"),
      r = unit(0.018, "snpc"),
      gp = gpar(
        col = NA,
        fill = grDevices::adjustcolor(
          workflow_shadow,
          alpha.f = 0.42
        )
      )
    )
  }
  
  grid.roundrect(
    x = unit(x + width / 2, "native"),
    y = unit(y + height / 2, "native"),
    width = unit(width, "native"),
    height = unit(height, "native"),
    r = unit(0.018, "snpc"),
    gp = gpar(
      col = edge_color,
      fill = face_color,
      lwd = 2
    )
  )
  
  if (!is.null(title) && nzchar(title)) {
    grid.text(
      label = title,
      x = unit(x + width / 2, "native"),
      y = unit(y + height - 0.20, "native"),
      just = c("center", "top"),
      gp = gpar(
        fontfamily = font_family,
        fontface = "bold",
        fontsize = title_size,
        col = edge_color
      )
    )
  }
  
  body_y <- if (
    !is.null(title) &&
    nzchar(title)
  ) {
    y + height / 2 - 0.08
  } else {
    y + height / 2
  }
  
  grid.text(
    label = body,
    x = unit(x + width / 2, "native"),
    y = unit(body_y, "native"),
    just = "center",
    gp = gpar(
      fontfamily = font_family,
      fontsize = body_size,
      col = workflow_gray,
      lineheight = 1.1 # فاصله خطوط برای جا شدن متن‌های بزرگتر فشرده‌تر شد
    )
  )
}

add_workflow_arrow <- function(
    x1,
    y1,
    x2,
    y2,
    color = workflow_gray,
    line_width = 2.2
) {
  
  grid.lines(
    x = unit(c(x1, x2), "native"),
    y = unit(c(y1, y2), "native"),
    arrow = arrow(
      type = "closed",
      ends = "last",
      length = unit(0.11, "inches")
    ),
    gp = gpar(
      col = color,
      lwd = line_width,
      lineend = "round"
    )
  )
}

add_workflow_path_arrow <- function(
    x,
    y,
    color = workflow_gray,
    line_width = 2.2
) {
  
  grid.lines(
    x = unit(x, "native"),
    y = unit(y, "native"),
    arrow = arrow(
      type = "closed",
      ends = "last",
      length = unit(0.11, "inches")
    ),
    gp = gpar(
      col = color,
      lwd = line_width,
      lineend = "round",
      linejoin = "round"
    )
  )
}

draw_workflow <- function() {
  
  grid.newpage()
  
  pushViewport(
    viewport(
      xscale = c(0, 12),
      yscale = c(0, 9)
    )
  )
  
  grid.text(
    label = paste0(
      "Workflow of Biomarker Discovery and ",
      "Internal Machine Learning Prioritization"
    ),
    x = unit(6, "native"),
    y = unit(8.85, "native"),
    gp = gpar(
      fontfamily = font_family,
      fontface = "bold",
      fontsize = 18,
      col = "#111111"
    )
  )
  
  add_phase_background(
    x = 0.20,
    y = 6.15,
    width = 11.60,
    height = 2.45,
    edge_color = workflow_blue$border,
    face_color = workflow_blue$face,
    title = "Phase 1: Bioinformatics Discovery",
    title_color = workflow_blue$title
  )
  
  add_phase_background(
    x = 0.20,
    y = 3.25,
    width = 11.60,
    height = 2.45,
    edge_color = workflow_orange$border,
    face_color = workflow_orange$face,
    title = "Phase 2: Internal Machine Learning Analysis",
    title_color = workflow_orange$title
  )
  
  add_phase_background(
    x = 0.20,
    y = 0.35,
    width = 11.60,
    height = 2.45,
    edge_color = workflow_green$border,
    face_color = workflow_green$face,
    title = "Phase 3: Candidate Biomarker Prioritization and Reporting",
    title_color = workflow_green$title
  )
  
  # Phase 1
  add_workflow_box(
    x = 0.60,
    y = 6.55,
    width = 2.80,
    height = 1.40,
    title = "Public GEO Datasets",
    body = paste0(
      "GSE29315: FA (n=17), FVPTC (n=13)\n",
      "GSE54958: FA (n=7), FVPTC (n=6)\n",
      "Total: n=43"
    ),
    edge_color = workflow_blue$border
  )
  
  add_workflow_box(
    x = 4.60,
    y = 6.55,
    width = 2.80,
    height = 1.40,
    title = "Preprocessing",
    body = paste0(
      "RMA Normalization\n",
      "Probe-to-Gene Annotation\n",
      "Shared-Gene Alignment"
    ),
    edge_color = workflow_blue$border
  )
  
  add_workflow_box(
    x = 8.60,
    y = 6.55,
    width = 2.80,
    height = 1.40,
    title = "Integrated Discovery",
    body = paste0(
      "ComBat Batch Correction\n",
      "limma DEG Analysis\n",
      "Candidate Hub Genes (n=8)"
    ),
    edge_color = workflow_blue$border,
    face_color = "#D4E6F1"
  )
  
  add_workflow_arrow(
    x1 = 3.40,
    y1 = 7.25,
    x2 = 4.60,
    y2 = 7.25
  )
  
  add_workflow_arrow(
    x1 = 7.40,
    y1 = 7.25,
    x2 = 8.60,
    y2 = 7.25
  )
  
  # Phase 2
  add_workflow_box(
    x = 0.60,
    y = 3.65,
    width = 2.80,
    height = 1.40,
    title = "Initial Screening",
    body = "7 Classifiers\n10 × 5-Fold CV",
    edge_color = workflow_orange$border
  )
  
  add_workflow_box(
    x = 4.60,
    y = 3.65,
    width = 2.80,
    height = 1.40,
    title = "Final Top-3 Models",
    body = paste0(
      final_models_text,
      "\nOuter 5-Fold CV × 20\n",
      "Inner 4-Fold Tuning"
    ),
    edge_color = workflow_orange$border,
    body_size = 12.5 # سایز اختصاصی متناسب با ارتقا
  )
  
  add_workflow_box(
    x = 8.60,
    y = 3.65,
    width = 2.80,
    height = 1.40,
    title = "ML Evidence",
    body = paste0(
      "Nested ROC-AUC and Metrics\n",
      "Nadeau–Bengio Tests\n",
      "Held-Out Permutation ΔAUC"
    ),
    edge_color = workflow_orange$border,
    face_color = "#FAD7A0",
    body_size = 12.5 # سایز اختصاصی متناسب با ارتقا
  )
  
  add_workflow_path_arrow(
    x = c(10.00, 10.00, 2.00, 2.00),
    y = c(6.55, 5.60, 5.60, 5.05)
  )
  
  add_workflow_arrow(
    x1 = 3.40,
    y1 = 4.35,
    x2 = 4.60,
    y2 = 4.35
  )
  
  add_workflow_arrow(
    x1 = 7.40,
    y1 = 4.35,
    x2 = 8.60,
    y2 = 4.35
  )
  
  # Phase 3
  add_workflow_box(
    x = 0.60,
    y = 0.85,
    width = 2.80,
    height = 1.40,
    title = "Consensus Prioritization",
    body = paste0(
      "Model-Normalized Borda\n",
      "Rank Stability\n",
      "Panel-Size Sensitivity"
    ),
    edge_color = workflow_green$border,
    body_size = 12.5 # سایز اختصاصی متناسب با ارتقا
  )
  
  add_workflow_box(
    x = 4.60,
    y = 0.85,
    width = 2.80,
    height = 1.40,
    title = "Expression Evidence",
    body = paste0(
      "Wilcoxon BH-FDR\n",
      "Direction-Aligned ROC\n",
      "All-8-Gene Correlation"
    ),
    edge_color = workflow_green$border,
    body_size = 12.5 # سایز اختصاصی متناسب با ارتقا
  )
  
  add_workflow_box(
    x = 8.35,
    y = 0.90,
    width = 3.15,
    height = 1.35,
    title = "Primary Panel",
    body = primary_panel_text,
    edge_color = workflow_green$border,
    face_color = "#A9DFBF",
    body_size = 13.5 # سایز اختصاصی متناسب با ارتقا
  )
  
  add_workflow_box(
    x = 8.35,
    y = 0.10,
    width = 3.15,
    height = 0.65,
    title = NULL,
    body = paste0(
      "Secondary: ",
      gsub("\n", " ", secondary_text),
      "\nExpression-only: ",
      gsub("\n", " ", expression_only_text)
    ),
    edge_color = workflow_light_gray,
    face_color = "#F8F9F9",
    body_size = 11.5, # سایز اختصاصی متناسب با ارتقا
    shadow = FALSE
  )
  
  add_workflow_path_arrow(
    x = c(10.00, 10.00, 2.00, 2.00),
    y = c(3.65, 2.75, 2.75, 2.25)
  )
  
  add_workflow_arrow(
    x1 = 3.40,
    y1 = 1.55,
    x2 = 4.60,
    y2 = 1.55
  )
  
  add_workflow_arrow(
    x1 = 7.40,
    y1 = 1.55,
    x2 = 8.35,
    y2 = 1.55
  )
  
  grid.text(
    label = paste0(
      "Note: The five-gene panel represents internal ",
      "prioritization and requires independent external validation."
    ),
    x = unit(0.35, "native"),
    y = unit(0.03, "native"),
    just = c("left", "bottom"),
    gp = gpar(
      fontfamily = font_family,
      fontface = "italic",
      fontsize = 12, # سایز پانویس متناسب با ارتقا
      col = workflow_light_gray
    )
  )
  
  popViewport()
}

# ------------------------------------------------------------------------------
# 6. Display and save corrected Figure 1
# ------------------------------------------------------------------------------

draw_workflow()

workflow_png <- file.path(
  figure_dir,
  "Figure_6_Analytical_Workflow.png"
)

workflow_tiff <- file.path(
  figure_dir,
  "Figure_6_Analytical_Workflow.tiff"
)

workflow_pdf <- file.path(
  figure_dir,
  "Figure_6_Analytical_Workflow.pdf"
)

png(
  filename = workflow_png,
  width = 16,
  height = 10,
  units = "in",
  res = 600,
  bg = "white"
)

draw_workflow()
invisible(dev.off())

tiff(
  filename = workflow_tiff,
  width = 16,
  height = 10,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

draw_workflow()
invisible(dev.off())

grDevices::cairo_pdf(
  filename = workflow_pdf,
  width = 16,
  height = 10,
  family = font_family
)

draw_workflow()
invisible(dev.off())

cat(
  "\nCorrected Figure 1 generated successfully:\n",
  workflow_png,
  "\n"
)















# ==============================================================================
# FIGURE 2
# MODEL SCREENING, CORRECTED NESTED-CV PERFORMANCE, AND MODEL COMPARISONS
# ==============================================================================

cat("\n--- Generating FIGURE 2: Corrected model performance and comparisons ---\n")

# ------------------------------------------------------------------------------
# 7. Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

screening_run_file <- file.path(
  stage1_dir,
  "Table_S0_Initial_7_Models_Run_Level_ROC.csv"
)

screening_summary_file <- file.path(
  stage1_dir,
  "Table_S0_Initial_7_Models_Numeric_Summary.csv"
)

corrected_performance_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Final_Model_Performance_Aggregated.csv"
)

corrected_outer_auc_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Outer_Resample_AUC.csv"
)

corrected_nb_test_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Nadeau_Bengio_Corrected_Tests.csv"
)

figure2_required_files <- c(
  screening_run_file,
  screening_summary_file,
  selected_models_file,
  corrected_performance_file,
  corrected_outer_auc_file,
  corrected_nb_test_file
)

figure2_missing_files <- figure2_required_files[
  !file.exists(figure2_required_files)
]

if (length(figure2_missing_files) > 0) {
  stop(
    "These Figure 2 source files are missing:\n",
    paste(figure2_missing_files, collapse = "\n")
  )
}

screening_runs <- read.csv(
  screening_run_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

screening_summary <- read.csv(
  screening_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

corrected_performance_table <- read.csv(
  corrected_performance_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

corrected_outer_auc <- read.csv(
  corrected_outer_auc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

corrected_nb_tests <- read.csv(
  corrected_nb_test_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

top3_models <- final_model_names

if (length(top3_models) != 3) {
  stop(
    "Exactly three final models are required. Found:\n",
    paste(top3_models, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 8. Validate columns
# ------------------------------------------------------------------------------

required_screening_columns <- c(
  "Run",
  "Model",
  "ROC"
)

missing_screening_columns <- setdiff(
  required_screening_columns,
  colnames(screening_runs)
)

if (length(missing_screening_columns) > 0) {
  stop(
    "Missing columns in the screening-run table:\n",
    paste(missing_screening_columns, collapse = ", ")
  )
}

required_summary_columns <- c(
  "Model",
  "MeanROC"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  colnames(screening_summary)
)

if (length(missing_summary_columns) > 0) {
  stop(
    "Missing columns in the screening summary:\n",
    paste(missing_summary_columns, collapse = ", ")
  )
}

corrected_metric_columns <- c(
  "AUC",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "F1",
  "MCC"
)

missing_corrected_metric_columns <- setdiff(
  c("Model", corrected_metric_columns),
  colnames(corrected_performance_table)
)

if (length(missing_corrected_metric_columns) > 0) {
  stop(
    "Missing columns in the corrected performance table:\n",
    paste(missing_corrected_metric_columns, collapse = ", ")
  )
}

required_outer_auc_columns <- c(
  "Model",
  "Resample",
  "AUC"
)

missing_outer_auc_columns <- setdiff(
  required_outer_auc_columns,
  colnames(corrected_outer_auc)
)

if (length(missing_outer_auc_columns) > 0) {
  stop(
    "Missing columns in the corrected outer-resample AUC table:\n",
    paste(missing_outer_auc_columns, collapse = ", ")
  )
}

required_nb_columns <- c(
  "Comparison",
  "Mean_AUC_difference",
  "Corrected_p_value",
  "N_paired_resamples"
)

missing_nb_columns <- setdiff(
  required_nb_columns,
  colnames(corrected_nb_tests)
)

if (length(missing_nb_columns) > 0) {
  stop(
    "Missing columns in the corrected Nadeau-Bengio table:\n",
    paste(missing_nb_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 9. Model colors
# ------------------------------------------------------------------------------

unknown_models <- setdiff(
  top3_models,
  names(default_model_colors)
)

if (length(unknown_models) > 0) {
  
  fallback_colors <- c(
    "#CC79A7",
    "#56B4E9",
    "#E69F00",
    "#000000"
  )
  
  default_model_colors <- c(
    default_model_colors,
    setNames(
      fallback_colors[
        seq_along(unknown_models)
      ],
      unknown_models
    )
  )
}

top3_colors <- default_model_colors[
  top3_models
]

# ------------------------------------------------------------------------------
# 10. Shared figure theme
# ------------------------------------------------------------------------------

publication_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    axis.text = element_text(
      color = "black",
      size = 10
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "gray30",
      margin = margin(b = 8)
    ),
    legend.text = element_text(
      size = 9.3
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------------------------
# 11. Panel A — Initial screening of seven classifiers
# ------------------------------------------------------------------------------

model_order <- screening_summary %>%
  arrange(MeanROC) %>%
  pull(Model) %>%
  as.character()

screening_runs <- screening_runs %>%
  mutate(
    Model = factor(
      Model,
      levels = model_order
    )
  )

all_model_colors <- setNames(
  rep("#BDBDBD", length(model_order)),
  model_order
)

matched_final_models <- intersect(
  top3_models,
  names(all_model_colors)
)

all_model_colors[
  matched_final_models
] <- top3_colors[
  matched_final_models
]

screening_minimum <- min(
  screening_runs$ROC,
  na.rm = TRUE
)

panel_a_lower_limit <- max(
  0,
  floor(
    (screening_minimum - 0.03) * 20
  ) / 20
)

p_screening <- ggplot(
  screening_runs,
  aes(
    x = Model,
    y = ROC,
    fill = Model,
    color = Model
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.42,
    width = 0.80,
    linewidth = 0.55
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    color = "#333333",
    linewidth = 0.50
  ) +
  geom_jitter(
    width = 0.07,
    height = 0,
    size = 1.75,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.6,
    fill = "black",
    color = "white",
    stroke = 0.65
  ) +
  scale_fill_manual(
    values = all_model_colors
  ) +
  scale_color_manual(
    values = all_model_colors
  ) +
  scale_y_continuous(
    limits = c(
      panel_a_lower_limit,
      1.015
    ),
    breaks = seq(
      panel_a_lower_limit,
      1,
      by = 0.05
    ),
    labels = scales::label_number(
      accuracy = 0.01
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  coord_flip(
    clip = "off"
  ) +
  labs(
    title = "A. Initial screening of seven classifiers",
    subtitle = paste0(
      "Ten runs of stratified 5-fold cross-validation; ",
      "diamonds indicate mean ROC-AUC"
    ),
    x = NULL,
    y = "Cross-validated ROC-AUC"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  publication_theme +
  theme(
    axis.text.y = element_text(
      size = 10.2
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.40
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 12. Panel B — Corrected nested-CV multi-metric performance
# ------------------------------------------------------------------------------

metric_labels <- c(
  "AUC" = "ROC-AUC",
  "Accuracy" = "Accuracy",
  "Sensitivity" = "Sensitivity",
  "Specificity" = "Specificity",
  "F1" = "F1 score",
  "MCC" = "MCC"
)

performance_long <- corrected_performance_table %>%
  filter(
    Model %in% top3_models
  ) %>%
  select(
    Model,
    all_of(corrected_metric_columns)
  ) %>%
  pivot_longer(
    cols = all_of(corrected_metric_columns),
    names_to = "Metric_code",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = unname(
      metric_labels[
        Metric_code
      ]
    ),
    Model = factor(
      Model,
      levels = top3_models
    ),
    Metric = factor(
      Metric,
      levels = rev(
        unname(metric_labels)
      )
    )
  )

p_metrics <- ggplot(
  performance_long,
  aes(
    x = Value,
    y = Metric,
    color = Model
  )
) +
  geom_point(
    position = position_dodge(
      width = 0.55
    ),
    size = 3.8
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Value
      )
    ),
    position = position_dodge(
      width = 0.55
    ),
    hjust = -0.42,
    size = 3.20,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = top3_colors
  ) +
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "B. Multi-metric performance",
    subtitle = paste0(
      "Nested repeated CV: outer 5-fold × 20; ",
      "inner 4-fold tuning"
    ),
    x = "Performance estimate",
    y = NULL,
    color = NULL
  ) +
  publication_theme +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(
      t = 4,
      b = 0
    ),
    panel.grid.major.y = element_line(
      color = "gray91",
      linewidth = 0.45
    ),
    panel.grid.major.x = element_line(
      color = "gray95",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 10.2
    )
  )

# ------------------------------------------------------------------------------
# 13. Panel C — Nadeau-Bengio corrected pairwise comparisons
# ------------------------------------------------------------------------------

compute_corrected_pairwise_auc <- function(
    outer_auc_table,
    model_1,
    model_2,
    k_folds = 5,
    confidence_level = 0.95
) {
  
  auc_1 <- outer_auc_table %>%
    filter(
      Model == model_1
    ) %>%
    transmute(
      Resample = as.character(Resample),
      AUC_1 = as.numeric(AUC)
    )
  
  auc_2 <- outer_auc_table %>%
    filter(
      Model == model_2
    ) %>%
    transmute(
      Resample = as.character(Resample),
      AUC_2 = as.numeric(AUC)
    )
  
  paired_auc <- inner_join(
    auc_1,
    auc_2,
    by = "Resample"
  )
  
  differences <- paired_auc$AUC_1 - paired_auc$AUC_2
  differences <- differences[
    is.finite(differences)
  ]
  
  n_resamples <- length(differences)
  
  if (n_resamples < 2) {
    stop(
      "Insufficient paired outer-resample AUCs for ",
      model_1,
      " versus ",
      model_2,
      "."
    )
  }
  
  mean_difference <- mean(
    differences
  )
  
  difference_variance <- var(
    differences
  )
  
  correction_factor <- (
    1 / n_resamples
  ) + (
    1 / (k_folds - 1)
  )
  
  if (
    !is.finite(difference_variance) ||
    difference_variance <= .Machine$double.eps
  ) {
    
    corrected_se <- 0
    
    corrected_t <- ifelse(
      abs(mean_difference) <= .Machine$double.eps,
      0,
      Inf
    )
    
    corrected_p_value <- ifelse(
      abs(mean_difference) <= .Machine$double.eps,
      1,
      0
    )
    
    ci_lower <- mean_difference
    ci_upper <- mean_difference
    
  } else {
    
    corrected_se <- sqrt(
      difference_variance *
        correction_factor
    )
    
    degrees_of_freedom <- n_resamples - 1
    
    corrected_t <- mean_difference /
      corrected_se
    
    corrected_p_value <- 2 * pt(
      -abs(corrected_t),
      df = degrees_of_freedom
    )
    
    alpha <- 1 - confidence_level
    
    critical_t <- qt(
      1 - alpha / 2,
      df = degrees_of_freedom
    )
    
    ci_lower <- mean_difference -
      critical_t * corrected_se
    
    ci_upper <- mean_difference +
      critical_t * corrected_se
  }
  
  data.frame(
    Model_1 = model_1,
    Model_2 = model_2,
    Comparison = paste(
      model_1,
      "vs",
      model_2
    ),
    Mean_AUC_difference = mean_difference,
    Corrected_SE = corrected_se,
    Corrected_t = corrected_t,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    Corrected_p_value = corrected_p_value,
    N_paired_resamples = n_resamples,
    stringsAsFactors = FALSE
  )
}

model_pairs <- combn(
  top3_models,
  2,
  simplify = FALSE
)

pairwise_statistics <- bind_rows(
  lapply(
    model_pairs,
    function(pair) {
      compute_corrected_pairwise_auc(
        outer_auc_table = corrected_outer_auc,
        model_1 = pair[1],
        model_2 = pair[2],
        k_folds = 5,
        confidence_level = 0.95
      )
    }
  )
)

# Preserve the exact corrected p-values and mean differences reported by 3B.
pairwise_statistics <- pairwise_statistics %>%
  select(
    -Corrected_p_value,
    -Mean_AUC_difference
  ) %>%
  left_join(
    corrected_nb_tests %>%
      select(
        Comparison,
        Mean_AUC_difference,
        Corrected_p_value,
        N_paired_resamples_reported =
          N_paired_resamples
      ),
    by = "Comparison"
  ) %>%
  mutate(
    P_label = case_when(
      Corrected_p_value < 0.001 ~
        "p < 0.001",
      TRUE ~ paste0(
        "p = ",
        sprintf(
          "%.3f",
          Corrected_p_value
        )
      )
    ),
    Display_label = paste0(
      P_label,
      ifelse(
        Corrected_p_value < 0.05,
        "; significant",
        "; ns"
      )
    ),
    Significance = ifelse(
      Corrected_p_value < 0.05,
      "Significant",
      "Not significant"
    ),
    Comparison = factor(
      Comparison,
      levels = rev(
        as.character(Comparison)
      )
    )
  )

maximum_ci_absolute <- max(
  abs(
    c(
      pairwise_statistics$CI_lower,
      pairwise_statistics$CI_upper
    )
  ),
  na.rm = TRUE
)

if (
  !is.finite(maximum_ci_absolute) ||
  maximum_ci_absolute < 0.01
) {
  maximum_ci_absolute <- 0.01
}

left_axis_limit <- -1.25 *
  maximum_ci_absolute

right_axis_limit <- 1.85 *
  maximum_ci_absolute

pairwise_statistics <- pairwise_statistics %>%
  mutate(
    Label_x = 1.12 *
      maximum_ci_absolute
  )

p_nadeau <- ggplot(
  pairwise_statistics,
  aes(
    x = Mean_AUC_difference,
    y = Comparison
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.65,
    color = "gray45"
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_upper,
      y = Comparison,
      yend = Comparison
    ),
    linewidth = 0.90,
    color = "#333333"
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_lower,
      y = as.numeric(Comparison) - 0.08,
      yend = as.numeric(Comparison) + 0.08
    ),
    linewidth = 0.75,
    color = "#333333"
  ) +
  geom_segment(
    aes(
      x = CI_upper,
      xend = CI_upper,
      y = as.numeric(Comparison) - 0.08,
      yend = as.numeric(Comparison) + 0.08
    ),
    linewidth = 0.75,
    color = "#333333"
  ) +
  geom_point(
    aes(
      fill = Significance
    ),
    shape = 21,
    size = 4.5,
    color = "black",
    stroke = 0.70
  ) +
  geom_text(
    aes(
      x = Label_x,
      label = Display_label
    ),
    hjust = 0,
    size = 3.30,
    family = font_family,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Significant" = "#B2182B",
      "Not significant" = "#737373"
    )
  ) +
  scale_x_continuous(
    limits = c(
      left_axis_limit,
      right_axis_limit
    ),
    labels = function(x) {
      sprintf("%+.2f", x)
    },
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "C. Pairwise ROC-AUC differences",
    subtitle = paste0(
      "Outer-resample AUC differences; bars show ",
      "Nadeau-Bengio-adjusted 95% CIs"
    ),
    x = paste0(
      "Mean ROC-AUC difference ",
      "(first model - second model)"
    ),
    y = NULL,
    fill = NULL
  ) +
  publication_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      size = 9.8
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.40
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 25,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------------------------
# 14. Assemble and save corrected Figure 2
# ------------------------------------------------------------------------------

figure_2 <- (
  p_screening |
    p_metrics |
    p_nadeau
) +
  plot_layout(
    widths = c(
      1.12,
      1.15,
      1.28
    )
  )

print(
  figure_2
)
save_publication_plot(
  plot = figure_2,
  filename_base =
    "Figure_7_Model_Selection_and_Comparison",
  width = 18,
  height = 6.25,
  dpi = 600
)

# ------------------------------------------------------------------------------
# 15. Save exact plotted data for reproducibility
# ------------------------------------------------------------------------------

write.csv(
  screening_runs,
  file.path(
    figure_dir,
    "Figure_2A_Corrected_Screening_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  performance_long,
  file.path(
    figure_dir,
    "Figure_2B_Corrected_Performance_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_statistics,
  file.path(
    figure_dir,
    "Figure_2C_Corrected_Nadeau_Bengio_Data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Figure 2 generated successfully.\n",
  "All outputs were saved in:\n",
  figure_dir,
  "\n"
)

cat(
  "\nCorrected nested-CV model performance:\n"
)

print(
  corrected_performance_table %>%
    select(
      Model,
      AUC,
      Accuracy,
      Sensitivity,
      Specificity,
      F1,
      MCC,
      TP,
      TN,
      FP,
      FN
    ),
  digits = 8
)

cat(
  "\nCorrected Nadeau-Bengio comparisons:\n"
)

print(
  pairwise_statistics %>%
    mutate(
      Comparison = as.character(
        Comparison
      )
    ) %>%
    select(
      Comparison,
      Mean_AUC_difference,
      CI_lower,
      CI_upper,
      Corrected_p_value,
      N_paired_resamples
    ),
  digits = 8
)

















# ==============================================================================
# FIGURE 3
# CORRECTED DIAGNOSTIC PERFORMANCE AND SAMPLE-LEVEL AGREEMENT
# ==============================================================================

cat("\n--- Generating FIGURE 3: Corrected diagnostic performance and agreement ---\n")

# ------------------------------------------------------------------------------
# 16. Packages required for Figure 3
# ------------------------------------------------------------------------------

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Install the pROC package first.")
}

library(pROC)

# ------------------------------------------------------------------------------
# 17. Input files from corrected 3B analysis
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

corrected_aggregated_prediction_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Aggregated_Predictions.csv"
)

corrected_performance_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Final_Model_Performance_Aggregated.csv"
)

figure3_required_files <- c(
  corrected_aggregated_prediction_file,
  corrected_performance_file
)

figure3_missing_files <- figure3_required_files[
  !file.exists(figure3_required_files)
]

if (length(figure3_missing_files) > 0) {
  stop(
    "These Figure 3 source files are missing:\n",
    paste(figure3_missing_files, collapse = "\n")
  )
}

corrected_predictions <- read.csv(
  corrected_aggregated_prediction_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

corrected_performance_table <- read.csv(
  corrected_performance_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# 18. Validate columns
# ------------------------------------------------------------------------------

required_prediction_columns <- c(
  "Model",
  "rowIndex",
  "Sample",
  "obs",
  "Mean_Probability_FVPTC",
  "Number_held_out_predictions"
)

missing_prediction_columns <- setdiff(
  required_prediction_columns,
  colnames(corrected_predictions)
)

if (length(missing_prediction_columns) > 0) {
  stop(
    "Missing columns in corrected aggregated predictions:\n",
    paste(missing_prediction_columns, collapse = ", ")
  )
}

required_performance_columns <- c(
  "Model",
  "AUC",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "F1",
  "MCC",
  "TP",
  "TN",
  "FP",
  "FN"
)

missing_performance_columns <- setdiff(
  required_performance_columns,
  colnames(corrected_performance_table)
)

if (length(missing_performance_columns) > 0) {
  stop(
    "Missing columns in corrected performance table:\n",
    paste(missing_performance_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 19. Model order and colors
# ------------------------------------------------------------------------------

if (!exists("final_model_names")) {
  
  selected_models_file <- file.path(
    stage2_old_dir,
    "Table_Final_Selected_Models_From_Stage1.csv"
  )
  
  selected_models_table <- read.csv(
    selected_models_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  final_model_names <- selected_models_table %>%
    arrange(Final_model_rank) %>%
    pull(Model) %>%
    as.character()
}

top3_models <- final_model_names

corrected_predictions <- corrected_predictions %>%
  filter(
    Model %in% top3_models
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = top3_models
    ),
    obs = factor(
      obs,
      levels = c(
        "FA",
        "FVPTC"
      )
    ),
    Mean_Probability_FVPTC = as.numeric(
      Mean_Probability_FVPTC
    ),
    pred = factor(
      ifelse(
        Mean_Probability_FVPTC >= 0.50,
        "FVPTC",
        "FA"
      ),
      levels = c(
        "FA",
        "FVPTC"
      )
    )
  )

corrected_performance_table <- corrected_performance_table %>%
  filter(
    Model %in% top3_models
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = top3_models
    )
  )

if (!exists("default_model_colors")) {
  default_model_colors <- c(
    "XGBoost" = "#0072B2",
    "Naive Bayes" = "#D55E00",
    "KNN" = "#009E73"
  )
}

model_colors <- default_model_colors[
  top3_models
]

if (any(is.na(model_colors))) {
  fallback_colors <- c(
    "#0072B2",
    "#D55E00",
    "#009E73"
  )
  
  model_colors <- setNames(
    fallback_colors[
      seq_along(top3_models)
    ],
    top3_models
  )
}

if (!exists("group_colors")) {
  group_colors <- c(
    "FA" = "#3B78A8",
    "FVPTC" = "#D9534F"
  )
}

agreement_colors <- c(
  "Unanimous correct" = "#1B9E77",
  "Majority correct" = "#66A61E",
  "Majority incorrect" = "#E6AB02",
  "Unanimous incorrect" = "#D95F02"
)

# ------------------------------------------------------------------------------
# 20. Theme
# ------------------------------------------------------------------------------

figure3_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    axis.text = element_text(
      color = "black",
      size = 9.8
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.3,
      color = "gray30",
      margin = margin(b = 8)
    ),
    legend.title = element_text(
      face = "bold",
      size = 9.5
    ),
    legend.text = element_text(
      size = 9.2
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------------------------
# 21. Panel A — Confusion matrices
# ------------------------------------------------------------------------------

confusion_data <- corrected_predictions %>%
  count(
    Model,
    obs,
    pred,
    name = "Count"
  ) %>%
  complete(
    Model,
    obs,
    pred,
    fill = list(
      Count = 0
    )
  ) %>%
  mutate(
    Cell = case_when(
      obs == "FA" & pred == "FA" ~ "TN",
      obs == "FA" & pred == "FVPTC" ~ "FP",
      obs == "FVPTC" & pred == "FA" ~ "FN",
      obs == "FVPTC" & pred == "FVPTC" ~ "TP",
      TRUE ~ ""
    ),
    Label = paste0(
      Cell,
      "\n",
      Count
    ),
    obs = factor(
      obs,
      levels = rev(
        c(
          "FA",
          "FVPTC"
        )
      )
    ),
    pred = factor(
      pred,
      levels = c(
        "FA",
        "FVPTC"
      )
    )
  )

p_confusion <- ggplot(
  confusion_data,
  aes(
    x = pred,
    y = obs,
    fill = Count
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.9
  ) +
  geom_text(
    aes(
      label = Label
    ),
    size = 4.1,
    fontface = "bold",
    family = font_family,
    color = "black"
  ) +
  facet_wrap(
    ~ Model,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2C7FB8"
  ) +
  labs(
    title = "A. Confusion matrices",
    subtitle = "Aggregated outer held-out predictions; threshold = 0.50",
    x = "Predicted class",
    y = "Observed class",
    fill = "Count"
  ) +
  figure3_theme +
  theme(
    strip.background = element_rect(
      fill = "gray95",
      color = "gray70"
    ),
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    panel.border = element_rect(
      fill = NA,
      color = "gray70",
      linewidth = 0.5
    ),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# 22. Panel B — Corrected ROC curves
# ------------------------------------------------------------------------------

roc_data <- bind_rows(
  lapply(
    top3_models,
    function(model_name) {
      
      model_data <- corrected_predictions %>%
        filter(
          as.character(Model) == model_name
        )
      
      roc_object <- pROC::roc(
        response = model_data$obs,
        predictor = model_data$Mean_Probability_FVPTC,
        levels = c(
          "FA",
          "FVPTC"
        ),
        direction = "<",
        quiet = TRUE
      )
      
      roc_coordinates <- as.data.frame(
        pROC::coords(
          roc_object,
          x = "all",
          ret = c(
            "specificity",
            "sensitivity"
          ),
          transpose = FALSE
        )
      )
      
      data.frame(
        Model = model_name,
        FPR = 1 - roc_coordinates$specificity,
        Sensitivity = roc_coordinates$sensitivity,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  mutate(
    Model = factor(
      Model,
      levels = top3_models
    )
  )

roc_summary <- corrected_performance_table %>%
  transmute(
    Model = as.character(Model),
    AUC = as.numeric(AUC),
    Legend_label = paste0(
      Model,
      " (AUC = ",
      sprintf(
        "%.3f",
        AUC
      ),
      ")"
    )
  )

roc_label_map <- setNames(
  roc_summary$Legend_label,
  roc_summary$Model
)

roc_color_map <- setNames(
  model_colors[
    roc_summary$Model
  ],
  roc_summary$Legend_label
)

roc_data <- roc_data %>%
  mutate(
    Legend_label = factor(
      roc_label_map[
        as.character(Model)
      ],
      levels = roc_summary$Legend_label
    )
  )

p_roc <- ggplot(
  roc_data,
  aes(
    x = FPR,
    y = Sensitivity,
    color = Legend_label
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray65",
    linewidth = 0.65
  ) +
  geom_line(
    linewidth = 1.15
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  scale_color_manual(
    values = roc_color_map
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      1,
      by = 0.25
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.25
    )
  ) +
  labs(
    title = "B. ROC curves",
    subtitle = "Curves from aggregated outer held-out probabilities",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  figure3_theme +
  theme(
    legend.position = c(
      0.63,
      0.18
    ),
    legend.background = element_rect(
      fill = grDevices::adjustcolor(
        "white",
        alpha.f = 0.85
      ),
      color = "gray80"
    ),
    panel.grid.major = element_line(
      color = "gray92",
      linewidth = 0.35
    )
  )

# ------------------------------------------------------------------------------
# 23. Panel C — Sample-level probability heatmap
# ------------------------------------------------------------------------------

sample_level_summary <- corrected_predictions %>%
  group_by(
    rowIndex,
    Sample,
    obs
  ) %>%
  summarise(
    Mean_probability_across_models = mean(
      Mean_Probability_FVPTC,
      na.rm = TRUE
    ),
    Correct_votes = sum(
      pred == obs,
      na.rm = TRUE
    ),
    Incorrect_votes = sum(
      pred != obs,
      na.rm = TRUE
    ),
    N_models = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Agreement_category = case_when(
      Correct_votes == N_models ~ "Unanimous correct",
      Correct_votes >= 2 ~ "Majority correct",
      Correct_votes == 1 ~ "Majority incorrect",
      Correct_votes == 0 ~ "Unanimous incorrect",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    obs,
    Mean_probability_across_models
  ) %>%
  mutate(
    Sample_display = paste0(
      "S",
      sprintf(
        "%02d",
        row_number()
      )
    ),
    Sample_display = factor(
      Sample_display,
      levels = Sample_display
    )
  )

heatmap_prediction_data <- corrected_predictions %>%
  left_join(
    sample_level_summary %>%
      select(
        rowIndex,
        Sample_display
      ),
    by = "rowIndex"
  ) %>%
  transmute(
    rowIndex,
    Sample,
    Sample_display,
    Row = as.character(Model),
    Value = Mean_Probability_FVPTC,
    obs
  )

observed_class_row <- sample_level_summary %>%
  transmute(
    rowIndex,
    Sample,
    Sample_display,
    Row = "Observed class",
    Value = ifelse(
      obs == "FVPTC",
      1,
      0
    ),
    obs
  )

heatmap_data <- bind_rows(
  heatmap_prediction_data,
  observed_class_row
) %>%
  mutate(
    Row = factor(
      Row,
      levels = c(
        "Observed class",
        rev(top3_models)
      )
    )
  )

class_boundary <- sample_level_summary %>%
  count(
    obs
  ) %>%
  arrange(
    obs
  )

fa_count <- class_boundary$n[
  class_boundary$obs == "FA"
]

p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Sample_display,
    y = Row,
    fill = Value
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = fa_count + 0.5,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.60
  ) +
  scale_fill_gradient2(
    low = "#3B78A8",
    mid = "white",
    high = "#D9534F",
    midpoint = 0.50,
    limits = c(0, 1),
    breaks = c(
      0,
      0.5,
      1
    ),
    labels = c(
      "0",
      "0.5",
      "1"
    )
  ) +
  labs(
    title = "C. Sample-level probability landscape",
    subtitle = "Samples ordered by observed class and mean FVPTC probability",
    x = "Samples",
    y = NULL,
    fill = "FVPTC\nprobability"
  ) +
  figure3_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# 24. Panel D — Model agreement summary
# ------------------------------------------------------------------------------

agreement_summary <- sample_level_summary %>%
  count(
    Agreement_category,
    name = "Sample_count"
  ) %>%
  complete(
    Agreement_category = names(agreement_colors),
    fill = list(
      Sample_count = 0
    )
  ) %>%
  mutate(
    Agreement_category = factor(
      Agreement_category,
      levels = names(agreement_colors)
    )
  )

p_agreement <- ggplot(
  agreement_summary,
  aes(
    x = Agreement_category,
    y = Sample_count,
    fill = Agreement_category
  )
) +
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Sample_count
    ),
    vjust = -0.35,
    size = 3.7,
    fontface = "bold",
    family = font_family
  ) +
  scale_fill_manual(
    values = agreement_colors
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        agreement_summary$Sample_count,
        na.rm = TRUE
      ) + 4
    ),
    breaks = scales::pretty_breaks(
      n = 5
    )
  ) +
  labs(
    title = "D. Sample-level agreement",
    subtitle = "Agreement across XGBoost, Naive Bayes, and KNN",
    x = NULL,
    y = "Number of samples"
  ) +
  figure3_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 9.2
    ),
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.35
    )
  )

# ------------------------------------------------------------------------------
# 25. Assemble and save corrected Figure 3
# ------------------------------------------------------------------------------

figure_3 <- (
  p_confusion |
    p_roc
) / (
  p_heatmap |
    p_agreement
) +
  plot_layout(
    widths = c(
      1.35,
      1.00
    ),
    heights = c(
      1.00,
      1.08
    )
  )

print(
  figure_3
)

save_publication_plot(
  plot = figure_3,
  filename_base =
    "Figure_8_Diagnostic_Performance_and_Agreement",
  width = 18,
  height = 10.8,
  dpi = 600
)

# ------------------------------------------------------------------------------
# 26. Save exact plotted data for reproducibility
# ------------------------------------------------------------------------------

write.csv(
  confusion_data,
  file.path(
    figure_dir,
    "Figure_3A_Corrected_Confusion_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  roc_summary,
  file.path(
    figure_dir,
    "Figure_3B_Corrected_ROC_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  heatmap_data,
  file.path(
    figure_dir,
    "Figure_3C_Corrected_Probability_Heatmap_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_level_summary,
  file.path(
    figure_dir,
    "Figure_3C_Corrected_Sample_Level_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  agreement_summary,
  file.path(
    figure_dir,
    "Figure_3D_Corrected_Agreement_Summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Figure 3 generated successfully.\n",
  "All outputs were saved in:\n",
  figure_dir,
  "\n"
)

cat(
  "\nCorrected aggregated performance used in Figure 3:\n"
)

print(
  corrected_performance_table %>%
    select(
      Model,
      AUC,
      Accuracy,
      Sensitivity,
      Specificity,
      F1,
      MCC,
      TP,
      TN,
      FP,
      FN
    ),
  digits = 8
)

cat(
  "\nSample-level agreement summary:\n"
)

print(
  agreement_summary,
  row.names = FALSE
)




















# ==============================================================================
# FIGURE 4
# CORRECTED INTEGRATED BIOMARKER PRIORITIZATION
# ==============================================================================

cat("\n--- Generating FIGURE 4: Corrected integrated biomarker prioritization ---\n")

# ------------------------------------------------------------------------------
# 27. Input files from corrected 3B analysis
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

corrected_pi_summary_file <- file.path(
  stage2_corrected_dir,
  "Table_CrossFitted_PI_Summary_All_Models.csv"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

corrected_selected_panels_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

figure4_required_files <- c(
  corrected_pi_summary_file,
  corrected_consensus_file,
  corrected_selected_panels_file
)

figure4_missing_files <- figure4_required_files[
  !file.exists(figure4_required_files)
]

if (length(figure4_missing_files) > 0) {
  stop(
    "These Figure 4 source files are missing:\n",
    paste(figure4_missing_files, collapse = "\n")
  )
}

pi_summary <- read.csv(
  corrected_pi_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

consensus_table <- read.csv(
  corrected_consensus_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

selected_panels_corrected <- read.csv(
  corrected_selected_panels_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
# 28. Validate columns
# ------------------------------------------------------------------------------

required_pi_columns <- c(
  "Model",
  "Gene",
  "Mean",
  "SD",
  "Positive_fraction",
  "Positive_permutation_fraction"
)

missing_pi_columns <- setdiff(
  required_pi_columns,
  colnames(pi_summary)
)

if (length(missing_pi_columns) > 0) {
  stop(
    "Missing columns in corrected PI summary:\n",
    paste(missing_pi_columns, collapse = ", ")
  )
}

required_consensus_columns <- c(
  "Consensus_rank",
  "Gene",
  "Normalized_consensus_score",
  "Score_drop",
  "Positive_models",
  "Mean_delta_AUC_across_models",
  "Mean_positive_frequency_across_models",
  "Expression_direction",
  "Wilcoxon_BH_FDR",
  "Expression_confirmed"
)

missing_consensus_columns <- setdiff(
  required_consensus_columns,
  colnames(consensus_table)
)

if (length(missing_consensus_columns) > 0) {
  stop(
    "Missing columns in corrected consensus table:\n",
    paste(missing_consensus_columns, collapse = ", ")
  )
}

if (!all(c("Panel", "Genes") %in% colnames(selected_panels_corrected))) {
  stop(
    "The corrected selected-panel table must contain ",
    "'Panel' and 'Genes'."
  )
}

# ------------------------------------------------------------------------------
# 29. Extract corrected panels
# ------------------------------------------------------------------------------

extract_corrected_genes <- function(panel_name) {
  
  selected_row <- selected_panels_corrected %>%
    filter(
      Panel == panel_name
    )
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

primary_panel <- extract_corrected_genes(
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- extract_corrected_genes(
  "Corrected_secondary_candidates"
)

expression_only_genes <- extract_corrected_genes(
  "Expression_confirmed_not_ML_supported"
)

zero_consensus_genes <- extract_corrected_genes(
  "Zero_normalized_consensus_genes"
)

if (length(primary_panel) == 0) {
  stop("No corrected primary panel genes were found.")
}

# ------------------------------------------------------------------------------
# 30. Model order, gene order, and colors
# ------------------------------------------------------------------------------

if (!exists("final_model_names")) {
  
  selected_models_file <- file.path(
    stage2_old_dir,
    "Table_Final_Selected_Models_From_Stage1.csv"
  )
  
  selected_models_table <- read.csv(
    selected_models_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  final_model_names <- selected_models_table %>%
    arrange(Final_model_rank) %>%
    pull(Model) %>%
    as.character()
}

top3_models <- final_model_names

consensus_table <- consensus_table %>%
  arrange(
    Consensus_rank
  )

gene_order <- as.character(
  consensus_table$Gene
)

if (!exists("default_model_colors")) {
  default_model_colors <- c(
    "XGBoost" = "#0072B2",
    "Naive Bayes" = "#D55E00",
    "KNN" = "#009E73"
  )
}

model_colors <- default_model_colors[
  top3_models
]

if (any(is.na(model_colors))) {
  fallback_colors <- c(
    "#0072B2",
    "#D55E00",
    "#009E73"
  )
  
  model_colors <- setNames(
    fallback_colors[
      seq_along(top3_models)
    ],
    top3_models
  )
}

panel_colors <- c(
  "Primary corrected panel" = "#1B9E77",
  "Secondary candidate" = "#E6AB02",
  "Weak/model-limited secondary" = "#D95F02",
  "Not ML-supported" = "#BDBDBD"
)

# ------------------------------------------------------------------------------
# 31. Shared Figure 4 theme
# ------------------------------------------------------------------------------

figure4_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    axis.text = element_text(
      color = "black",
      size = 9.8
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 9.3,
      color = "gray30",
      margin = margin(b = 8)
    ),
    legend.title = element_text(
      face = "bold",
      size = 9.4
    ),
    legend.text = element_text(
      size = 9.0
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------------------------
# 32. Panel A — Cross-fitted permutation importance matrix
# ------------------------------------------------------------------------------

pi_matrix <- pi_summary %>%
  filter(
    Model %in% top3_models,
    Gene %in% gene_order
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = top3_models
    ),
    Gene = factor(
      Gene,
      levels = rev(gene_order)
    ),
    Mean = as.numeric(Mean),
    Mean_label = sprintf(
      "%.3f",
      Mean
    )
  )

pi_limit <- max(
  abs(pi_matrix$Mean),
  na.rm = TRUE
)

if (
  !is.finite(pi_limit) ||
  pi_limit <= 0
) {
  pi_limit <- 0.01
}

p_pi <- ggplot(
  pi_matrix,
  aes(
    x = Model,
    y = Gene,
    fill = Mean
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Mean_label
    ),
    size = 4.0,              # سایز بزرگتر شد
    fontface = "bold",       # بُولد شد
    family = font_family,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B78A8",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -pi_limit,
      pi_limit
    ),
    name = "Mean\nΔAUC"
  ) +
  labs(
    title = "A. Cross-fitted permutation importance",
    subtitle = "Mean ΔAUC on outer held-out samples",
    x = NULL,
    y = NULL
  ) +
  figure4_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# 33. Panel B — Model-normalized Borda consensus ranking
# ------------------------------------------------------------------------------

consensus_plot <- consensus_table %>%
  mutate(
    Gene = factor(
      Gene,
      levels = rev(gene_order)
    ),
    Evidence_group = case_when(
      as.character(Gene) %in% primary_panel ~
        "Primary corrected panel",
      as.character(Gene) == "SPP1" ~
        "Weak/model-limited secondary",
      as.character(Gene) %in% secondary_candidates ~
        "Secondary candidate",
      TRUE ~ "Not ML-supported"
    ),
    Evidence_group = factor(
      Evidence_group,
      levels = names(panel_colors)
    ),
    Score_label = sprintf(
      "%.3f",
      Normalized_consensus_score
    ),
    Rank_label = paste0(
      "#",
      Consensus_rank
    )
  )

p_consensus <- ggplot(
  consensus_plot,
  aes(
    x = Normalized_consensus_score,
    y = Gene,
    fill = Evidence_group
  )
) +
  geom_col(
    width = 0.70,
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = Score_label
    ),
    hjust = -0.10,
    size = 4.2,              # سایز بزرگتر شد
    fontface = "bold",       # بُولد شد
    family = font_family
  ) +
  scale_fill_manual(
    values = panel_colors,
    labels = c(
      "Primary corrected panel" = "Primary panel",
      "Secondary candidate" = "Secondary candidate",
      "Weak/model-limited secondary" = "Weak/model-limited secondary",
      "Not ML-supported" = "Not ML-supported"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(
        consensus_plot$Normalized_consensus_score,
        na.rm = TRUE
      ) * 1.25
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "B. Model-normalized consensus ranking",
    subtitle = "Each classifier contributes equal total Borda weight",
    x = "Normalized consensus score",
    y = NULL,
    fill = NULL
  ) +
  figure4_theme +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.35
    )
  )

# ------------------------------------------------------------------------------
# 34. Panel C — Dynamic score-drop panel selection
# ------------------------------------------------------------------------------

primary_cutoff <- max(
  consensus_table$Consensus_rank[
    consensus_table$Gene %in% primary_panel
  ],
  na.rm = TRUE
)

elbow_plot <- consensus_table %>%
  mutate(
    Panel_status = ifelse(
      Gene %in% primary_panel,
      "Primary corrected panel",
      "Secondary or weak"
    ),
    Panel_status = factor(
      Panel_status,
      levels = c(
        "Primary corrected panel",
        "Secondary or weak"
      )
    )
  )

p_elbow <- ggplot(
  elbow_plot,
  aes(
    x = Consensus_rank,
    y = Normalized_consensus_score
  )
) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = primary_cutoff + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#1B9E77",
    alpha = 0.08
  ) +
  geom_line(
    linewidth = 0.85,
    color = "gray35"
  ) +
  geom_point(
    aes(
      fill = Panel_status
    ),
    shape = 21,
    size = 4.1,
    color = "black",
    stroke = 0.65
  ) +
  geom_text(
    aes(
      label = Gene
    ),
    nudge_y = 0.012,
    size = 4.2,              # سایز بزرگتر شد
    fontface = "bold",       # بُولد شد
    family = font_family,
    check_overlap = TRUE
  ) +
  geom_vline(
    xintercept = primary_cutoff + 0.5,
    linetype = "dashed",
    linewidth = 0.70,
    color = "gray35"
  ) +
  scale_fill_manual(
    values = c(
      "Primary corrected panel" = "#1B9E77",
      "Secondary or weak" = "#BDBDBD"
    )
  ) +
  scale_x_continuous(
    breaks = consensus_table$Consensus_rank
  ) +
  labs(
    title = "C. Dynamic score-drop panel selection",
    subtitle = paste0(
      "Primary panel: ",
      paste(primary_panel, collapse = ", ")
    ),
    x = "Consensus rank",
    y = "Normalized consensus score",
    fill = NULL
  ) +
  figure4_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.35
    )
  )

# ------------------------------------------------------------------------------
# 35. Panel D — Integrated expression and ML evidence
# ------------------------------------------------------------------------------

evidence_summary <- consensus_table %>%
  mutate(
    Evidence_group = case_when(
      Gene %in% primary_panel ~
        "Primary corrected panel",
      Gene == "SPP1" ~
        "Weak/model-limited secondary",
      Gene %in% secondary_candidates ~
        "Secondary candidate",
      TRUE ~ "Not ML-supported"
    ),
    Evidence_group = factor(
      Evidence_group,
      levels = names(panel_colors)
    ),
    Gene = factor(
      Gene,
      levels = rev(gene_order)
    ),
    ML_models_label = paste0(
      Positive_models,
      "/3 models"
    ),
    Stability_label = sprintf(
      "%.2f",
      Mean_positive_frequency_across_models
    ),
    FDR_label = paste0(
      "FDR=",
      formatC(
        Wilcoxon_BH_FDR,
        format = "e",
        digits = 1
      )
    )
  )

evidence_long <- bind_rows(
  evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "ML-positive models",
      Label = ML_models_label
    ),
  evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "Mean positive\nresample frequency",
      Label = Stability_label
    ),
  evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "Expression\nBH-FDR",
      Label = FDR_label
    )
) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "ML-positive models",
        "Mean positive\nresample frequency",
        "Expression\nBH-FDR"
      )
    )
  )

p_evidence <- ggplot(
  evidence_long,
  aes(
    x = Metric,
    y = Gene,
    fill = Evidence_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Label
    ),
    size = 3.8,              # سایز بزرگتر شد
    fontface = "bold",       # بُولد شد
    family = font_family,
    color = "black",
    lineheight = 0.90
  ) +
  scale_fill_manual(
    values = panel_colors,
    labels = c(
      "Primary corrected panel" = "Primary panel",
      "Secondary candidate" = "Secondary candidate",
      "Weak/model-limited secondary" = "Weak/model-limited secondary",
      "Not ML-supported" = "Not ML-supported"
    ),
    drop = FALSE
  ) +
  labs(
    title = "D. Integrated expression and ML evidence",
    subtitle = "Expression support uses BH-FDR adjustment; ML support uses cross-fitted PI",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  figure4_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      face = "bold",
      size = 9.0
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

# ------------------------------------------------------------------------------
# 36. Assemble and save corrected Figure 4
# ------------------------------------------------------------------------------

figure_4 <- (
  p_pi |
    p_consensus
) / (
  p_elbow |
    p_evidence
) +
  plot_layout(
    heights = c(
      1.05,
      1.00
    ),
    widths = c(
      1.05,
      1.00
    )
  )

print(
  figure_4
)

save_publication_plot(
  plot = figure_4,
  filename_base =
    "Figure_9_Biomarker_Prioritization",
  width = 18,
  height = 12.3,
  dpi = 600
)

# ------------------------------------------------------------------------------
# 37. Save exact plotted data for reproducibility
# ------------------------------------------------------------------------------

write.csv(
  pi_matrix,
  file.path(
    figure_dir,
    "Figure_4A_Corrected_PI_Matrix_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  consensus_plot,
  file.path(
    figure_dir,
    "Figure_4B_Corrected_Consensus_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  elbow_plot,
  file.path(
    figure_dir,
    "Figure_4C_Corrected_Elbow_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  evidence_long,
  file.path(
    figure_dir,
    "Figure_4D_Corrected_Evidence_Data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Figure 4 generated successfully.\n",
  "All outputs were saved in:\n",
  figure_dir,
  "\n"
)

cat(
  "\nCorrected primary panel:\n",
  paste(primary_panel, collapse = ", "),
  "\n"
)

cat(
  "\nCorrected secondary candidates:\n",
  paste(secondary_candidates, collapse = ", "),
  "\n"
)

cat(
  "\nCorrected normalized consensus ranking:\n"
)

print(
  consensus_table %>%
    select(
      Consensus_rank,
      Gene,
      Normalized_consensus_score,
      Score_drop,
      Positive_models,
      Mean_delta_AUC_across_models,
      Mean_positive_frequency_across_models,
      Expression_direction,
      Wilcoxon_BH_FDR
    ),
  digits = 8
)
































# ==============================================================================
# FIGURE 5
# CORRECTED EXPRESSION LANDSCAPE AND UNIVARIATE EVIDENCE
# ==============================================================================

cat("\n--- Generating FIGURE 5: Corrected expression landscape and univariate evidence ---\n")

# ------------------------------------------------------------------------------
# 38. Packages required for Figure 5
# ------------------------------------------------------------------------------

required_figure5_packages <- c(
  "pROC",
  "pheatmap",
  "RColorBrewer",
  "tibble"
)

missing_figure5_packages <- required_figure5_packages[
  !vapply(
    required_figure5_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_figure5_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_figure5_packages, collapse = ", ")
  )
}

library(pROC)
library(pheatmap)
library(RColorBrewer)
library(tibble)

if (!exists("read_required_csv")) {
  
  read_required_csv <- function(
    path,
    required_columns = NULL
  ) {
    
    if (!file.exists(path)) {
      stop(
        "Required file was not found:\n",
        path
      )
    }
    
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!is.null(required_columns)) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(output)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Missing columns in:\n",
          path,
          "\nMissing: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }
    
    output
  }
}

if (!exists("save_publication_plot")) {
  
  save_publication_plot <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
  ) {
    
    output_base <- file.path(
      figure_dir,
      filename_base
    )
    
    ggsave(
      filename = paste0(output_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".tiff"),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".pdf"),
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

format_fdr <- function(x) {
  
  ifelse(
    is.na(x),
    "FDR = NA",
    ifelse(
      x < 0.001,
      "FDR < 0.001",
      paste0(
        "FDR = ",
        sprintf("%.3f", x)
      )
    )
  )
}

# ------------------------------------------------------------------------------
# 39. Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

expression_file <- file.path(
  project_dir,
  "final_data_clean.rds"
)

sample_label_file <- file.path(
  stage2_old_dir,
  "Report_Sample_Label_Order_Used_By_ML.csv"
)

corrected_prediction_file <- file.path(
  stage2_corrected_dir,
  "Table_NestedCV_Aggregated_Predictions.csv"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

corrected_univariate_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Univariate_Expression_Wilcoxon.csv"
)

corrected_selected_panels_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

figure5_required_files <- c(
  expression_file,
  sample_label_file,
  corrected_prediction_file,
  corrected_consensus_file,
  corrected_univariate_file,
  corrected_selected_panels_file
)

figure5_missing_files <- figure5_required_files[
  !file.exists(figure5_required_files)
]

if (length(figure5_missing_files) > 0) {
  stop(
    "These Figure 5 source files are missing:\n",
    paste(figure5_missing_files, collapse = "\n")
  )
}

final_data <- readRDS(
  expression_file
)

if (
  !is.matrix(final_data) &&
  !is.data.frame(final_data)
) {
  stop(
    "final_data_clean.rds must contain a matrix or data frame."
  )
}

final_data <- as.matrix(
  final_data
)

sample_labels <- read_required_csv(
  sample_label_file,
  required_columns = c(
    "Sample_order",
    "Sample_name",
    "Label"
  )
) %>%
  arrange(
    Sample_order
  )

corrected_predictions <- read_required_csv(
  corrected_prediction_file,
  required_columns = c(
    "Model",
    "rowIndex",
    "Sample",
    "obs",
    "Mean_Probability_FVPTC"
  )
)

consensus_table <- read_required_csv(
  corrected_consensus_file,
  required_columns = c(
    "Consensus_rank",
    "Gene",
    "Normalized_consensus_score",
    "Expression_direction",
    "Wilcoxon_BH_FDR"
  )
)

univariate_gene_report <- read_required_csv(
  corrected_univariate_file,
  required_columns = c(
    "Gene",
    "Median_FA",
    "Median_FVPTC",
    "Expression_direction",
    "Wilcoxon_BH_FDR"
  )
)

selected_panels_corrected <- read_required_csv(
  corrected_selected_panels_file,
  required_columns = c(
    "Panel",
    "Genes"
  )
)

# ------------------------------------------------------------------------------
# 40. Extract corrected primary and secondary panels
# ------------------------------------------------------------------------------

extract_corrected_genes <- function(panel_name) {
  
  selected_row <- selected_panels_corrected %>%
    filter(
      Panel == panel_name
    )
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

primary_panel <- extract_corrected_genes(
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- extract_corrected_genes(
  "Corrected_secondary_candidates"
)

if (length(primary_panel) == 0) {
  stop(
    "No corrected primary-panel genes were found."
  )
}

gene_list <- consensus_table %>%
  arrange(
    Consensus_rank
  ) %>%
  pull(
    Gene
  ) %>%
  as.character()

missing_expression_genes <- setdiff(
  gene_list,
  rownames(final_data)
)

if (length(missing_expression_genes) > 0) {
  stop(
    "These candidate genes are absent from final_data_clean.rds:\n",
    paste(missing_expression_genes, collapse = ", ")
  )
}

missing_primary_expression_genes <- setdiff(
  primary_panel,
  rownames(final_data)
)

if (length(missing_primary_expression_genes) > 0) {
  stop(
    "These corrected primary-panel genes are absent from final_data_clean.rds:\n",
    paste(missing_primary_expression_genes, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 41. Sample labels and corrected sample ordering
# ------------------------------------------------------------------------------

if (!identical(
  as.character(sample_labels$Sample_name),
  colnames(final_data)
)) {
  stop(
    "Sample order mismatch between final_data_clean.rds and ",
    "Report_Sample_Label_Order_Used_By_ML.csv."
  )
}

if (nrow(sample_labels) != 43) {
  stop(
    "The expected integrated dataset contains 43 samples, but ",
    nrow(sample_labels),
    " samples were loaded."
  )
}

batch_vector <- c(
  rep("GSE29315", 30),
  rep("GSE54958", 13)
)

sample_labels <- sample_labels %>%
  mutate(
    Label = factor(
      Label,
      levels = c(
        "FA",
        "FVPTC"
      )
    ),
    Batch = factor(
      batch_vector,
      levels = c(
        "GSE29315",
        "GSE54958"
      )
    )
  )

if (any(is.na(sample_labels$Label))) {
  stop(
    "At least one sample label is not FA or FVPTC."
  )
}

mean_probability_by_sample <- corrected_predictions %>%
  group_by(
    Sample
  ) %>%
  summarise(
    Mean_model_probability = mean(
      Mean_Probability_FVPTC,
      na.rm = TRUE
    ),
    Number_models = n_distinct(
      Model
    ),
    .groups = "drop"
  )

sample_level_order <- sample_labels %>%
  transmute(
    Sample = as.character(
      Sample_name
    ),
    Label,
    Batch
  ) %>%
  left_join(
    mean_probability_by_sample,
    by = "Sample"
  )

if (any(is.na(sample_level_order$Mean_model_probability))) {
  stop(
    "At least one sample could not be matched to corrected aggregated predictions."
  )
}

sample_level_order <- sample_level_order %>%
  arrange(
    Label,
    Mean_model_probability
  ) %>%
  mutate(
    Display_order = row_number()
  )

sample_expression_order <- sample_level_order$Sample

# ------------------------------------------------------------------------------
# 42. Construct sample-by-gene expression data
# ------------------------------------------------------------------------------

gene_data <- as.data.frame(
  t(
    final_data[
      gene_list,
      ,
      drop = FALSE
    ]
  ),
  stringsAsFactors = FALSE
)

gene_data$Sample <- rownames(
  gene_data
)

gene_data <- gene_data %>%
  left_join(
    sample_labels %>%
      transmute(
        Sample = as.character(
          Sample_name
        ),
        label = Label,
        Batch
      ),
    by = "Sample"
  )

if (any(is.na(gene_data$label))) {
  stop(
    "At least one expression sample could not be matched to its label."
  )
}

# ------------------------------------------------------------------------------
# 43. Figure 5 theme
# ------------------------------------------------------------------------------

figure5_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.5,
      color = "gray30",
      margin = margin(
        b = 8
      )
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 8
    )
  )

if (!exists("group_colors")) {
  group_colors <- c(
    "FA" = "#3B78A8",
    "FVPTC" = "#D9534F"
  )
}

batch_colors <- c(
  "GSE29315" = "#56B4E9",
  "GSE54958" = "#CC79A7"
)

# ==============================================================================
# PANEL A
# STANDARDIZED EXPRESSION LANDSCAPE OF THE EIGHT CANDIDATE GENES
# ==============================================================================

expression_gene_order <- gene_list

expression_matrix <- final_data[
  expression_gene_order,
  sample_expression_order,
  drop = FALSE
]

expression_matrix_scaled <- t(
  scale(
    t(
      expression_matrix
    )
  )
)

expression_matrix_scaled[
  !is.finite(expression_matrix_scaled)
] <- 0

expression_annotation <- sample_level_order %>%
  arrange(
    Display_order
  ) %>%
  transmute(
    Class = factor(
      Label,
      levels = c(
        "FA",
        "FVPTC"
      )
    ),
    Dataset = factor(
      Batch,
      levels = c(
        "GSE29315",
        "GSE54958"
      )
    )
  )

rownames(expression_annotation) <- sample_expression_order

expression_annotation_colors <- list(
  Class = group_colors,
  Dataset = batch_colors
)

class_gap <- sum(
  expression_annotation$Class == "FA"
)

heatmap_object <- pheatmap(
  expression_matrix_scaled,
  color = colorRampPalette(
    rev(
      RColorBrewer::brewer.pal(
        11,
        "RdBu"
      )
    )
  )(255),
  breaks = seq(
    -3,
    3,
    length.out = 256
  ),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = expression_annotation,
  annotation_colors = expression_annotation_colors,
  gaps_col = class_gap,
  fontsize = 10,
  fontsize_row = 10,
  main = paste0(
    "A. Standardized expression of the eight ",
    "candidate genes"
  ),
  silent = TRUE
)

p_expression_heatmap <- wrap_elements(
  full = heatmap_object$gtable
)

# ==============================================================================
# PANEL B
# EXPRESSION DISTRIBUTIONS OF THE CORRECTED PRIMARY PANEL
# ==============================================================================

primary_expression_long <- gene_data %>%
  select(
    Sample,
    label,
    all_of(primary_panel)
  ) %>%
  pivot_longer(
    cols = all_of(primary_panel),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = primary_panel
    )
  )

primary_fdr_labels <- primary_expression_long %>%
  group_by(
    Gene
  ) %>%
  summarise(
    expression_min = min(
      Expression,
      na.rm = TRUE
    ),
    expression_max = max(
      Expression,
      na.rm = TRUE
    ),
    expression_range =
      expression_max -
      expression_min,
    y_position =
      expression_max +
      ifelse(
        expression_range > 0,
        0.12 * expression_range,
        0.25
      ),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = as.character(Gene)
  ) %>%
  left_join(
    univariate_gene_report %>%
      filter(
        Gene %in% primary_panel
      ) %>%
      transmute(
        Gene = as.character(Gene),
        Wilcoxon_BH_FDR = as.numeric(
          Wilcoxon_BH_FDR
        )
      ),
    by = "Gene"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = primary_panel
    ),
    FDR_label = format_fdr(
      Wilcoxon_BH_FDR
    )
  )

p_primary_expression <- ggplot(
  primary_expression_long,
  aes(
    x = label,
    y = Expression,
    fill = label,
    color = label
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.28,
    width = 0.88,
    linewidth = 0.55
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.92,
    linewidth = 0.50
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.35,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  geom_text(
    data = primary_fdr_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = FDR_label
    ),
    inherit.aes = FALSE,
    family = font_family,
    size = 2.75
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    title = "B. Expression distributions of the primary panel",
    subtitle = paste0(
      "Points represent individual samples; significance is ",
      "BH-adjusted across the eight candidate genes"
    ),
    x = NULL,
    y = "Expression level (log2)"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  figure5_theme +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    panel.grid.major.y = element_line(
      color = "gray93",
      linewidth = 0.35
    )
  )

# ==============================================================================
# PANEL C
# DESCRIPTIVE DIRECTION-CORRECTED SINGLE-GENE ROC-AUC
# ==============================================================================

single_gene_roc_summary <- bind_rows(
  lapply(
    gene_list,
    function(gene_name) {
      
      expression_values <- as.numeric(
        gene_data[[gene_name]]
      )
      
      direction_result <- univariate_gene_report %>%
        filter(
          Gene == gene_name
        ) %>%
        pull(
          Expression_direction
        )
      
      direction_label <- if (
        length(direction_result) > 0 &&
        !is.na(direction_result[1])
      ) {
        as.character(direction_result[1])
      } else {
        
        median_fvptc <- median(
          expression_values[
            gene_data$label == "FVPTC"
          ],
          na.rm = TRUE
        )
        
        median_fa <- median(
          expression_values[
            gene_data$label == "FA"
          ],
          na.rm = TRUE
        )
        
        ifelse(
          median_fvptc > median_fa,
          "Higher in FVPTC",
          "Higher in FA"
        )
      }
      
      roc_predictor <- if (
        identical(
          direction_label,
          "Higher in FA"
        )
      ) {
        -expression_values
      } else {
        expression_values
      }
      
      roc_object <- pROC::roc(
        response = gene_data$label,
        predictor = roc_predictor,
        levels = c(
          "FA",
          "FVPTC"
        ),
        direction = "<",
        quiet = TRUE
      )
      
      auc_value <- as.numeric(
        pROC::auc(
          roc_object
        )
      )
      
      auc_ci <- tryCatch(
        {
          as.numeric(
            pROC::ci.auc(
              roc_object,
              method = "delong"
            )
          )
        },
        error = function(e) {
          c(
            NA_real_,
            NA_real_,
            NA_real_
          )
        }
      )
      
      fdr_result <- univariate_gene_report %>%
        filter(
          Gene == gene_name
        ) %>%
        pull(
          Wilcoxon_BH_FDR
        )
      
      gene_fdr <- if (
        length(fdr_result) > 0
      ) {
        as.numeric(fdr_result[1])
      } else {
        NA_real_
      }
      
      data.frame(
        Gene = gene_name,
        AUC = auc_value,
        CI_lower = auc_ci[1],
        CI_upper = auc_ci[3],
        Expression_direction = direction_label,
        Wilcoxon_BH_FDR = gene_fdr,
        Primary_panel = gene_name %in% primary_panel,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  arrange(
    AUC
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = Gene
    ),
    FDR_label = format_fdr(
      Wilcoxon_BH_FDR
    ),
    Direction_short = ifelse(
      Expression_direction == "Higher in FVPTC",
      "↑ FVPTC",
      "↑ FA"
    ),
    Display_label = paste0(
      sprintf("%.2f", AUC),
      "; ",
      Direction_short,
      "; ",
      FDR_label
    )
  )

single_gene_roc_summary$Label_x <- 1.02

p_single_gene_auc <- ggplot(
  single_gene_roc_summary,
  aes(
    x = AUC,
    y = Gene
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "gray55"
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_upper,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.75,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_lower,
      y = as.numeric(Gene) - 0.08,
      yend = as.numeric(Gene) + 0.08
    ),
    linewidth = 0.70,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = CI_upper,
      xend = CI_upper,
      y = as.numeric(Gene) - 0.08,
      yend = as.numeric(Gene) + 0.08
    ),
    linewidth = 0.70,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      fill = Primary_panel
    ),
    shape = 21,
    size = 4.2,
    color = "black",
    stroke = 0.70
  ) +
  geom_text(
    aes(
      x = Label_x,
      label = Display_label
    ),
    hjust = 0,
    family = font_family,
    size = 2.75
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#1B9E77",
      "FALSE" = "#BDBDBD"
    )
  ) +
  scale_x_continuous(
    limits = c(
      0.45,
      1.27
    ),
    breaks = seq(
      0.5,
      1,
      by = 0.1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "C. Descriptive single-gene discrimination",
    subtitle = paste0(
      "Direction-aligned ROC-AUC with DeLong 95% CI; ",
      "green points denote primary-panel genes"
    ),
    x = "Single-gene ROC-AUC",
    y = NULL
  ) +
  guides(
    fill = "none"
  ) +
  figure5_theme +
  theme(
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      color = "gray93",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 35,
      b = 8,
      l = 8
    )
  )

# ==============================================================================
# 44. Assemble corrected Figure 5
# ==============================================================================

expression_bottom_row <- (
  p_primary_expression |
    p_single_gene_auc
) +
  plot_layout(
    widths = c(
      1.70,
      1.00
    ),
    guides = "keep"
  )

figure_5 <- (
  p_expression_heatmap /
    expression_bottom_row
) +
  plot_layout(
    heights = c(
      1.10,
      1.00
    ),
    guides = "keep"
  )

print(
  figure_5
)

# ==============================================================================
# 45. Save corrected Figure 5 and exact plotted data
# ==============================================================================

save_publication_plot(
  plot = figure_5,
  filename_base =
    "Figure_10_Expression_Landscape",
  width = 18,
  height = 11.5,
  dpi = 600
)

write.csv(
  sample_level_order,
  file.path(
    figure_dir,
    "Figure_5A_Corrected_Sample_Order.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    expression_matrix_scaled
  ) %>%
    tibble::rownames_to_column(
      "Gene"
    ),
  file.path(
    figure_dir,
    "Figure_5A_Corrected_Standardized_Expression_Matrix.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_expression_long,
  file.path(
    figure_dir,
    "Figure_5B_Corrected_Primary_Expression_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_fdr_labels,
  file.path(
    figure_dir,
    "Figure_5B_Corrected_FDR_Labels.csv"
  ),
  row.names = FALSE
)

write.csv(
  single_gene_roc_summary,
  file.path(
    figure_dir,
    "Figure_5C_Corrected_Single_Gene_ROC.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Figure 5 generated successfully.\n",
  "Output directory:\n",
  figure_dir,
  "\n\nCorrected primary panel:\n",
  paste(
    primary_panel,
    collapse = " | "
  ),
  "\n"
)













# ==============================================================================
# SUPPLEMENTARY FIGURE S1
# CORRECTED SECONDARY AND WEAK/MODEL-LIMITED CANDIDATE EXPRESSION
# ==============================================================================

cat("\n--- Generating Supplementary Figure S1: Corrected secondary candidate expression ---\n")

# ------------------------------------------------------------------------------
# S1.1 Packages
# ------------------------------------------------------------------------------

required_s1_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales",
  "pROC"
)

missing_s1_packages <- required_s1_packages[
  !vapply(
    required_s1_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_s1_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_s1_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(pROC)

# ------------------------------------------------------------------------------
# S1.2 Helper functions
# ------------------------------------------------------------------------------

if (!exists("read_required_csv")) {
  
  read_required_csv <- function(
    path,
    required_columns = NULL
  ) {
    
    if (!file.exists(path)) {
      stop(
        "Required file was not found:\n",
        path
      )
    }
    
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!is.null(required_columns)) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(output)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Missing columns in:\n",
          path,
          "\nMissing: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }
    
    output
  }
}

if (!exists("save_publication_plot")) {
  
  save_publication_plot <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
  ) {
    
    output_base <- file.path(
      figure_dir,
      filename_base
    )
    
    ggsave(
      filename = paste0(output_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".tiff"),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".pdf"),
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

format_fdr_s1 <- function(x) {
  
  ifelse(
    is.na(x),
    "FDR = NA",
    ifelse(
      x < 0.001,
      "FDR < 0.001",
      paste0(
        "FDR = ",
        sprintf("%.3f", x)
      )
    )
  )
}

extract_corrected_genes_s1 <- function(
    selected_panels_table,
    panel_name
) {
  
  selected_row <- selected_panels_table %>%
    filter(
      Panel == panel_name
    )
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

# ------------------------------------------------------------------------------
# S1.3 Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

expression_file <- file.path(
  project_dir,
  "final_data_clean.rds"
)

sample_label_file <- file.path(
  stage2_old_dir,
  "Report_Sample_Label_Order_Used_By_ML.csv"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

corrected_univariate_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Univariate_Expression_Wilcoxon.csv"
)

corrected_selected_panels_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

s1_required_files <- c(
  expression_file,
  sample_label_file,
  corrected_consensus_file,
  corrected_univariate_file,
  corrected_selected_panels_file
)

s1_missing_files <- s1_required_files[
  !file.exists(s1_required_files)
]

if (length(s1_missing_files) > 0) {
  stop(
    "These Supplementary Figure S1 source files are missing:\n",
    paste(s1_missing_files, collapse = "\n")
  )
}

final_data <- readRDS(
  expression_file
)

final_data <- as.matrix(
  final_data
)

sample_labels <- read_required_csv(
  sample_label_file,
  required_columns = c(
    "Sample_order",
    "Sample_name",
    "Label"
  )
) %>%
  arrange(
    Sample_order
  )

consensus_table <- read_required_csv(
  corrected_consensus_file,
  required_columns = c(
    "Consensus_rank",
    "Gene",
    "Normalized_consensus_score",
    "Positive_models",
    "Mean_positive_frequency_across_models",
    "Expression_direction",
    "Wilcoxon_BH_FDR"
  )
)

univariate_gene_report <- read_required_csv(
  corrected_univariate_file,
  required_columns = c(
    "Gene",
    "Median_FA",
    "Median_FVPTC",
    "Expression_direction",
    "Wilcoxon_BH_FDR"
  )
)

selected_panels_corrected <- read_required_csv(
  corrected_selected_panels_file,
  required_columns = c(
    "Panel",
    "Genes"
  )
)

# ------------------------------------------------------------------------------
# S1.4 Define corrected non-primary genes
# ------------------------------------------------------------------------------

primary_panel <- extract_corrected_genes_s1(
  selected_panels_corrected,
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- extract_corrected_genes_s1(
  selected_panels_corrected,
  "Corrected_secondary_candidates"
)

expression_only_genes <- extract_corrected_genes_s1(
  selected_panels_corrected,
  "Expression_confirmed_not_ML_supported"
)

non_primary_genes <- unique(
  c(
    secondary_candidates,
    expression_only_genes
  )
)

non_primary_genes <- non_primary_genes[
  nzchar(non_primary_genes)
]

if (length(non_primary_genes) == 0) {
  non_primary_genes <- setdiff(
    consensus_table$Gene,
    primary_panel
  )
}

non_primary_genes <- consensus_table %>%
  filter(
    Gene %in% non_primary_genes
  ) %>%
  arrange(
    Consensus_rank
  ) %>%
  pull(
    Gene
  ) %>%
  as.character()

if (length(non_primary_genes) == 0) {
  stop(
    "No corrected secondary or expression-only genes were found for S1."
  )
}

missing_non_primary_genes <- setdiff(
  non_primary_genes,
  rownames(final_data)
)

if (length(missing_non_primary_genes) > 0) {
  stop(
    "These non-primary genes are absent from final_data_clean.rds:\n",
    paste(missing_non_primary_genes, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# S1.5 Sample labels
# ------------------------------------------------------------------------------

if (!identical(
  as.character(sample_labels$Sample_name),
  colnames(final_data)
)) {
  stop(
    "Sample order mismatch between final_data_clean.rds and ",
    "Report_Sample_Label_Order_Used_By_ML.csv."
  )
}

sample_labels <- sample_labels %>%
  mutate(
    Label = factor(
      Label,
      levels = c(
        "FA",
        "FVPTC"
      )
    )
  )

if (any(is.na(sample_labels$Label))) {
  stop(
    "At least one sample label is not FA or FVPTC."
  )
}

gene_data <- as.data.frame(
  t(
    final_data[
      non_primary_genes,
      ,
      drop = FALSE
    ]
  ),
  stringsAsFactors = FALSE
)

gene_data$Sample <- rownames(
  gene_data
)

gene_data <- gene_data %>%
  left_join(
    sample_labels %>%
      transmute(
        Sample = as.character(
          Sample_name
        ),
        label = Label
      ),
    by = "Sample"
  )

if (any(is.na(gene_data$label))) {
  stop(
    "At least one expression sample could not be matched to its label."
  )
}

# ------------------------------------------------------------------------------
# S1.6 Theme and colors
# ------------------------------------------------------------------------------

if (!exists("font_family")) {
  font_family <- "Arial"
}

if (!exists("group_colors")) {
  group_colors <- c(
    "FA" = "#3B78A8",
    "FVPTC" = "#D9534F"
  )
}
s1_group_colors <- c(
  "Secondary candidate" = "#E6AB02",
  "Weak/model-limited secondary" = "#D95F02"
)

s1_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "gray30",
      margin = margin(
        b = 8
      )
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# ==============================================================================
# PANEL A
# EXPRESSION DISTRIBUTIONS OF NON-PRIMARY CANDIDATES
# ==============================================================================

secondary_expression_long <- gene_data %>%
  select(
    Sample,
    label,
    all_of(non_primary_genes)
  ) %>%
  pivot_longer(
    cols = all_of(non_primary_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = non_primary_genes
    )
  )

secondary_fdr_labels <- secondary_expression_long %>%
  group_by(
    Gene
  ) %>%
  summarise(
    expression_min = min(
      Expression,
      na.rm = TRUE
    ),
    expression_max = max(
      Expression,
      na.rm = TRUE
    ),
    expression_range =
      expression_max -
      expression_min,
    y_position =
      expression_max +
      ifelse(
        expression_range > 0,
        0.12 * expression_range,
        0.25
      ),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = as.character(Gene)
  ) %>%
  left_join(
    univariate_gene_report %>%
      filter(
        Gene %in% non_primary_genes
      ) %>%
      transmute(
        Gene = as.character(Gene),
        Wilcoxon_BH_FDR = as.numeric(
          Wilcoxon_BH_FDR
        )
      ),
    by = "Gene"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = non_primary_genes
    ),
    FDR_label = format_fdr_s1(
      Wilcoxon_BH_FDR
    )
  )

p_s1_expression <- ggplot(
  secondary_expression_long,
  aes(
    x = label,
    y = Expression,
    fill = label,
    color = label
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.28,
    width = 0.88,
    linewidth = 0.55
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.92,
    linewidth = 0.50
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.45,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  geom_text(
    data = secondary_fdr_labels,
    aes(
      x = 1.5,
      y = y_position,
      label = FDR_label
    ),
    inherit.aes = FALSE,
    family = font_family,
    size = 3.05
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    title = "A. Expression distributions of non-primary candidates",
    subtitle = paste0(
      "Secondary and weak/model-limited candidates after ",
      "normalized Borda prioritization"
    ),
    x = NULL,
    y = "Expression level (log2)"
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  s1_theme +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.major.y = element_line(
      color = "gray93",
      linewidth = 0.35
    )
  )

# ==============================================================================
# PANEL B
# SINGLE-GENE ROC-AUC OF NON-PRIMARY CANDIDATES
# ==============================================================================

s1_single_gene_roc <- bind_rows(
  lapply(
    non_primary_genes,
    function(gene_name) {
      
      expression_values <- as.numeric(
        gene_data[[gene_name]]
      )
      
      direction_result <- univariate_gene_report %>%
        filter(
          Gene == gene_name
        ) %>%
        pull(
          Expression_direction
        )
      
      direction_label <- if (
        length(direction_result) > 0 &&
        !is.na(direction_result[1])
      ) {
        as.character(direction_result[1])
      } else {
        
        median_fvptc <- median(
          expression_values[
            gene_data$label == "FVPTC"
          ],
          na.rm = TRUE
        )
        
        median_fa <- median(
          expression_values[
            gene_data$label == "FA"
          ],
          na.rm = TRUE
        )
        
        ifelse(
          median_fvptc > median_fa,
          "Higher in FVPTC",
          "Higher in FA"
        )
      }
      
      roc_predictor <- if (
        identical(
          direction_label,
          "Higher in FA"
        )
      ) {
        -expression_values
      } else {
        expression_values
      }
      
      roc_object <- pROC::roc(
        response = gene_data$label,
        predictor = roc_predictor,
        levels = c(
          "FA",
          "FVPTC"
        ),
        direction = "<",
        quiet = TRUE
      )
      
      auc_value <- as.numeric(
        pROC::auc(
          roc_object
        )
      )
      
      auc_ci <- tryCatch(
        {
          as.numeric(
            pROC::ci.auc(
              roc_object,
              method = "delong"
            )
          )
        },
        error = function(e) {
          c(
            NA_real_,
            NA_real_,
            NA_real_
          )
        }
      )
      
      fdr_result <- univariate_gene_report %>%
        filter(
          Gene == gene_name
        ) %>%
        pull(
          Wilcoxon_BH_FDR
        )
      
      gene_fdr <- if (
        length(fdr_result) > 0
      ) {
        as.numeric(fdr_result[1])
      } else {
        NA_real_
      }
      
      group_label <- case_when(
        gene_name == "SPP1" ~
          "Weak/model-limited secondary",
        gene_name %in% secondary_candidates ~
          "Secondary candidate",
        TRUE ~
          "Expression-only"
      )
      
      data.frame(
        Gene = gene_name,
        AUC = auc_value,
        CI_lower = auc_ci[1],
        CI_upper = auc_ci[3],
        Expression_direction = direction_label,
        Wilcoxon_BH_FDR = gene_fdr,
        Evidence_group = group_label,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  arrange(
    AUC
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = Gene
    ),
    Evidence_group = factor(
      Evidence_group,
      levels = names(s1_group_colors)
    ),
    Direction_short = ifelse(
      Expression_direction == "Higher in FVPTC",
      "↑ FVPTC",
      "↑ FA"
    ),
    Display_label = paste0(
      sprintf("%.2f", AUC),
      "; ",
      Direction_short,
      "; ",
      format_fdr_s1(Wilcoxon_BH_FDR)
    ),
    Label_x = 1.02
  )

p_s1_auc <- ggplot(
  s1_single_gene_roc,
  aes(
    x = AUC,
    y = Gene
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "gray55"
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_upper,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.80,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_lower,
      y = as.numeric(Gene) - 0.08,
      yend = as.numeric(Gene) + 0.08
    ),
    linewidth = 0.70,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = CI_upper,
      xend = CI_upper,
      y = as.numeric(Gene) - 0.08,
      yend = as.numeric(Gene) + 0.08
    ),
    linewidth = 0.70,
    color = "#333333",
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      fill = Evidence_group
    ),
    shape = 21,
    size = 4.3,
    color = "black",
    stroke = 0.70
  ) +
  geom_text(
    aes(
      x = Label_x,
      label = Display_label
    ),
    hjust = 0,
    family = font_family,
    size = 3.05
  ) +
  scale_fill_manual(
    values = s1_group_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      0.45,
      1.25
    ),
    breaks = seq(
      0.5,
      1,
      by = 0.1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "B. Descriptive single-gene discrimination",
    subtitle = "Direction-aligned ROC-AUC with DeLong 95% CI",
    x = "Single-gene ROC-AUC",
    y = NULL,
    fill = NULL
  ) +
  s1_theme +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      color = "gray93",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 35,
      b = 8,
      l = 8
    )
  )

# ==============================================================================
# PANEL C
# CORRECTED ML AND EXPRESSION EVIDENCE SUMMARY
# ==============================================================================

s1_evidence_summary <- consensus_table %>%
  filter(
    Gene %in% non_primary_genes
  ) %>%
  arrange(
    Consensus_rank
  ) %>%
  mutate(
    Evidence_group = case_when(
      Gene == "SPP1" ~
        "Weak/model-limited secondary",
      Gene %in% secondary_candidates ~
        "Secondary candidate",
      TRUE ~
        "Expression-only"
    ),
    Evidence_group = factor(
      Evidence_group,
      levels = names(s1_group_colors)
    ),
    Gene = factor(
      Gene,
      levels = rev(non_primary_genes)
    ),
    Consensus_label = sprintf(
      "%.3f",
      Normalized_consensus_score
    ),
    Models_label = paste0(
      Positive_models,
      "/3"
    ),
    Frequency_label = sprintf(
      "%.2f",
      Mean_positive_frequency_across_models
    ),
    FDR_label = format_fdr_s1(
      Wilcoxon_BH_FDR
    )
  )

s1_evidence_long <- bind_rows(
  s1_evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "Normalized\nconsensus",
      Label = Consensus_label
    ),
  s1_evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "ML-positive\nmodels",
      Label = Models_label
    ),
  s1_evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "Mean positive\nfrequency",
      Label = Frequency_label
    ),
  s1_evidence_summary %>%
    transmute(
      Gene,
      Evidence_group,
      Metric = "BH-FDR",
      Label = FDR_label
    )
) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "Normalized\nconsensus",
        "ML-positive\nmodels",
        "Mean positive\nfrequency",
        "BH-FDR"
      )
    )
  )

p_s1_evidence <- ggplot(
  s1_evidence_long,
  aes(
    x = Metric,
    y = Gene,
    fill = Evidence_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Label
    ),
    family = font_family,
    size = 3.15,
    color = "black",
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = s1_group_colors,
    drop = FALSE
  ) +
  labs(
    title = "C. Non-primary evidence summary",
    subtitle = "SPP1 is retained as weak/model-limited rather than primary-panel supported",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  s1_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      face = "bold",
      size = 9.2
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

# ------------------------------------------------------------------------------
# S1.7 Assemble and save Supplementary Figure S1
# ------------------------------------------------------------------------------

supplementary_figure_s1 <- (
  p_s1_expression /
    (
      p_s1_auc |
        p_s1_evidence
    )
) +
  plot_layout(
    heights = c(
      1.05,
      1.00
    ),
    widths = c(
      1.00,
      1.05
    )
  )

print(
  supplementary_figure_s1
)

save_publication_plot(
  plot = supplementary_figure_s1,
  filename_base =
    "Figure_S1_Secondary_Candidates",
  width = 16,
  height = 9.5,
  dpi = 600
)

write.csv(
  secondary_expression_long,
  file.path(
    figure_dir,
    "Supplementary_Figure_S1A_Corrected_Secondary_Expression_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  secondary_fdr_labels,
  file.path(
    figure_dir,
    "Supplementary_Figure_S1A_Corrected_FDR_Labels.csv"
  ),
  row.names = FALSE
)

write.csv(
  s1_single_gene_roc,
  file.path(
    figure_dir,
    "Supplementary_Figure_S1B_Corrected_Secondary_ROC_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  s1_evidence_long,
  file.path(
    figure_dir,
    "Supplementary_Figure_S1C_Corrected_Evidence_Data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Supplementary Figure S1 generated successfully.\n",
  "Genes shown:\n",
  paste(
    non_primary_genes,
    collapse = ", "
  ),
  "\nOutput directory:\n",
  figure_dir,
  "\n"
)





















# ==============================================================================
# SUPPLEMENTARY FIGURE S2
# CORRECTED ALL-EIGHT-GENE SPEARMAN CORRELATION
# ==============================================================================

cat("\n--- Generating Supplementary Figure S2: All-8-gene Spearman correlation ---\n")

# ------------------------------------------------------------------------------
# S2.1 Packages
# ------------------------------------------------------------------------------

required_s2_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales",
  "tibble"
)

missing_s2_packages <- required_s2_packages[
  !vapply(
    required_s2_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_s2_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_s2_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(tibble)

# ------------------------------------------------------------------------------
# S2.2 Helper functions
# ------------------------------------------------------------------------------

if (!exists("read_required_csv")) {
  
  read_required_csv <- function(
    path,
    required_columns = NULL
  ) {
    
    if (!file.exists(path)) {
      stop(
        "Required file was not found:\n",
        path
      )
    }
    
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!is.null(required_columns)) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(output)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Missing columns in:\n",
          path,
          "\nMissing: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }
    
    output
  }
}

if (!exists("save_publication_plot")) {
  
  save_publication_plot <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
  ) {
    
    output_base <- file.path(
      figure_dir,
      filename_base
    )
    
    ggsave(
      filename = paste0(output_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".tiff"),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".pdf"),
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

# ------------------------------------------------------------------------------
# S2.3 Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

expression_file <- file.path(
  project_dir,
  "final_data_clean.rds"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

s2_required_files <- c(
  expression_file,
  corrected_consensus_file
)

s2_missing_files <- s2_required_files[
  !file.exists(s2_required_files)
]

if (length(s2_missing_files) > 0) {
  stop(
    "These Supplementary Figure S2 source files are missing:\n",
    paste(s2_missing_files, collapse = "\n")
  )
}

final_data <- readRDS(
  expression_file
)

final_data <- as.matrix(
  final_data
)

consensus_table <- read_required_csv(
  corrected_consensus_file,
  required_columns = c(
    "Consensus_rank",
    "Gene",
    "Normalized_consensus_score"
  )
)

# ------------------------------------------------------------------------------
# S2.4 Gene order and correlation calculation
# ------------------------------------------------------------------------------

candidate_genes <- consensus_table %>%
  arrange(
    Consensus_rank
  ) %>%
  pull(
    Gene
  ) %>%
  as.character()

missing_candidate_genes <- setdiff(
  candidate_genes,
  rownames(final_data)
)

if (length(missing_candidate_genes) > 0) {
  stop(
    "These candidate genes are absent from final_data_clean.rds:\n",
    paste(missing_candidate_genes, collapse = ", ")
  )
}

candidate_expression <- t(
  final_data[
    candidate_genes,
    ,
    drop = FALSE
  ]
)

spearman_matrix <- cor(
  candidate_expression,
  method = "spearman",
  use = "pairwise.complete.obs"
)

spearman_matrix <- spearman_matrix[
  candidate_genes,
  candidate_genes,
  drop = FALSE
]

spearman_long <- as.data.frame(
  as.table(
    spearman_matrix
  ),
  stringsAsFactors = FALSE
) %>%
  rename(
    Gene_1 = Var1,
    Gene_2 = Var2,
    Spearman_rho = Freq
  ) %>%
  mutate(
    Gene_1 = factor(
      Gene_1,
      levels = candidate_genes
    ),
    Gene_2 = factor(
      Gene_2,
      levels = rev(candidate_genes)
    ),
    Rho_label = sprintf(
      "%.2f",
      Spearman_rho
    )
  )

pairwise_correlation <- bind_rows(
  lapply(
    seq_along(candidate_genes),
    function(i) {
      
      if (i == length(candidate_genes)) {
        return(NULL)
      }
      
      bind_rows(
        lapply(
          (i + 1):length(candidate_genes),
          function(j) {
            
            gene_1 <- candidate_genes[i]
            gene_2 <- candidate_genes[j]
            rho_value <- spearman_matrix[
              gene_1,
              gene_2
            ]
            
            data.frame(
              Gene_1 = gene_1,
              Gene_2 = gene_2,
              Pair = paste(
                gene_1,
                gene_2,
                sep = " - "
              ),
              Spearman_rho = as.numeric(
                rho_value
              ),
              Abs_Spearman_rho = abs(
                as.numeric(
                  rho_value
                )
              ),
              High_correlation_abs_ge_0.80 =
                abs(
                  as.numeric(
                    rho_value
                  )
                ) >= 0.80,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
) %>%
  arrange(
    desc(
      Abs_Spearman_rho
    )
  ) %>%
  mutate(
    Pair = factor(
      Pair,
      levels = rev(
        Pair
      )
    ),
    Correlation_group = ifelse(
      High_correlation_abs_ge_0.80,
      "|rho| >= 0.80",
      "|rho| < 0.80"
    ),
    Correlation_group = factor(
      Correlation_group,
      levels = c(
        "|rho| >= 0.80",
        "|rho| < 0.80"
      )
    ),
    Label = sprintf(
      "%.3f",
      Spearman_rho
    )
  )

high_correlation_count <- sum(
  pairwise_correlation$High_correlation_abs_ge_0.80,
  na.rm = TRUE
)

maximum_pair <- pairwise_correlation[
  1,
  ,
  drop = FALSE
]

maximum_pair_label <- paste0(
  as.character(maximum_pair$Gene_1),
  " - ",
  as.character(maximum_pair$Gene_2),
  ": rho = ",
  sprintf(
    "%.3f",
    maximum_pair$Spearman_rho
  )
)

# ------------------------------------------------------------------------------
# S2.5 Theme and colors
# ------------------------------------------------------------------------------

if (!exists("font_family")) {
  font_family <- "Arial"
}

s2_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "gray30",
      margin = margin(
        b = 8
      )
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

correlation_colors <- c(
  "|rho| >= 0.80" = "#B2182B",
  "|rho| < 0.80" = "#4D4D4D"
)

# ==============================================================================
# PANEL A
# SPEARMAN CORRELATION MATRIX
# ==============================================================================

p_s2_matrix <- ggplot(
  spearman_long,
  aes(
    x = Gene_1,
    y = Gene_2,
    fill = Spearman_rho
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Rho_label
    ),
    family = font_family,
    size = 3.1,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B78A8",
    mid = "white",
    high = "#D9534F",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    breaks = c(
      -1,
      -0.5,
      0,
      0.5,
      1
    ),
    name = "Spearman\nrho"
  ) +
  coord_fixed() +
  labs(
    title = "A. Spearman correlation matrix of all eight candidate genes",
    subtitle = "Genes are ordered by normalized consensus rank",
    x = NULL,
    y = NULL
  ) +
  s2_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ==============================================================================
# PANEL B
# PAIRWISE ABSOLUTE CORRELATION RANKING
# ==============================================================================

p_s2_pairs <- ggplot(
  pairwise_correlation,
  aes(
    x = Abs_Spearman_rho,
    y = Pair,
    color = Correlation_group
  )
) +
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed",
    color = "#B2182B",
    linewidth = 0.75
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = Abs_Spearman_rho,
      y = Pair,
      yend = Pair
    ),
    linewidth = 0.75,
    color = "gray70"
  ) +
  geom_point(
    size = 3.4
  ) +
  geom_text(
    aes(
      label = Label
    ),
    hjust = -0.18,
    family = font_family,
    size = 2.9,
    color = "black"
  ) +
  scale_color_manual(
    values = correlation_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(
        0.86,
        max(
          pairwise_correlation$Abs_Spearman_rho,
          na.rm = TRUE
        ) + 0.08
      )
    ),
    breaks = seq(
      0,
      0.8,
      by = 0.2
    ),
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "B. Pairwise absolute Spearman correlation",
    subtitle = paste0(
      "High-correlation threshold: |rho| >= 0.80; ",
      "pairs detected = ",
      high_correlation_count,
      "; maximum pair: ",
      maximum_pair_label
    ),
    x = "|Spearman rho|",
    y = NULL,
    color = NULL
  ) +
  s2_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      size = 8.6
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 35,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------------------------
# S2.6 Assemble and save Supplementary Figure S2
# ------------------------------------------------------------------------------

supplementary_figure_s2 <- (
  p_s2_matrix |
    p_s2_pairs
) +
  plot_layout(
    widths = c(
      1.00,
      1.35
    )
  )

print(
  supplementary_figure_s2
)

save_publication_plot(
  plot = supplementary_figure_s2,
  filename_base =
    "Figure_S2_Spearman_Correlation",
  width = 17,
  height = 8.5,
  dpi = 600
)

write.csv(
  spearman_long,
  file.path(
    figure_dir,
    "Supplementary_Figure_S2A_Corrected_All8_Spearman_Matrix_Long.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    spearman_matrix
  ) %>%
    tibble::rownames_to_column(
      "Gene"
    ),
  file.path(
    figure_dir,
    "Supplementary_Figure_S2A_Corrected_All8_Spearman_Matrix_Wide.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_correlation,
  file.path(
    figure_dir,
    "Supplementary_Figure_S2B_Corrected_All8_Spearman_Pairs.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Supplementary Figure S2 generated successfully.\n",
  "High-correlation pairs with |rho| >= 0.80: ",
  high_correlation_count,
  "\nMaximum pair: ",
  maximum_pair_label,
  "\nOutput directory:\n",
  figure_dir,
  "\n"
)




















# ==============================================================================
# SUPPLEMENTARY FIGURE S3
# CORRECTED CROSS-FITTED PERMUTATION-IMPORTANCE RANK STABILITY
# ==============================================================================

cat("\n--- Generating Supplementary Figure S3: Corrected PI rank stability ---\n")

# ------------------------------------------------------------------------------
# S3.1 Packages
# ------------------------------------------------------------------------------

required_s3_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales"
)

missing_s3_packages <- required_s3_packages[
  !vapply(
    required_s3_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_s3_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_s3_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ------------------------------------------------------------------------------
# S3.2 Helper functions
# ------------------------------------------------------------------------------

if (!exists("read_required_csv")) {
  
  read_required_csv <- function(
    path,
    required_columns = NULL
  ) {
    
    if (!file.exists(path)) {
      stop(
        "Required file was not found:\n",
        path
      )
    }
    
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!is.null(required_columns)) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(output)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Missing columns in:\n",
          path,
          "\nMissing: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }
    
    output
  }
}

if (!exists("save_publication_plot")) {
  
  save_publication_plot <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
  ) {
    
    output_base <- file.path(
      figure_dir,
      filename_base
    )
    
    ggsave(
      filename = paste0(output_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".tiff"),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".pdf"),
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

extract_corrected_genes_s3 <- function(
    selected_panels_table,
    panel_name
) {
  
  selected_row <- selected_panels_table %>%
    filter(
      Panel == panel_name
    )
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

# ------------------------------------------------------------------------------
# S3.3 Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

rank_stability_file <- file.path(
  stage2_corrected_dir,
  "Table_CrossFitted_PI_Rank_Stability.csv"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

corrected_selected_panels_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

s3_required_files <- c(
  rank_stability_file,
  corrected_consensus_file,
  corrected_selected_panels_file
)

s3_missing_files <- s3_required_files[
  !file.exists(s3_required_files)
]

if (length(s3_missing_files) > 0) {
  stop(
    "These Supplementary Figure S3 source files are missing:\n",
    paste(s3_missing_files, collapse = "\n")
  )
}

rank_stability <- read_required_csv(
  rank_stability_file,
  required_columns = c(
    "Model",
    "Gene",
    "Positive_resample_frequency",
    "Top4_resample_frequency",
    "Median_positive_rank",
    "Rank_IQR",
    "Number_outer_resamples"
  )
)

consensus_table <- read_required_csv(
  corrected_consensus_file,
  required_columns = c(
    "Consensus_rank",
    "Gene",
    "Normalized_consensus_score"
  )
)

selected_panels_corrected <- read_required_csv(
  corrected_selected_panels_file,
  required_columns = c(
    "Panel",
    "Genes"
  )
)

# ------------------------------------------------------------------------------
# S3.4 Gene order, model order, and evidence groups
# ------------------------------------------------------------------------------

if (!exists("final_model_names")) {
  
  selected_models_file <- file.path(
    stage2_old_dir,
    "Table_Final_Selected_Models_From_Stage1.csv"
  )
  
  selected_models_table <- read_required_csv(
    selected_models_file,
    required_columns = c(
      "Final_model_rank",
      "Model"
    )
  )
  
  final_model_names <- selected_models_table %>%
    arrange(
      Final_model_rank
    ) %>%
    pull(
      Model
    ) %>%
    as.character()
}

top3_models <- final_model_names

candidate_genes <- consensus_table %>%
  arrange(
    Consensus_rank
  ) %>%
  pull(
    Gene
  ) %>%
  as.character()

primary_panel <- extract_corrected_genes_s3(
  selected_panels_corrected,
  "Corrected_primary_dynamic_elbow_panel"
)

secondary_candidates <- extract_corrected_genes_s3(
  selected_panels_corrected,
  "Corrected_secondary_candidates"
)

expression_only_genes <- extract_corrected_genes_s3(
  selected_panels_corrected,
  "Expression_confirmed_not_ML_supported"
)

if (!exists("font_family")) {
  font_family <- "Arial"
}

if (!exists("default_model_colors")) {
  default_model_colors <- c(
    "XGBoost" = "#0072B2",
    "Naive Bayes" = "#D55E00",
    "KNN" = "#009E73"
  )
}

model_colors <- default_model_colors[
  top3_models
]

if (any(is.na(model_colors))) {
  
  fallback_colors <- c(
    "#0072B2",
    "#D55E00",
    "#009E73"
  )
  
  model_colors <- setNames(
    fallback_colors[
      seq_along(top3_models)
    ],
    top3_models
  )
}

panel_colors <- c(
  "Primary corrected panel" = "#1B9E77",
  "Secondary candidate" = "#E6AB02",
  "Weak/model-limited secondary" = "#D95F02",
  "Expression-only" = "#BDBDBD"
)

rank_stability <- rank_stability %>%
  filter(
    Model %in% top3_models,
    Gene %in% candidate_genes
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = top3_models
    ),
    Gene = factor(
      Gene,
      levels = rev(candidate_genes)
    ),
    Positive_resample_frequency = as.numeric(
      Positive_resample_frequency
    ),
    Top4_resample_frequency = as.numeric(
      Top4_resample_frequency
    ),
    Median_positive_rank = as.numeric(
      Median_positive_rank
    ),
    Rank_IQR = as.numeric(
      Rank_IQR
    ),
    Positive_label = sprintf(
      "%.2f",
      Positive_resample_frequency
    ),
    Top4_label = sprintf(
      "%.2f",
      Top4_resample_frequency
    ),
    Median_rank_label = ifelse(
      is.na(Median_positive_rank),
      "NA",
      sprintf("%.1f", Median_positive_rank)
    )
  )

outer_resample_count <- unique(
  rank_stability$Number_outer_resamples
)

outer_resample_label <- if (
  length(outer_resample_count) == 1
) {
  paste0(
    outer_resample_count,
    " outer resamples"
  )
} else {
  "outer resamples"
}

# ------------------------------------------------------------------------------
# S3.5 Theme
# ------------------------------------------------------------------------------

s3_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "gray30",
      margin = margin(
        b = 8
      )
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 9.0
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# ==============================================================================
# PANEL A
# POSITIVE OUTER-RESAMPLE FREQUENCY
# ==============================================================================

p_s3_positive <- ggplot(
  rank_stability,
  aes(
    x = Model,
    y = Gene,
    fill = Positive_resample_frequency
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Positive_label
    ),
    family = font_family,
    size = 3.15,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#1B9E77",
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    name = "Positive\nfrequency"
  ) +
  labs(
    title = "A. Positive outer-resample frequency",
    subtitle = paste0(
      "Fraction of ",
      outer_resample_label,
      " with positive mean ΔAUC"
    ),
    x = NULL,
    y = NULL
  ) +
  s3_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ==============================================================================
# PANEL B
# TOP-FOUR OUTER-RESAMPLE FREQUENCY
# ==============================================================================

p_s3_top4 <- ggplot(
  rank_stability,
  aes(
    x = Model,
    y = Gene,
    fill = Top4_resample_frequency
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Top4_label
    ),
    family = font_family,
    size = 3.15,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF7EC",
    high = "#D95F02",
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    name = "Top-four\nfrequency"
  ) +
  labs(
    title = "B. Top-four rank frequency",
    subtitle = paste0(
      "Fraction of ",
      outer_resample_label,
      " in which the gene ranked among top four"
    ),
    x = NULL,
    y = NULL
  ) +
  s3_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ==============================================================================
# PANEL C
# GENE-LEVEL STABILITY SUMMARY
# ==============================================================================

gene_stability_summary <- rank_stability %>%
  mutate(
    Gene_character = as.character(Gene)
  ) %>%
  group_by(
    Gene_character
  ) %>%
  summarise(
    Mean_positive_frequency = mean(
      Positive_resample_frequency,
      na.rm = TRUE
    ),
    Mean_top4_frequency = mean(
      Top4_resample_frequency,
      na.rm = TRUE
    ),
    Median_positive_rank_across_models = median(
      Median_positive_rank,
      na.rm = TRUE
    ),
    Mean_rank_IQR = mean(
      Rank_IQR,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  rename(
    Gene = Gene_character
  ) %>%
  left_join(
    consensus_table %>%
      transmute(
        Gene = as.character(Gene),
        Consensus_rank,
        Normalized_consensus_score
      ),
    by = "Gene"
  ) %>%
  mutate(
    Evidence_group = case_when(
      Gene %in% primary_panel ~
        "Primary corrected panel",
      Gene == "SPP1" ~
        "Weak/model-limited secondary",
      Gene %in% secondary_candidates ~
        "Secondary candidate",
      Gene %in% expression_only_genes ~
        "Expression-only",
      TRUE ~
        "Secondary candidate"
    ),
    Evidence_group = factor(
      Evidence_group,
      levels = names(panel_colors)
    ),
    Gene = factor(
      Gene,
      levels = rev(candidate_genes)
    )
  ) %>%
  arrange(
    Consensus_rank
  )

stability_long <- gene_stability_summary %>%
  select(
    Gene,
    Evidence_group,
    Mean_positive_frequency,
    Mean_top4_frequency
  ) %>%
  pivot_longer(
    cols = c(
      Mean_positive_frequency,
      Mean_top4_frequency
    ),
    names_to = "Metric",
    values_to = "Frequency"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Mean_positive_frequency" = "Mean positive frequency",
      "Mean_top4_frequency" = "Mean top-four frequency"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Mean positive frequency",
        "Mean top-four frequency"
      )
    ),
    Label = sprintf(
      "%.2f",
      Frequency
    )
  )

p_s3_gene_summary <- ggplot(
  stability_long,
  aes(
    x = Frequency,
    y = Gene,
    color = Metric,
    shape = Metric
  )
) +
  geom_point(
    size = 3.6,
    position = position_dodge(
      width = 0.55
    )
  ) +
  geom_text(
    aes(
      label = Label
    ),
    hjust = -0.15,
    size = 2.9,
    family = font_family,
    position = position_dodge(
      width = 0.55
    ),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Mean positive frequency" = "#1B9E77",
      "Mean top-four frequency" = "#D95F02"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Mean positive frequency" = 16,
      "Mean top-four frequency" = 17
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1.10
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    title = "C. Gene-level stability summary",
    subtitle = "Frequencies averaged across the three final classifiers",
    x = "Frequency",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  s3_theme +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# PANEL D
# MEDIAN POSITIVE RANK BY MODEL
# ==============================================================================

p_s3_median_rank <- ggplot(
  rank_stability,
  aes(
    x = Model,
    y = Gene,
    fill = Median_positive_rank
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = Median_rank_label
    ),
    family = font_family,
    size = 3.15,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#1B9E77",
    high = "#FEE08B",
    na.value = "gray90",
    name = "Median\npositive\nrank"
  ) +
  labs(
    title = "D. Median positive rank",
    subtitle = "Lower values indicate more stable high-ranking importance",
    x = NULL,
    y = NULL
  ) +
  s3_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# S3.6 Assemble and save Supplementary Figure S3
# ------------------------------------------------------------------------------

supplementary_figure_s3 <- (
  p_s3_positive |
    p_s3_top4
) / (
  p_s3_gene_summary |
    p_s3_median_rank
) +
  plot_layout(
    heights = c(
      1.00,
      1.05
    ),
    widths = c(
      1.00,
      1.00
    )
  )

print(
  supplementary_figure_s3
)

save_publication_plot(
  plot = supplementary_figure_s3,
  filename_base =
    "Figure_S3_Rank_Stability",
  width = 17,
  height = 11.5,
  dpi = 600
)

write.csv(
  rank_stability,
  file.path(
    figure_dir,
    "Supplementary_Figure_S3A_B_D_Corrected_Rank_Stability_Data.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_stability_summary,
  file.path(
    figure_dir,
    "Supplementary_Figure_S3C_Corrected_Gene_Level_Stability_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  stability_long,
  file.path(
    figure_dir,
    "Supplementary_Figure_S3C_Corrected_Gene_Level_Stability_Long.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Supplementary Figure S3 generated successfully.\n",
  "Output directory:\n",
  figure_dir,
  "\n"
)

cat(
  "\nGene-level stability summary:\n"
)

print(
  gene_stability_summary %>%
    select(
      Gene,
      Consensus_rank,
      Normalized_consensus_score,
      Mean_positive_frequency,
      Mean_top4_frequency,
      Median_positive_rank_across_models,
      Mean_rank_IQR,
      Evidence_group
    ),
  row.names = FALSE,
  digits = 8
)
















# ==============================================================================
# SUPPLEMENTARY FIGURE S4
# CORRECTED PANEL-SIZE SENSITIVITY ANALYSIS
# ==============================================================================

cat("\n--- Generating Supplementary Figure S4: Corrected panel-size sensitivity ---\n")

# ------------------------------------------------------------------------------
# S4.1 Packages
# ------------------------------------------------------------------------------

required_s4_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales"
)

missing_s4_packages <- required_s4_packages[
  !vapply(
    required_s4_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_s4_packages) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_s4_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ------------------------------------------------------------------------------
# S4.2 Helper functions
# ------------------------------------------------------------------------------

if (!exists("read_required_csv")) {
  
  read_required_csv <- function(
    path,
    required_columns = NULL
  ) {
    
    if (!file.exists(path)) {
      stop(
        "Required file was not found:\n",
        path
      )
    }
    
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (!is.null(required_columns)) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(output)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Missing columns in:\n",
          path,
          "\nMissing: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }
    
    output
  }
}

if (!exists("save_publication_plot")) {
  
  save_publication_plot <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
  ) {
    
    output_base <- file.path(
      figure_dir,
      filename_base
    )
    
    ggsave(
      filename = paste0(output_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".tiff"),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
    
    ggsave(
      filename = paste0(output_base, ".pdf"),
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}

extract_corrected_genes_s4 <- function(
    selected_panels_table,
    panel_name
) {
  
  selected_row <- selected_panels_table %>%
    filter(
      Panel == panel_name
    )
  
  if (nrow(selected_row) == 0) {
    return(character(0))
  }
  
  genes_text <- selected_row$Genes[1]
  
  if (
    is.na(genes_text) ||
    !nzchar(trimws(genes_text))
  ) {
    return(character(0))
  }
  
  trimws(
    unlist(
      strsplit(
        genes_text,
        ",",
        fixed = TRUE
      )
    )
  )
}

wrap_s4_text <- function(
    x,
    width = 32
) {
  
  vapply(
    x,
    function(one_text) {
      paste(
        strwrap(
          one_text,
          width = width
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}

# ------------------------------------------------------------------------------
# S4.3 Input files
# ------------------------------------------------------------------------------

if (!exists("stage2_corrected_dir")) {
  stage2_corrected_dir <- file.path(
    project_dir,
    "ml-output",
    "stage2_methodologically_corrected"
  )
}

panel_size_file <- file.path(
  stage2_corrected_dir,
  "Table_Panel_Size_Sensitivity_Performance.csv"
)

corrected_consensus_file <- file.path(
  stage2_corrected_dir,
  "Table_Corrected_Consensus_Normalized_Borda.csv"
)

corrected_selected_panels_file <- file.path(
  stage2_corrected_dir,
  "Table_Selected_Panels_Corrected.csv"
)

s4_required_files <- c(
  panel_size_file,
  corrected_consensus_file,
  corrected_selected_panels_file
)

s4_missing_files <- s4_required_files[
  !file.exists(s4_required_files)
]

if (length(s4_missing_files) > 0) {
  stop(
    "These Supplementary Figure S4 source files are missing:\n",
    paste(s4_missing_files, collapse = "\n")
  )
}

panel_size_performance <- read_required_csv(
  panel_size_file,
  required_columns = c(
    "Number_genes",
    "Genes",
    "Model",
    "AUC",
    "Model_maximum_AUC",
    "Difference_from_model_maximum",
    "Within_0_01_of_model_maximum"
  )
)

consensus_table <- read_required_csv(
  corrected_consensus_file,
  required_columns = c(
    "Consensus_rank",
    "Gene",
    "Normalized_consensus_score"
  )
)

selected_panels_corrected <- read_required_csv(
  corrected_selected_panels_file,
  required_columns = c(
    "Panel",
    "Genes"
  )
)

# ------------------------------------------------------------------------------
# S4.4 Model order, gene order, and corrected primary-panel size
# ------------------------------------------------------------------------------

if (!exists("final_model_names")) {
  
  selected_models_file <- file.path(
    stage2_old_dir,
    "Table_Final_Selected_Models_From_Stage1.csv"
  )
  
  selected_models_table <- read_required_csv(
    selected_models_file,
    required_columns = c(
      "Final_model_rank",
      "Model"
    )
  )
  
  final_model_names <- selected_models_table %>%
    arrange(
      Final_model_rank
    ) %>%
    pull(
      Model
    ) %>%
    as.character()
}

top3_models <- final_model_names

candidate_genes <- consensus_table %>%
  arrange(
    Consensus_rank
  ) %>%
  pull(
    Gene
  ) %>%
  as.character()

primary_panel <- extract_corrected_genes_s4(
  selected_panels_corrected,
  "Corrected_primary_dynamic_elbow_panel"
)

if (length(primary_panel) == 0) {
  stop(
    "No corrected primary panel was found in Table_Selected_Panels_Corrected.csv."
  )
}

primary_panel_size <- length(
  primary_panel
)

if (!exists("font_family")) {
  font_family <- "Arial"
}

if (!exists("default_model_colors")) {
  default_model_colors <- c(
    "XGBoost" = "#0072B2",
    "Naive Bayes" = "#D55E00",
    "KNN" = "#009E73"
  )
}

model_colors <- default_model_colors[
  top3_models
]

if (any(is.na(model_colors))) {
  
  fallback_colors <- c(
    "#0072B2",
    "#D55E00",
    "#009E73"
  )
  
  model_colors <- setNames(
    fallback_colors[
      seq_along(top3_models)
    ],
    top3_models
  )
}

panel_size_performance <- panel_size_performance %>%
  filter(
    Model %in% top3_models
  ) %>%
  mutate(
    Number_genes = as.integer(
      Number_genes
    ),
    AUC = as.numeric(
      AUC
    ),
    Model_maximum_AUC = as.numeric(
      Model_maximum_AUC
    ),
    Difference_from_model_maximum = as.numeric(
      Difference_from_model_maximum
    ),
    Model = factor(
      Model,
      levels = top3_models
    ),
    Panel_label = paste0(
      "Top ",
      Number_genes
    ),
    Panel_label = factor(
      Panel_label,
      levels = paste0(
        "Top ",
        sort(
          unique(
            Number_genes
          )
        )
      )
    ),
    Primary_status = ifelse(
      Number_genes == primary_panel_size,
      "Corrected primary panel",
      "Other panel sizes"
    ),
    Primary_status = factor(
      Primary_status,
      levels = c(
        "Corrected primary panel",
        "Other panel sizes"
      )
    )
  )

panel_size_summary <- panel_size_performance %>%
  group_by(
    Number_genes,
    Genes
  ) %>%
  summarise(
    Mean_AUC = mean(
      AUC,
      na.rm = TRUE
    ),
    SD_AUC = sd(
      AUC,
      na.rm = TRUE
    ),
    Min_AUC = min(
      AUC,
      na.rm = TRUE
    ),
    Max_AUC = max(
      AUC,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Number_genes
  ) %>%
  mutate(
    Delta_from_previous = Mean_AUC -
      dplyr::lag(
        Mean_AUC
      ),
    Delta_from_previous = ifelse(
      is.na(Delta_from_previous),
      0,
      Delta_from_previous
    ),
    Best_mean_AUC = max(
      Mean_AUC,
      na.rm = TRUE
    ),
    Difference_from_best_mean = Best_mean_AUC -
      Mean_AUC,
    Within_0_01_of_best_mean =
      Difference_from_best_mean <= 0.01,
    Primary_status = ifelse(
      Number_genes == primary_panel_size,
      "Corrected primary panel",
      "Other panel sizes"
    ),
    Primary_status = factor(
      Primary_status,
      levels = c(
        "Corrected primary panel",
        "Other panel sizes"
      )
    ),
    Panel_label = paste0(
      "Top ",
      Number_genes
    ),
    Panel_label = factor(
      Panel_label,
      levels = paste0(
        "Top ",
        Number_genes
      )
    ),
    Genes_wrapped = wrap_s4_text(
      Genes,
      width = 38
    )
  )

best_panel_summary <- panel_size_summary[
  which.max(
    panel_size_summary$Mean_AUC
  ),
  ,
  drop = FALSE
]

best_panel_size <- best_panel_summary$Number_genes[1]
best_panel_auc <- best_panel_summary$Mean_AUC[1]

primary_panel_mean_auc <- panel_size_summary %>%
  filter(
    Number_genes == primary_panel_size
  ) %>%
  pull(
    Mean_AUC
  )

if (length(primary_panel_mean_auc) == 0) {
  primary_panel_mean_auc <- NA_real_
}

primary_panel_label <- paste0(
  "Primary panel = top ",
  primary_panel_size,
  " genes"
)

# ------------------------------------------------------------------------------
# S4.5 Theme
# ------------------------------------------------------------------------------

s4_theme <- theme_classic(
  base_size = 11,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "gray30",
      margin = margin(
        b = 8
      )
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 9.0
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

primary_colors <- c(
  "Corrected primary panel" = "#1B9E77",
  "Other panel sizes" = "#737373"
)

# ==============================================================================
# PANEL A
# MODEL-SPECIFIC AUC ACROSS PANEL SIZES
# ==============================================================================

p_s4_model_auc <- ggplot(
  panel_size_performance,
  aes(
    x = Number_genes,
    y = AUC,
    color = Model,
    group = Model
  )
) +
  annotate(
    "rect",
    xmin = primary_panel_size - 0.5,
    xmax = primary_panel_size + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#1B9E77",
    alpha = 0.08
  ) +
  geom_line(
    linewidth = 0.95
  ) +
  geom_point(
    size = 3.2
  ) +
  geom_vline(
    xintercept = primary_panel_size,
    linetype = "dashed",
    color = "#1B9E77",
    linewidth = 0.75
  ) +
  scale_color_manual(
    values = model_colors
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(
        panel_size_performance$Number_genes
      )
    )
  ) +
  scale_y_continuous(
    limits = c(
      max(
        0.70,
        min(
          panel_size_performance$AUC,
          na.rm = TRUE
        ) - 0.03
      ),
      min(
        1.00,
        max(
          panel_size_performance$AUC,
          na.rm = TRUE
        ) + 0.04
      )
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    title = "A. Model-specific panel-size sensitivity",
    subtitle = paste0(
      "Genes are added according to normalized consensus rank; ",
      primary_panel_label
    ),
    x = "Number of top-ranked genes included",
    y = "Outer held-out ROC-AUC",
    color = NULL
  ) +
  s4_theme +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# PANEL B
# MEAN AUC ACROSS THREE MODELS
# ==============================================================================

p_s4_mean_auc <- ggplot(
  panel_size_summary,
  aes(
    x = Number_genes,
    y = Mean_AUC
  )
) +
  annotate(
    "rect",
    xmin = primary_panel_size - 0.5,
    xmax = primary_panel_size + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#1B9E77",
    alpha = 0.08
  ) +
  geom_errorbar(
    aes(
      ymin = Min_AUC,
      ymax = Max_AUC
    ),
    width = 0.12,
    linewidth = 0.65,
    color = "gray35"
  ) +
  geom_line(
    linewidth = 0.90,
    color = "gray35"
  ) +
  geom_point(
    aes(
      fill = Primary_status
    ),
    shape = 21,
    size = 4.2,
    color = "black",
    stroke = 0.65
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Mean_AUC
      )
    ),
    nudge_y = 0.014,
    size = 3.0,
    family = font_family,
    check_overlap = TRUE
  ) +
  geom_vline(
    xintercept = primary_panel_size,
    linetype = "dashed",
    color = "#1B9E77",
    linewidth = 0.75
  ) +
  scale_fill_manual(
    values = primary_colors,
    labels = c(
      "Corrected primary panel" = "Primary panel",
      "Other panel sizes" = "Other panel sizes"
    )
  ) +
  scale_x_continuous(
    breaks = panel_size_summary$Number_genes
  ) +
  scale_y_continuous(
    limits = c(
      max(
        0.70,
        min(
          panel_size_summary$Min_AUC,
          na.rm = TRUE
        ) - 0.03
      ),
      min(
        1.00,
        max(
          panel_size_summary$Max_AUC,
          na.rm = TRUE
        ) + 0.05
      )
    ),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    title = "B. Mean AUC across the three final classifiers",
    subtitle = paste0(
      "Best mean AUC at top ",
      best_panel_size,
      " genes: ",
      sprintf(
        "%.3f",
        best_panel_auc
      )
    ),
    x = "Number of top-ranked genes included",
    y = "Mean ROC-AUC",
    fill = NULL
  ) +
  s4_theme +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# PANEL C
# MODEL-BY-PANEL-SIZE AUC HEATMAP
# ==============================================================================

heatmap_panel_data <- panel_size_performance %>%
  mutate(
    Model = factor(
      Model,
      levels = rev(
        top3_models
      )
    ),
    Panel_label = factor(
      paste0(
        "Top ",
        Number_genes
      ),
      levels = paste0(
        "Top ",
        sort(
          unique(
            Number_genes
          )
        )
      )
    ),
    AUC_label = sprintf(
      "%.3f",
      AUC
    )
  )

p_s4_heatmap <- ggplot(
  heatmap_panel_data,
  aes(
    x = Panel_label,
    y = Model,
    fill = AUC
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = AUC_label
    ),
    family = font_family,
    size = 3.0,
    color = "black"
  ) +
  geom_vline(
    xintercept = primary_panel_size + 0.5,
    linetype = "dashed",
    color = "#1B9E77",
    linewidth = 0.70
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2C7FB8",
    limits = c(
      min(
        heatmap_panel_data$AUC,
        na.rm = TRUE
      ),
      max(
        heatmap_panel_data$AUC,
        na.rm = TRUE
      )
    ),
    name = "ROC-AUC"
  ) +
  labs(
    title = "C. Model-by-panel-size ROC-AUC",
    subtitle = "Dashed line marks the five-gene primary panel",
    x = NULL,
    y = NULL
  ) +
  s4_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

# ==============================================================================
# PANEL D
# DIFFERENCE FROM BEST MEAN PANEL
# ==============================================================================

p_s4_difference <- ggplot(
  panel_size_summary,
  aes(
    x = Difference_from_best_mean,
    y = Panel_label,
    fill = Primary_status
  )
) +
  geom_col(
    width = 0.70,
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Difference_from_best_mean
      )
    ),
    hjust = -0.12,
    size = 3.0,
    family = font_family
  ) +
  scale_fill_manual(
    values = primary_colors,
    labels = c(
      "Corrected primary panel" = "Primary panel",
      "Other panel sizes" = "Other panel sizes"
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(
        panel_size_summary$Difference_from_best_mean,
        na.rm = TRUE
      ) * 1.25 + 0.005
    ),
    labels = scales::label_number(
      accuracy = 0.001
    )
  ) +
  labs(
    title = "D. Difference from best mean panel performance",
    subtitle = "Lower values indicate performance closer to the best mean ROC-AUC",
    x = "Best mean ROC-AUC - panel mean ROC-AUC",
    y = NULL,
    fill = NULL
  ) +
  s4_theme +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# S4.6 Assemble and save Supplementary Figure S4
# ------------------------------------------------------------------------------

supplementary_figure_s4 <- (
  p_s4_model_auc |
    p_s4_mean_auc
) / (
  p_s4_heatmap |
    p_s4_difference
) +
  plot_layout(
    heights = c(
      1.00,
      1.05
    ),
    widths = c(
      1.08,
      1.00
    )
  )

print(
  supplementary_figure_s4
)

save_publication_plot(
  plot = supplementary_figure_s4,
  filename_base =
    "Figure_S4_Panel_Size_Sensitivity",
  width = 18,
  height = 11.5,
  dpi = 600
)

write.csv(
  panel_size_performance,
  file.path(
    figure_dir,
    "Supplementary_Figure_S4_Model_Level_Panel_Size_Performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  panel_size_summary,
  file.path(
    figure_dir,
    "Supplementary_Figure_S4_Mean_Panel_Size_Performance.csv"
  ),
  row.names = FALSE
)

write.csv(
  heatmap_panel_data,
  file.path(
    figure_dir,
    "Supplementary_Figure_S4_Heatmap_Data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nCorrected Supplementary Figure S4 generated successfully.\n",
  "Output directory:\n",
  figure_dir,
  "\n"
)

cat(
  "\nBest mean panel-size result:\n",
  "Top ",
  best_panel_size,
  " genes; mean AUC = ",
  sprintf(
    "%.6f",
    best_panel_auc
  ),
  "\n",
  sep = ""
)

cat(
  "\nCorrected primary panel:\n",
  paste(
    primary_panel,
    collapse = ", "
  ),
  "\nMean AUC of corrected primary panel = ",
  sprintf(
    "%.6f",
    primary_panel_mean_auc
  ),
  "\n",
  sep = ""
)

cat(
  "\nPanel-size summary:\n"
)

print(
  panel_size_summary %>%
    select(
      Number_genes,
      Genes,
      Mean_AUC,
      SD_AUC,
      Min_AUC,
      Max_AUC,
      Difference_from_best_mean,
      Within_0_01_of_best_mean
    ),
  row.names = FALSE,
  digits = 8
)



















# ==============================================================================
# SUPPLEMENTARY FIGURE S5
# GRAPHICAL SUMMARY OF DEG GO AND KEGG ENRICHMENT
# Add this block after the Supplementary Figure S4 section in code 4.
# The input gene list is DEGS.txt, produced by code 2 from the limma results.
# ==============================================================================

cat("\nGenerating Supplementary Figure S5: DEG GO and KEGG enrichment\n")

# S5.1 Required packages

required_s5_packages <- c(
  "enrichR",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "stringr",
  "scales",
  "viridis"
)

missing_s5_packages <- required_s5_packages[
  !vapply(
    required_s5_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_s5_packages) > 0) {
  stop(
    "Install these packages before generating Supplementary Figure S5:\n",
    paste(missing_s5_packages, collapse = ", "),
    "\nExample:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_s5_packages), collapse = ", "),
    "))"
  )
}

library(enrichR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(stringr)
library(scales)
library(viridis)

# S5.2 Input genes from code 2

deg_file <- file.path(
  project_dir,
  "DEGS.txt"
)

if (!file.exists(deg_file)) {
  stop(
    "DEGS.txt was not found:\n",
    deg_file,
    "\nRun code 2 first. Code 2 should create DEGS.txt in project_dir."
  )
}

deg_genes <- readLines(
  deg_file,
  warn = FALSE,
  encoding = "UTF-8"
)

deg_genes <- unique(
  trimws(
    deg_genes
  )
)

deg_genes <- deg_genes[
  !is.na(deg_genes) &
    nzchar(deg_genes)
]

deg_genes <- gsub(
  '^"|"$',
  "",
  deg_genes
)

if (length(deg_genes) == 0) {
  stop("No valid gene symbols were read from DEGS.txt.")
}

if (length(deg_genes) != 74) {
  warning(
    "The manuscript reports 74 DEGs, but DEGS.txt contains ",
    length(deg_genes),
    " unique nonempty gene symbols. Check DEGS.txt before using the figure."
  )
}

write.csv(
  data.frame(
    Gene = deg_genes,
    stringsAsFactors = FALSE
  ),
  file.path(
    figure_dir,
    "Supplementary_Figure_S5_Input_DEGs.csv"
  ),
  row.names = FALSE
)

# S5.3 Select available Enrichr databases

available_enrichr_databases <- enrichR::listEnrichrDbs()

if (
  is.null(available_enrichr_databases) ||
  !"libraryName" %in% colnames(available_enrichr_databases)
) {
  stop(
    "The Enrichr database list could not be retrieved. ",
    "Check the internet connection and try again."
  )
}

available_database_names <- as.character(
  available_enrichr_databases$libraryName
)

select_enrichr_database <- function(
    candidate_names,
    category_name
) {
  
  selected_name <- candidate_names[
    candidate_names %in% available_database_names
  ]
  
  if (length(selected_name) == 0) {
    stop(
      "No supported Enrichr database was found for ",
      category_name,
      ". Tried:\n",
      paste(candidate_names, collapse = ", ")
    )
  }
  
  selected_name[1]
}

database_map <- c(
  "GO Biological Process" = select_enrichr_database(
    c(
      "GO_Biological_Process_2025",
      "GO_Biological_Process_2023",
      "GO_Biological_Process_2021"
    ),
    "GO Biological Process"
  ),
  "GO Cellular Component" = select_enrichr_database(
    c(
      "GO_Cellular_Component_2025",
      "GO_Cellular_Component_2023",
      "GO_Cellular_Component_2021"
    ),
    "GO Cellular Component"
  ),
  "GO Molecular Function" = select_enrichr_database(
    c(
      "GO_Molecular_Function_2025",
      "GO_Molecular_Function_2023",
      "GO_Molecular_Function_2021"
    ),
    "GO Molecular Function"
  ),
  "KEGG pathways" = select_enrichr_database(
    c(
      "KEGG_2021_Human",
      "KEGG_2019_Human",
      "KEGG_2016"
    ),
    "KEGG pathways"
  )
)

write.csv(
  data.frame(
    Category = names(database_map),
    Enrichr_database = unname(database_map),
    stringsAsFactors = FALSE
  ),
  file.path(
    figure_dir,
    "Supplementary_Figure_S5_Databases_Used.csv"
  ),
  row.names = FALSE
)

cat("\nEnrichr databases used:\n")
print(database_map)

# S5.4 Run enrichment

enrichment_results_raw <- tryCatch(
  enrichR::enrichr(
    deg_genes,
    databases = unname(database_map)
  ),
  error = function(error_object) {
    stop(
      "Enrichr analysis failed. Check the internet connection.\n",
      conditionMessage(error_object)
    )
  }
)

standardize_enrichr_result <- function(
    enrichment_table,
    category_name,
    database_name
) {
  
  if (
    is.null(enrichment_table) ||
    nrow(enrichment_table) == 0
  ) {
    return(
      data.frame(
        Category = character(0),
        Database = character(0),
        Term = character(0),
        Overlap = character(0),
        Overlap_Count = numeric(0),
        Gene_Set_Size = numeric(0),
        Gene_Ratio = numeric(0),
        P_value = numeric(0),
        Adjusted_P = numeric(0),
        Odds_Ratio = numeric(0),
        Combined_Score = numeric(0),
        Genes = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  required_columns <- c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(enrichment_table)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "The Enrichr result for ",
      category_name,
      " is missing columns:\n",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  enrichment_table %>%
    dplyr::transmute(
      Category = category_name,
      Database = database_name,
      Term = as.character(.data$Term),
      Overlap = as.character(.data$Overlap),
      Overlap_Count = suppressWarnings(
        as.numeric(
          sub(
            "/.*$",
            "",
            .data$Overlap
          )
        )
      ),
      Gene_Set_Size = suppressWarnings(
        as.numeric(
          sub(
            "^.*/",
            "",
            .data$Overlap
          )
        )
      ),
      Gene_Ratio = dplyr::if_else(
        is.finite(.data$Overlap_Count) &
          is.finite(.data$Gene_Set_Size) &
          .data$Gene_Set_Size > 0,
        .data$Overlap_Count / .data$Gene_Set_Size,
        NA_real_
      ),
      P_value = as.numeric(.data$P.value),
      Adjusted_P = as.numeric(.data$Adjusted.P.value),
      Odds_Ratio = as.numeric(.data$Odds.Ratio),
      Combined_Score = as.numeric(.data$Combined.Score),
      Genes = as.character(.data$Genes)
    )
}

enrichment_all <- dplyr::bind_rows(
  lapply(
    names(database_map),
    function(category_name) {
      
      database_name <- unname(
        database_map[
          category_name
        ]
      )
      
      standardize_enrichr_result(
        enrichment_table = enrichment_results_raw[[database_name]],
        category_name = category_name,
        database_name = database_name
      )
    }
  )
)

write.csv(
  enrichment_all,
  file.path(
    figure_dir,
    "Supplementary_Figure_S5_All_Enrichment_Results.csv"
  ),
  row.names = FALSE
)

# S5.5 Select the plotted terms

number_terms_per_panel <- 8
fdr_cutoff <- 0.05

enrichment_significant <- enrichment_all %>%
  dplyr::filter(
    is.finite(.data$Adjusted_P),
    .data$Adjusted_P < fdr_cutoff,
    is.finite(.data$Combined_Score),
    is.finite(.data$Overlap_Count),
    .data$Overlap_Count > 0
  ) %>%
  dplyr::group_by(
    .data$Category
  ) %>%
  dplyr::arrange(
    .data$Adjusted_P,
    dplyr::desc(.data$Combined_Score),
    .by_group = TRUE
  ) %>%
  dplyr::distinct(
    .data$Term,
    .keep_all = TRUE
  ) %>%
  dplyr::slice_head(
    n = number_terms_per_panel
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Minus_log10_FDR = -log10(
      pmax(
        .data$Adjusted_P,
        .Machine$double.xmin
      )
    ),
    Term_clean = gsub(
      "\\s*\\(GO:[0-9]+\\)\\s*$",
      "",
      .data$Term
    ),
    Term_clean = gsub(
      "_Homo sapiens$",
      "",
      .data$Term_clean
    )
  )

categories_without_significant_terms <- setdiff(
  names(database_map),
  unique(enrichment_significant$Category)
)

if (length(categories_without_significant_terms) > 0) {
  warning(
    "No terms passed BH FDR < 0.05 for:\n",
    paste(categories_without_significant_terms, collapse = ", "),
    "\nThe corresponding panel will state that no term passed the threshold."
  )
}

write.csv(
  enrichment_significant,
  file.path(
    figure_dir,
    "Supplementary_Figure_S5_Plotted_Enrichment_Terms.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# REPLACEMENT BLOCK FOR SUPPLEMENTARY FIGURE S5
# Replace the previous S5.6 and S5.7 sections with this block.
# ==============================================================================

# S5.6 Prepare transformed values and shared visual limits

enrichment_significant <- enrichment_significant %>%
  dplyr::mutate(
    Log10_Combined_Score = log10(
      pmax(
        .data$Combined_Score,
        .Machine$double.xmin
      )
    )
  )

global_fill_limits <- range(
  enrichment_significant$Minus_log10_FDR,
  na.rm = TRUE
)

global_size_limits <- range(
  enrichment_significant$Overlap_Count,
  na.rm = TRUE
)

if (
  length(global_fill_limits) != 2 ||
  any(!is.finite(global_fill_limits)) ||
  diff(global_fill_limits) == 0
) {
  global_fill_limits <- c(0, 1)
}

if (
  length(global_size_limits) != 2 ||
  any(!is.finite(global_size_limits)) ||
  diff(global_size_limits) == 0
) {
  global_size_limits <- c(1, 2)
}

s5_theme <- theme_classic(
  base_size = 12,
  base_family = font_family
) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 12
    ),
    axis.text.x = element_text(
      color = "black",
      size = 10.5
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10.5,
      lineheight = 0.92
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size = 10.2,
      color = "gray30",
      margin = margin(b = 7)
    ),
    legend.title = element_text(
      face = "bold",
      size = 10.5
    ),
    legend.text = element_text(
      size = 10
    ),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.40
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 14,
      b = 10,
      l = 10
    )
  )

make_enrichment_panel <- function(
    category_name,
    panel_title
) {
  
  panel_data <- enrichment_significant %>%
    dplyr::filter(
      .data$Category == category_name
    ) %>%
    dplyr::arrange(
      .data$Log10_Combined_Score
    ) %>%
    dplyr::mutate(
      Term_wrapped = stringr::str_wrap(
        .data$Term_clean,
        width = 28
      ),
      Term_wrapped = factor(
        .data$Term_wrapped,
        levels = unique(.data$Term_wrapped)
      )
    )
  
  if (nrow(panel_data) == 0) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No term passed BH FDR < 0.05",
          family = font_family,
          size = 4.4
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
          title = panel_title
        ) +
        theme_void(
          base_family = font_family
        ) +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 14,
            hjust = 0
          )
        )
    )
  }
  
  ggplot(
    panel_data,
    aes(
      x = .data$Log10_Combined_Score,
      y = .data$Term_wrapped,
      size = .data$Overlap_Count,
      fill = .data$Minus_log10_FDR
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.50,
      alpha = 0.94
    ) +
    scale_size_continuous(
      limits = global_size_limits,
      range = c(4.0, 9.2),
      breaks = scales::pretty_breaks(
        n = 4
      )
    ) +
    scale_fill_viridis_c(
      option = "C",
      limits = global_fill_limits,
      end = 0.92
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 5
      ),
      labels = scales::label_number(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(0.03, 0.10)
      )
    ) +
    labs(
      title = panel_title,
      subtitle = paste0(
        "Top ",
        min(number_terms_per_panel, nrow(panel_data)),
        " terms with BH FDR < 0.05"
      ),
      x = expression(log[10]("Enrichr combined score")),
      y = NULL,
      size = "Gene count",
      fill = expression(-log[10]("adjusted p value"))
    ) +
    s5_theme
}

p_s5_bp <- make_enrichment_panel(
  "GO Biological Process",
  "A. GO Biological Process"
)

p_s5_cc <- make_enrichment_panel(
  "GO Cellular Component",
  "B. GO Cellular Component"
)

p_s5_mf <- make_enrichment_panel(
  "GO Molecular Function",
  "C. GO Molecular Function"
)

p_s5_kegg <- make_enrichment_panel(
  "KEGG pathways",
  "D. KEGG pathways"
)

# S5.7 Assemble and save the revised figure

supplementary_figure_s5 <- (
  p_s5_bp |
    p_s5_cc
) / (
  p_s5_mf |
    p_s5_kegg
) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(1, 1),
    widths = c(1, 1)
  ) &
  theme(
    legend.position = "bottom"
  )

print(
  supplementary_figure_s5
)

save_publication_plot(
  plot = supplementary_figure_s5,
  filename_base =
    "Figure_S5_Enrichment",
  width = 18,
  height = 13,
  dpi = 600
)

cat(
  "\nRevised Supplementary Figure S5 generated successfully.\n",
  "The x axis now shows log10 transformed Enrichr combined scores.\n",
  sep = ""
)