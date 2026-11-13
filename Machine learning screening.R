# ==============================================================================
# INTEGRATED REPRODUCIBLE ML PIPELINE
# FA vs FVPTC | Stage 1 screening + Stage 2 final top-3 modeling + complete one-file report
# ==============================================================================

# ---------------------------
# 0. REPRODUCIBILITY SETTINGS
# ---------------------------

SEED_GLOBAL <- 123
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(SEED_GLOBAL)

# Avoid hidden parallel/thread-related variation
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

# Keep cache disabled for the screening stage
use_cache <- FALSE

# ---------------------------
# 1. LIBRARIES
# ---------------------------

library(caret)
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)
library(randomForest)
library(kernlab)
library(xgboost)

if (!requireNamespace("naivebayes", quietly = TRUE)) {
  stop("Package 'naivebayes' is required for caret method='naive_bayes'. Install it first.")
}
library(naivebayes)

# Force sequential backend if foreach/doParallel are present
if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()
}

# ---------------------------
# 2. PATHS AND OUTPUT FOLDER
# ---------------------------

project_dir <- getwd()

# All stage-1 outputs go here
outdir <- file.path(project_dir, "ml-output", "stage1_screening_reproducible")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(outdir)) {
  stop("Output directory could not be created: ", outdir)
}

out_file <- function(filename) {
  file.path(outdir, filename)
}

# Check write permission before doing any modeling
test_file <- out_file(".write_test.tmp")
can_write <- tryCatch({
  ok <- file.create(test_file)
  if (ok && file.exists(test_file)) file.remove(test_file)
  ok
}, warning = function(w) FALSE, error = function(e) FALSE)

if (!isTRUE(can_write)) {
  stop(
    "R cannot write to the output directory: ", outdir,
    "\nClose files in this folder, check OneDrive sync/permissions, or choose a local folder."
  )
}

save_figure_outputs <- function(plot, filename_base, width, height, dpi = 300, bg = "white") {
  ggsave(
    filename = out_file(paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg
  )
  
  ggsave(
    filename = out_file(paste0(filename_base, ".tiff")),
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = bg
  )
}

cat("\nStage-1 screening outputs will be saved in:\n")
cat(outdir, "\n\n")

# ---------------------------
# 3. PREPARE DATASET
# ---------------------------

# Use final_data from the current R session if code 2 has already created it.
# Otherwise, try to read the saved RDS created by code 2.
if (!exists("final_data")) {
  final_data_rds <- file.path(project_dir, "final_data_clean.rds")
  
  if (!file.exists(final_data_rds)) {
    stop(
      "Object 'final_data' was not found in the R session, and this file was not found: ",
      final_data_rds,
      "\nRun code 1 and code 2 first, or save final_data as final_data_clean.rds."
    )
  }
  
  final_data <- readRDS(final_data_rds)
}

gene_list <- c(
  "FN1", "KRT19", "TIMP1", "CHI3L1",
  "CLU", "SPP1", "APOE", "EFEMP1"
)

Groups <- c(
  rep("FA", 17),
  rep("FVPTC", 13),
  rep("FA", 7),
  rep("FVPTC", 6)
)

missing_genes <- setdiff(gene_list, rownames(final_data))
if (length(missing_genes) > 0) {
  stop(
    paste(
      "These genes are missing from final_data:",
      paste(missing_genes, collapse = ", ")
    )
  )
}

if (ncol(final_data) != length(Groups)) {
  stop(
    "Number of columns in final_data does not match number of manual group labels. ",
    "ncol(final_data) = ", ncol(final_data),
    "; length(Groups) = ", length(Groups)
  )
}

gene_data <- final_data[gene_list, , drop = FALSE]
gene_data <- as.data.frame(t(gene_data), stringsAsFactors = FALSE)
gene_data$label <- factor(Groups, levels = c("FA", "FVPTC"))

for (col in colnames(gene_data)[!colnames(gene_data) %in% "label"]) {
  gene_data[[col]] <- as.numeric(gene_data[[col]])
}

data_check <- data.frame(
  Total_samples = nrow(gene_data),
  Total_genes = length(gene_list),
  FA = sum(gene_data$label == "FA"),
  FVPTC = sum(gene_data$label == "FVPTC")
)

write.csv(data_check, out_file("Report_Data_Check.csv"), row.names = FALSE)

sample_label_check <- data.frame(
  Sample_order = seq_len(nrow(gene_data)),
  Sample_name = rownames(gene_data),
  Label = gene_data$label,
  stringsAsFactors = FALSE
)

write.csv(sample_label_check, out_file("Report_Sample_Label_Order_Used_By_ML.csv"), row.names = FALSE)

cat("\n=======================================================\n")
cat("--- STAGE 1: DATA PREPARATION COMPLETED ---\n")
cat("Total samples:", nrow(gene_data), "\n")
cat("Total genes analyzed:", length(gene_list), "\n")
cat("Group distribution:\n")
print(table(gene_data$label))
cat("=======================================================\n")

# ---------------------------
# 4. REPRODUCIBLE INITIAL SCREENING OF 7 MODELS
# ---------------------------

model_methods <- list(
  "Logistic Regression" = "glm",
  "SVM" = "svmLinear",
  "Decision Tree" = "rpart",
  "Random Forest" = "rf",
  "Naive Bayes" = "naive_bayes",
  "KNN" = "knn",
  "XGBoost" = "xgbTree"
)

n_runs <- 10
k_folds <- 5

all_results <- data.frame()
all_resamples <- data.frame()
all_models <- list()

cat("\n=======================================================\n")
cat("--- STAGE 2: REPRODUCIBLE INITIAL 7-MODEL SCREENING ---\n")
cat("Runs:", n_runs, "| Folds per run:", k_folds, "\n")
cat("Cache disabled:", use_cache, "\n")
cat("=======================================================\n")

for (i in seq_len(n_runs)) {
  
  cat("\nRunning screening run", i, "of", n_runs, "\n")
  
  set.seed(SEED_GLOBAL + i)
  
  # Same folds are used for all models within this run.
  run_index <- caret::createFolds(
    gene_data$label,
    k = k_folds,
    returnTrain = TRUE
  )
  
  run_index_out <- lapply(
    run_index,
    function(train_idx) {
      setdiff(seq_len(nrow(gene_data)), train_idx)
    }
  )
  
  train_control_run <- trainControl(
    method = "cv",
    index = run_index,
    indexOut = run_index_out,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = FALSE
  )
  
  all_models[[i]] <- list()
  
  for (m in names(model_methods)) {
    
    cat("  Model:", m, "\n")
    
    set.seed(SEED_GLOBAL + i * 100 + match(m, names(model_methods)))
    
    # Keep initial screening close to the original code:
    # no model-specific custom tuning grids are imposed at this stage.
    model <- train(
      label ~ .,
      data = gene_data,
      method = model_methods[[m]],
      trControl = train_control_run,
      metric = "ROC"
    )
    
    res <- model$resample
    res$Run <- i
    res$Model <- m
    
    all_resamples <- rbind(
      all_resamples,
      res[, c("Run", "Model", setdiff(colnames(res), c("Run", "Model")))]
    )
    
    all_results <- rbind(
      all_results,
      data.frame(
        Run = i,
        Model = m,
        ROC = mean(model$resample$ROC, na.rm = TRUE),
        Sens = mean(model$resample$Sens, na.rm = TRUE),
        Spec = mean(model$resample$Spec, na.rm = TRUE),
        ROC_SD_within_run = sd(model$resample$ROC, na.rm = TRUE),
        BestTune = paste(names(model$bestTune), model$bestTune[1, ], sep = "=", collapse = "; "),
        stringsAsFactors = FALSE
      )
    )
    
    all_models[[i]][[m]] <- model
  }
}

# ---------------------------
# 5. NUMERIC REPORTS FOR SCREENING FIGURE
# ---------------------------

summary_results <- all_results %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    MeanROC = mean(ROC, na.rm = TRUE),
    SDROC = sd(ROC, na.rm = TRUE),
    MedianROC = median(ROC, na.rm = TRUE),
    IQRROC = IQR(ROC, na.rm = TRUE),
    MinROC = min(ROC, na.rm = TRUE),
    MaxROC = max(ROC, na.rm = TRUE),
    MeanSensitivity = mean(Sens, na.rm = TRUE),
    MeanSpecificity = mean(Spec, na.rm = TRUE),
    N_runs = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(MeanROC)) %>%
  dplyr::mutate(Screening_rank = dplyr::row_number())

top3_candidates <- summary_results %>%
  dplyr::slice_head(n = 3)

write.csv(
  all_results,
  out_file("Table_S0_Initial_7_Models_Run_Level_ROC.csv"),
  row.names = FALSE
)

write.csv(
  all_resamples,
  out_file("Table_S0_Initial_7_Models_Fold_Level_Resamples.csv"),
  row.names = FALSE
)

write.csv(
  summary_results,
  out_file("Table_S0_Initial_7_Models_Numeric_Summary.csv"),
  row.names = FALSE
)

write.csv(
  top3_candidates,
  out_file("Table_S0_Top3_Model_Candidates_FOR_REVIEW.csv"),
  row.names = FALSE
)

saveRDS(
  all_models,
  out_file("Stage1_All_Initial_Screening_Models.rds")
)

saveRDS(
  list(
    all_results = all_results,
    all_resamples = all_resamples,
    summary_results = summary_results,
    top3_candidates = top3_candidates,
    gene_list = gene_list,
    Groups = Groups,
    seed = SEED_GLOBAL,
    n_runs = n_runs,
    k_folds = k_folds
  ),
  out_file("Stage1_Screening_Results_Object.rds")
)

cat("\n--- SUMMARY OF INITIAL 7 MODELS ---\n")
print(summary_results)

cat("\n--- TOP 3 MODEL CANDIDATES; DO NOT CONTINUE UNTIL REVIEWED ---\n")
print(top3_candidates)

# ---------------------------
# 6. FIGURE S0 + NUMERIC REPORT
# ---------------------------

p_initial_roc <- ggplot(
  all_results,
  aes(x = Model, y = ROC, fill = Model)
) +
  geom_violin(trim = FALSE, alpha = 0.25, color = NA) +
  geom_jitter(
    aes(color = Model),
    width = 0.16,
    size = 2,
    alpha = 0.55
  ) +
  # ایجاد خط مشکی نازک برای نمایش میانگین
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.4,       # عرض خط میانگین
    color = "black",   # رنگ خط
    linewidth = 0.2    # ضخامت خط (نازک)
  ) +                  # <--- این علامت + در کد اصلی شما جا افتاده بود
  labs(
    title = "ROC (AUC) Distribution per Model (Initial Screening)",
    subtitle = paste0(
      "Reproducible screening: ", n_runs,
      " runs of stratified ", k_folds, "-fold CV"
    ),
    x = NULL,
    y = "ROC (AUC)"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "none"
  )

save_figure_outputs(
  plot = p_initial_roc,
  filename_base = "Supplementary_Figure_S0_Initial_Screening",
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


screening_report_lines <- c(
  "NUMERIC REPORT FOR SUPPLEMENTARY FIGURE S0",
  "==========================================",
  paste0("Seed: ", SEED_GLOBAL),
  paste0("RNGkind: ", paste(RNGkind(), collapse = " | ")),
  paste0("Runs: ", n_runs),
  paste0("Folds per run: ", k_folds),
  paste0("Total samples: ", nrow(gene_data)),
  paste0("Group distribution: FA=", sum(gene_data$label == "FA"),
         ", FVPTC=", sum(gene_data$label == "FVPTC")),
  "",
  "Model ranking by mean ROC-AUC:",
  paste(
    summary_results$Screening_rank,
    summary_results$Model,
    paste0("MeanROC=", round(summary_results$MeanROC, 4)),
    paste0("SD=", round(summary_results$SDROC, 4)),
    paste0("Median=", round(summary_results$MedianROC, 4)),
    paste0("Range=", round(summary_results$MinROC, 4), "-", round(summary_results$MaxROC, 4)),
    sep = ". "
  ),
  "",
  "Top 3 candidate models for review before final modeling:",
  paste(
    top3_candidates$Screening_rank,
    top3_candidates$Model,
    paste0("MeanROC=", round(top3_candidates$MeanROC, 4)),
    sep = ". "
  ),
  "",
  "Integrated run note: Stage 1 screening is completed and the script continues automatically to Stage 2 final modeling.",
  "Stage 1 key output files:",
  "- Table_S0_Initial_7_Models_Numeric_Summary.csv",
  "- Table_S0_Top3_Model_Candidates_FOR_REVIEW.csv",
  "- Supplementary_Figure_S0_Initial_Screening.png or .tiff"
)

writeLines(
  screening_report_lines,
  con = out_file("Report_Supplementary_Figure_S0_Initial_Screening.txt")
)

capture.output(
  sessionInfo(),
  file = out_file("SessionInfo_Stage1_Screening.txt")
)

# ---------------------------
# 7. STOP POINT
# ---------------------------

cat("\n=======================================================\n")
cat(">>> STAGE 1 SCREENING COMPLETED SUCCESSFULLY <<<\n")
cat(">>> STAGE 1 SCREENING COMPLETED; CONTINUING TO STAGE 2 FINAL MODELING <<<\n")
cat("=======================================================\n")
cat("Outputs saved in:\n", outdir, "\n")
cat("\nStage 1 key files generated:\n")
cat("1) Table_S0_Initial_7_Models_Numeric_Summary.csv\n")
cat("2) Table_S0_Top3_Model_Candidates_FOR_REVIEW.csv\n")
cat("3) Supplementary_Figure_S0_Initial_Screening.png or .tiff\n")
cat("4) Report_Supplementary_Figure_S0_Initial_Screening.txt\n")
cat("=======================================================\n")















# ==============================================================================
# STAGE 2: REPRODUCIBLE FINAL MODELING USING TOP-3 SCREENING MODELS
# FA vs FVPTC | XGBoost + Naive Bayes + KNN by default from Stage-1 screening
# Generates PNG + TIFF figures and numeric reports for every figure.
# ==============================================================================

# ---------------------------
# 0. REPRODUCIBILITY SETTINGS
# ---------------------------

SEED_GLOBAL <- 123
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(SEED_GLOBAL)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

use_cache <- FALSE

# ---------------------------
# 1. LIBRARIES
# ---------------------------

library(caret)
library(ggplot2)
library(pROC)
library(patchwork)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggpubr)
library(kernlab)
library(xgboost)

if (!requireNamespace("naivebayes", quietly = TRUE)) {
  stop("Package 'naivebayes' is required. Install it first: install.packages('naivebayes')")
}
library(naivebayes)

if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()
}

# ---------------------------
# 2. PATHS AND OUTPUT FOLDER
# ---------------------------
project_dir <- getwd()
stage1_dir <- file.path(project_dir, "ml-output", "stage1_screening_reproducible")
outdir <- file.path(project_dir, "ml-output", "stage2_top3_final_reproducible")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(outdir)) {
  stop("Output directory could not be created: ", outdir)
}

out_file <- function(filename) {
  file.path(outdir, filename)
}

# Check write permission
test_file <- out_file(".write_test.tmp")
can_write <- tryCatch({
  ok <- file.create(test_file)
  if (ok && file.exists(test_file)) file.remove(test_file)
  ok
}, warning = function(w) FALSE, error = function(e) FALSE)

if (!isTRUE(can_write)) {
  stop(
    "R cannot write to the output directory: ", outdir,
    "\nClose files in this folder, check OneDrive sync/permissions, or choose a local folder."
  )
}

save_figure_outputs <- function(plot, filename_base, width, height, dpi = 300, bg = "white") {
  ggsave(
    filename = out_file(paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg
  )
  
  ggsave(
    filename = out_file(paste0(filename_base, ".tiff")),
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = bg
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

cat("\nStage-2 final outputs will be saved in:\n")
cat(outdir, "\n\n")

# ---------------------------
# 3. LOAD STAGE-1 MODEL SELECTION
# ---------------------------

stage1_rds <- file.path(stage1_dir, "Stage1_Screening_Results_Object.rds")

if (file.exists(stage1_rds)) {
  stage1_obj <- readRDS(stage1_rds)
  final_model_names <- as.character(stage1_obj$top3_candidates$Model)
} else {
  warning(
    "Stage1_Screening_Results_Object.rds was not found. ",
    "Using manually specified top-3 models from the reviewed screening output."
  )
  final_model_names <- c("XGBoost", "Naive Bayes", "KNN")
}

allowed_models <- c("XGBoost", "Naive Bayes", "KNN", "Random Forest", "SVM", "Decision Tree", "Logistic Regression")
if (!all(final_model_names %in% allowed_models)) {
  stop("Unsupported model found in final_model_names: ", paste(setdiff(final_model_names, allowed_models), collapse = ", "))
}

if (length(final_model_names) != 3) {
  stop("Exactly three final models are expected. Found: ", paste(final_model_names, collapse = ", "))
}

write.csv(
  data.frame(Final_model_rank = 1:3, Model = final_model_names),
  out_file("Table_Final_Selected_Models_From_Stage1.csv"),
  row.names = FALSE
)

cat("Final models selected from Stage-1 screening:\n")
print(data.frame(Rank = 1:3, Model = final_model_names))

# ---------------------------
# 4. PREPARE DATASET
# ---------------------------

if (!exists("final_data")) {
  final_data_rds <- file.path(project_dir, "final_data_clean.rds")
  
  if (!file.exists(final_data_rds)) {
    stop(
      "Object 'final_data' was not found in the R session, and this file was not found: ",
      final_data_rds,
      "\nRun code 1 and code 2 first, or save final_data as final_data_clean.rds."
    )
  }
  
  final_data <- readRDS(final_data_rds)
}

gene_list <- c(
  "FN1", "KRT19", "TIMP1", "CHI3L1",
  "CLU", "SPP1", "APOE", "EFEMP1"
)

Groups <- c(
  rep("FA", 17),
  rep("FVPTC", 13),
  rep("FA", 7),
  rep("FVPTC", 6)
)

missing_genes <- setdiff(gene_list, rownames(final_data))
if (length(missing_genes) > 0) {
  stop("These genes are missing from final_data: ", paste(missing_genes, collapse = ", "))
}

if (ncol(final_data) != length(Groups)) {
  stop(
    "Number of columns in final_data does not match number of manual group labels. ",
    "ncol(final_data) = ", ncol(final_data), "; length(Groups) = ", length(Groups)
  )
}

gene_data <- final_data[gene_list, , drop = FALSE]
gene_data <- as.data.frame(t(gene_data), stringsAsFactors = FALSE)
gene_data$label <- factor(Groups, levels = c("FA", "FVPTC"))

for (col in colnames(gene_data)[!colnames(gene_data) %in% "label"]) {
  gene_data[[col]] <- as.numeric(gene_data[[col]])
}

write.csv(
  data.frame(
    Total_samples = nrow(gene_data),
    Total_genes = length(gene_list),
    FA = sum(gene_data$label == "FA"),
    FVPTC = sum(gene_data$label == "FVPTC")
  ),
  out_file("Report_Data_Check.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Sample_order = seq_len(nrow(gene_data)),
    Sample_name = rownames(gene_data),
    Label = gene_data$label,
    stringsAsFactors = FALSE
  ),
  out_file("Report_Sample_Label_Order_Used_By_ML.csv"),
  row.names = FALSE
)

cat("\n=======================================================\n")
cat("--- DATA PREPARATION COMPLETED ---\n")
cat("Total samples:", nrow(gene_data), "\n")
cat("Total genes analyzed:", length(gene_list), "\n")
cat("Group distribution:\n")
print(table(gene_data$label))
cat("=======================================================\n")

# ---------------------------
# 5. UNIFIED REPEATED 5-FOLD CV FOR TOP-3 FINAL MODELS
# ---------------------------

set.seed(SEED_GLOBAL)
fold_index <- createMultiFolds(gene_data$label, k = 5, times = 20)

fold_index_out <- lapply(
  fold_index,
  function(train_idx) {
    setdiff(seq_len(nrow(gene_data)), train_idx)
  }
)

unified_control <- trainControl(
  method = "cv",
  index = fold_index,
  indexOut = fold_index_out,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter = FALSE,
  allowParallel = FALSE
)

train_final_model <- function(model_name) {
  
  cat("\nTraining final model:", model_name, "\n")
  
  if (model_name == "XGBoost") {
    xgb_grid <- expand.grid(
      nrounds = c(50, 100, 150),
      max_depth = c(2, 3),
      eta = c(0.01, 0.05, 0.1),
      gamma = 0,
      colsample_bytree = c(0.8, 1.0),
      min_child_weight = 1,
      subsample = c(0.8, 1.0)
    )
    
    model <- train(
      label ~ .,
      data = gene_data,
      method = "xgbTree",
      metric = "ROC",
      trControl = unified_control,
      tuneGrid = xgb_grid,
      nthread = 1
    )
    
  } else if (model_name == "Naive Bayes") {
    nb_grid <- expand.grid(
      laplace = c(0, 1),
      usekernel = c(TRUE, FALSE),
      adjust = c(1, 2)
    )
    
    model <- train(
      label ~ .,
      data = gene_data,
      method = "naive_bayes",
      metric = "ROC",
      trControl = unified_control,
      tuneGrid = nb_grid
    )
    
  } else if (model_name == "KNN") {
    knn_grid <- expand.grid(k = seq(3, 15, by = 2))
    
    model <- train(
      label ~ .,
      data = gene_data,
      method = "knn",
      metric = "ROC",
      trControl = unified_control,
      tuneGrid = knn_grid,
      preProcess = c("center", "scale")
    )
    
  } else if (model_name == "Random Forest") {
    model <- train(
      label ~ .,
      data = gene_data,
      method = "rf",
      metric = "ROC",
      trControl = unified_control,
      tuneGrid = expand.grid(mtry = 2),
      ntree = 100
    )
    
  } else if (model_name == "SVM") {
    model <- train(
      label ~ .,
      data = gene_data,
      method = "svmLinear",
      metric = "ROC",
      trControl = unified_control,
      tuneGrid = expand.grid(C = 2^(-2:5)),
      preProcess = c("center", "scale")
    )
    
  } else if (model_name == "Decision Tree") {
    model <- train(
      label ~ .,
      data = gene_data,
      method = "rpart",
      metric = "ROC",
      trControl = unified_control
    )
    
  } else if (model_name == "Logistic Regression") {
    model <- train(
      label ~ .,
      data = gene_data,
      method = "glm",
      metric = "ROC",
      trControl = unified_control
    )
    
  } else {
    stop("Unsupported model: ", model_name)
  }
  
  model
}

final_models <- list()

for (i in seq_along(final_model_names)) {
  set.seed(SEED_GLOBAL + 1000 + i)
  model_name <- final_model_names[i]
  final_models[[model_name]] <- train_final_model(model_name)
  saveRDS(final_models[[model_name]], out_file(paste0("Final_Model_", safe_name(model_name), ".rds")))
}

best_tune_table <- dplyr::bind_rows(
  lapply(names(final_models), function(model_name) {
    data.frame(
      Model = model_name,
      final_models[[model_name]]$bestTune,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(best_tune_table, out_file("Table_Final_Model_Best_Tuning_Parameters.csv"), row.names = FALSE)

# ---------------------------
# 6. AGGREGATED PER-SAMPLE PREDICTIONS AND PERFORMANCE
# ---------------------------

get_best_preds <- function(model_object) {
  model_object$pred %>%
    dplyr::semi_join(model_object$bestTune, by = names(model_object$bestTune))
}

aggregate_sample_preds <- function(best_preds) {
  best_preds %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      obs = dplyr::first(obs),
      Mean_Probability_FVPTC = mean(FVPTC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Sample = rownames(gene_data)[rowIndex],
      pred = factor(
        ifelse(Mean_Probability_FVPTC >= 0.5, "FVPTC", "FA"),
        levels = c("FA", "FVPTC")
      )
    ) %>%
    dplyr::select(Sample, rowIndex, obs, Mean_Probability_FVPTC, pred)
}

make_cm <- function(pred, obs) {
  table(
    Predicted = factor(pred, levels = c("FA", "FVPTC")),
    Actual = factor(obs, levels = c("FA", "FVPTC"))
  )
}

calculate_metrics <- function(agg_preds, model_name) {
  cm <- make_cm(agg_preds$pred, agg_preds$obs)
  
  TP <- as.numeric(cm["FVPTC", "FVPTC"])
  TN <- as.numeric(cm["FA", "FA"])
  FP <- as.numeric(cm["FVPTC", "FA"])
  FN <- as.numeric(cm["FA", "FVPTC"])
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  f1 <- ifelse(is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0,
               NA, 2 * precision * sensitivity / (precision + sensitivity))
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denominator == 0, 0, ((TP * TN) - (FP * FN)) / denominator)
  
  roc_obj <- pROC::roc(
    agg_preds$obs,
    agg_preds$Mean_Probability_FVPTC,
    levels = c("FA", "FVPTC"),
    direction = "<",
    quiet = TRUE
  )
  
  data.frame(
    Model = model_name,
    AUC_aggregated = as.numeric(pROC::auc(roc_obj)),
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision_PPV = precision,
    NPV = npv,
    F1 = f1,
    Balanced_Accuracy = balanced_accuracy,
    MCC = mcc,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    stringsAsFactors = FALSE
  )
}

best_preds_list <- lapply(final_models, get_best_preds)
agg_preds_list <- lapply(best_preds_list, aggregate_sample_preds)

for (model_name in names(agg_preds_list)) {
  write.csv(
    agg_preds_list[[model_name]],
    out_file(paste0("Table_Aggregated_Predictions_", safe_name(model_name), ".csv")),
    row.names = FALSE
  )
}

performance_table <- dplyr::bind_rows(
  lapply(names(agg_preds_list), function(model_name) {
    calculate_metrics(agg_preds_list[[model_name]], model_name)
  })
)

write.csv(performance_table, out_file("Table_Final_Model_Performance_Aggregated.csv"), row.names = FALSE)

confusion_long <- dplyr::bind_rows(
  lapply(names(agg_preds_list), function(model_name) {
    cm <- make_cm(agg_preds_list[[model_name]]$pred, agg_preds_list[[model_name]]$obs)
    as.data.frame(cm) %>%
      dplyr::mutate(Model = model_name) %>%
      dplyr::select(Model, Predicted, Actual, Freq)
  })
)

write.csv(confusion_long, out_file("Table_Final_Model_Confusion_Matrices_Long.csv"), row.names = FALSE)

# Nadeau-Bengio corrected tests
nadeau_bengio_test <- function(res1, res2, k_folds = 5) {
  required_columns <- c("Resample", "ROC")
  
  if (!all(required_columns %in% names(res1))) stop("res1 must contain Resample and ROC columns.")
  if (!all(required_columns %in% names(res2))) stop("res2 must contain Resample and ROC columns.")
  
  merged <- dplyr::inner_join(
    res1 %>% dplyr::transmute(Resample = as.character(Resample), ROC_1 = as.numeric(ROC)),
    res2 %>% dplyr::transmute(Resample = as.character(Resample), ROC_2 = as.numeric(ROC)),
    by = "Resample"
  )
  
  differences <- merged$ROC_1 - merged$ROC_2
  differences <- differences[is.finite(differences)]
  
  if (length(differences) < 2) stop("Fewer than two finite paired ROC differences were available.")
  
  mean_diff <- mean(differences)
  var_diff <- stats::var(differences)
  
  if (is.na(var_diff) || var_diff <= .Machine$double.eps) {
    return(if (abs(mean_diff) <= .Machine$double.eps) 1 else 0)
  }
  
  total_runs <- length(differences)
  correction <- (1 / total_runs) + (1 / (k_folds - 1))
  corrected_var <- var_diff * correction
  t_stat <- mean_diff / sqrt(corrected_var)
  p_value <- 2 * stats::pt(-abs(t_stat), df = total_runs - 1)
  as.numeric(p_value)
}

pairwise_tests <- combn(names(final_models), 2, simplify = FALSE)
nadeau_table <- dplyr::bind_rows(
  lapply(pairwise_tests, function(pair) {
    data.frame(
      Comparison = paste(pair[1], "vs", pair[2]),
      P_value = nadeau_bengio_test(final_models[[pair[1]]]$resample, final_models[[pair[2]]]$resample),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(nadeau_table, out_file("Table_Nadeau_Bengio_Corrected_Tests.csv"), row.names = FALSE)

# ---------------------------
# 7. FIGURE 7: CONFUSION MATRICES AND ROC CURVES
# ---------------------------

base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(family = "Times New Roman")
  )

draw_cm <- function(cm_table, title) {
  ggplot(as.data.frame(cm_table), aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Freq), size = 7, color = "black") +
    scale_fill_gradient(low = "#EBF5FB", high = "#5DADE2") +
    labs(title = title, x = "Actual", y = "Predicted") +
    coord_fixed(ratio = 1) +
    guides(fill = "none") +
    base_theme +
    theme(axis.text = element_text(size = 12, face = "bold", color = "black"))
}

roc_objects <- lapply(agg_preds_list, function(df) {
  pROC::roc(
    df$obs,
    df$Mean_Probability_FVPTC,
    levels = c("FA", "FVPTC"),
    direction = "<",
    quiet = TRUE
  )
})

draw_single_roc <- function(roc_obj, model_name) {
  df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  ggplot(df, aes(x = fpr, y = tpr)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", color = "gray50") +
    coord_fixed(ratio = 1) +
    labs(
      title = model_name,
      subtitle = paste0("AUC = ", round(pROC::auc(roc_obj), 3)),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    base_theme +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12, face = "bold"))
}

cm_plots <- lapply(names(agg_preds_list), function(model_name) {
  draw_cm(make_cm(agg_preds_list[[model_name]]$pred, agg_preds_list[[model_name]]$obs), model_name)
})

roc_individual_plots <- lapply(names(roc_objects), function(model_name) {
  draw_single_roc(roc_objects[[model_name]], model_name)
})

roc_df <- dplyr::bind_rows(
  lapply(names(roc_objects), function(model_name) {
    tibble::tibble(
      Model = model_name,
      fpr = 1 - roc_objects[[model_name]]$specificities,
      tpr = roc_objects[[model_name]]$sensitivities
    )
  })
)

p_roc_overlay <- ggplot(roc_df, aes(x = fpr, y = tpr, color = Model)) +
  geom_line(linewidth = 1.15) +
  geom_abline(linetype = "dashed", color = "gray50") +
  coord_fixed(ratio = 1) +
  labs(
    title = "ROC Curves (Overlay)",
    subtitle = paste(
      paste0(names(roc_objects), " AUC = ", round(sapply(roc_objects, pROC::auc), 3)),
      collapse = " | "
    ),
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"))

fig7_final <- (
  wrap_plots(cm_plots, nrow = 1) /
    (wrap_plots(roc_individual_plots, nrow = 1) | p_roc_overlay) +
    plot_layout(heights = c(1, 1.15), widths = c(1.6, 1))
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20, family = "Times New Roman"))

save_figure_outputs(
  plot = fig7_final,
  filename_base = "Figure_7_CM_ROC_Top3_Final_Models",
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "NUMERIC REPORT FOR FIGURE 7",
    "===========================",
    "Final models:",
    paste("-", final_model_names),
    "",
    "Aggregated performance table:",
    capture.output(print(performance_table)),
    "",
    "Confusion matrices:",
    capture.output(print(confusion_long)),
    "",
    "Nadeau-Bengio corrected pairwise ROC tests:",
    capture.output(print(nadeau_table))
  ),
  con = out_file("Report_Figure_7_CM_ROC.txt")
)

# ---------------------------
# 8. PERMUTATION IMPORTANCE WITH SHARED PERMUTATIONS
# ---------------------------

cat("\n=======================================================\n")
cat("--- PERMUTATION FEATURE IMPORTANCE ---\n")
cat("=======================================================\n")

n_runs_perm <- 20

set.seed(SEED_GLOBAL + 3000)
perm_indices <- lapply(seq_len(n_runs_perm), function(r) {
  lapply(gene_list, function(g) sample(seq_len(nrow(gene_data))))
})
names(perm_indices) <- paste0("Perm", seq_len(n_runs_perm))
for (r in seq_len(n_runs_perm)) names(perm_indices[[r]]) <- gene_list

calc_permutation_importance <- function(model_object, model_name) {
  
  base_auc <- pROC::auc(
    pROC::roc(
      gene_data$label,
      predict(model_object, newdata = gene_data, type = "prob")[, "FVPTC"],
      levels = c("FA", "FVPTC"),
      direction = "<",
      quiet = TRUE
    )
  )
  
  perm_list <- list()
  
  for (r in seq_len(n_runs_perm)) {
    cat("Model:", model_name, "| permutation", r, "of", n_runs_perm, "\n")
    imp <- data.frame(Gene = gene_list, delta_auc = NA_real_)
    
    for (g in gene_list) {
      data_perm <- gene_data
      data_perm[, g] <- data_perm[perm_indices[[r]][[g]], g]
      
      prob_perm <- predict(model_object, newdata = data_perm, type = "prob")[, "FVPTC"]
      auc_perm <- pROC::auc(
        pROC::roc(
          data_perm$label,
          prob_perm,
          levels = c("FA", "FVPTC"),
          direction = "<",
          quiet = TRUE
        )
      )
      
      imp$delta_auc[imp$Gene == g] <- as.numeric(base_auc - auc_perm)
    }
    
    imp$Permutation_run <- r
    imp$Model <- model_name
    perm_list[[r]] <- imp
  }
  
  df <- do.call(rbind, perm_list)
  
  summary <- df %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Mean = mean(delta_auc, na.rm = TRUE),
      SD = sd(delta_auc, na.rm = TRUE),
      Median = median(delta_auc, na.rm = TRUE),
      Min = min(delta_auc, na.rm = TRUE),
      Max = max(delta_auc, na.rm = TRUE),
      Positive_runs = sum(delta_auc > 0, na.rm = TRUE),
      N_runs = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Mean))
  
  list(raw = df, summary = summary)
}

perm_results <- lapply(names(final_models), function(model_name) {
  calc_permutation_importance(final_models[[model_name]], model_name)
})
names(perm_results) <- names(final_models)

for (model_name in names(perm_results)) {
  write.csv(
    perm_results[[model_name]]$raw,
    out_file(paste0("Table_PI_Raw_", safe_name(model_name), ".csv")),
    row.names = FALSE
  )
  write.csv(
    perm_results[[model_name]]$summary,
    out_file(paste0("Table_PI_Summary_", safe_name(model_name), ".csv")),
    row.names = FALSE
  )
}

# ---------------------------
# 9. CONSENSUS SCORING AND PANEL SELECTION
# ---------------------------

pi_threshold <- 0

positive_pi_table <- dplyr::bind_rows(
  lapply(names(perm_results), function(model_name) {
    perm_results[[model_name]]$summary %>%
      dplyr::transmute(
        Model = model_name,
        Gene = as.character(Gene),
        Mean_delta_AUC = as.numeric(Mean),
        SD_delta_AUC = as.numeric(SD)
      ) %>%
      dplyr::filter(!is.na(Mean_delta_AUC), Mean_delta_AUC > pi_threshold) %>%
      dplyr::arrange(dplyr::desc(Mean_delta_AUC)) %>%
      dplyr::mutate(PI_rank = dplyr::row_number())
  })
)

if (nrow(positive_pi_table) == 0) {
  stop("No genes with positive permutation importance were found.")
}

positive_counts_by_model <- positive_pi_table %>%
  dplyr::count(Model, name = "Number_positive_PI_genes")

max_positive_pi <- max(positive_counts_by_model$Number_positive_PI_genes)

score_long <- positive_pi_table %>%
  dplyr::group_by(Model) %>%
  dplyr::arrange(dplyr::desc(Mean_delta_AUC), .by_group = TRUE) %>%
  dplyr::mutate(Rank_score = max_positive_pi - PI_rank + 1) %>%
  dplyr::ungroup()

score_wide <- score_long %>%
  dplyr::mutate(Model_ID = safe_name(Model)) %>%
  dplyr::select(Gene, Model_ID, Rank_score) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = Rank_score,
    values_fill = 0,
    names_glue = "{Model_ID}_score"
  )

rank_wide <- score_long %>%
  dplyr::mutate(Model_ID = safe_name(Model)) %>%
  dplyr::select(Gene, Model_ID, PI_rank) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = PI_rank,
    values_fill = NA,
    names_glue = "{Model_ID}_rank"
  )

positive_model_counts <- positive_pi_table %>%
  dplyr::distinct(Gene, Model) %>%
  dplyr::count(Gene, name = "Positive_models")

mean_auc_summary <- positive_pi_table %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Mean_positive_delta_AUC = mean(Mean_delta_AUC, na.rm = TRUE),
    .groups = "drop"
  )

consensus_table <- data.frame(Gene = gene_list) %>%
  dplyr::left_join(score_wide, by = "Gene") %>%
  dplyr::left_join(rank_wide, by = "Gene") %>%
  dplyr::left_join(positive_model_counts, by = "Gene") %>%
  dplyr::left_join(mean_auc_summary, by = "Gene")

score_cols <- grep("_score$", names(consensus_table), value = TRUE)
for (col in score_cols) consensus_table[[col]] <- tidyr::replace_na(consensus_table[[col]], 0)

consensus_table <- consensus_table %>%
  dplyr::mutate(
    Positive_models = tidyr::replace_na(Positive_models, 0),
    Mean_positive_delta_AUC = tidyr::replace_na(Mean_positive_delta_AUC, 0),
    Total_score = rowSums(dplyr::across(dplyr::all_of(score_cols)), na.rm = TRUE)
  ) %>%
  dplyr::arrange(dplyr::desc(Total_score), dplyr::desc(Positive_models), dplyr::desc(Mean_positive_delta_AUC)) %>%
  dplyr::mutate(
    Consensus_rank = dplyr::row_number(),
    Score_drop = Total_score - dplyr::lead(Total_score)
  ) %>%
  dplyr::select(Consensus_rank, Gene, Total_score, Score_drop, Positive_models, Mean_positive_delta_AUC, dplyr::everything())

wilcox_table <- data.frame(Gene = consensus_table$Gene, Wilcoxon_p = NA_real_)
for (g in wilcox_table$Gene) {
  stat_res <- ggpubr::compare_means(as.formula(paste(g, "~ label")), data = gene_data, method = "wilcox.test")
  wilcox_table$Wilcoxon_p[wilcox_table$Gene == g] <- stat_res$p
}

consensus_table <- consensus_table %>%
  dplyr::left_join(wilcox_table, by = "Gene") %>%
  dplyr::mutate(Expression_confirmed = Wilcoxon_p < 0.05)

# Dynamic elbow-supported set for reporting
min_elbow_size <- 3
eligible_drops <- consensus_table %>%
  dplyr::filter(Consensus_rank >= min_elbow_size, Total_score > 0, !is.na(Score_drop))

if (nrow(eligible_drops) > 0 && max(eligible_drops$Score_drop, na.rm = TRUE) > 0) {
  cut_rank <- eligible_drops$Consensus_rank[which.max(eligible_drops$Score_drop)]
} else {
  cut_rank <- min_elbow_size
}

elbow_supported_genes <- consensus_table %>%
  dplyr::filter(Consensus_rank <= cut_rank, Total_score > 0) %>%
  dplyr::pull(Gene)

# Primary dynamic elbow-supported panel for manuscript reporting
primary_panel <- elbow_supported_genes

secondary_candidates <- consensus_table %>%
  dplyr::filter(
    Total_score > 0,
    !(Gene %in% primary_panel)
  ) %>%
  dplyr::pull(Gene)

expression_confirmed_not_ml_supported <- consensus_table %>%
  dplyr::filter(
    Expression_confirmed == TRUE,
    Total_score == 0
  ) %>%
  dplyr::pull(Gene)

supplementary_genes <- secondary_candidates

selected_panels <- list(
  final_models = final_model_names,
  primary_dynamic_elbow_panel = primary_panel,
  secondary_candidates = secondary_candidates,
  expression_confirmed_not_ml_supported = expression_confirmed_not_ml_supported,
  dynamic_elbow_cut_rank = cut_rank
)

write.csv(positive_pi_table, out_file("Table_Positive_Permutation_Importance.csv"), row.names = FALSE)
write.csv(score_long, out_file("Table_Model_Wise_Dynamic_Borda_Scores.csv"), row.names = FALSE)
write.csv(consensus_table, out_file("Table_Consensus_Scores_With_Expression_Confirmation.csv"), row.names = FALSE)
saveRDS(selected_panels, out_file("Selected_Panels.rds"))

write.csv(
  data.frame(
    Panel = c(
      "Primary_dynamic_elbow_supported_panel",
      "Secondary_ML_supported_candidates",
      "Expression_confirmed_not_ML_supported"
    ),
    Genes = c(
      paste(primary_panel, collapse = ", "),
      paste(secondary_candidates, collapse = ", "),
      paste(expression_confirmed_not_ml_supported, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  ),
  out_file("Table_Selected_Panels.csv"),
  row.names = FALSE
)

# ---------------------------
# 10. FIGURE 8ABC: PI, HEATMAP, CONSENSUS TABLE
# ---------------------------

make_pi_plot <- function(df, title_txt, fill_col) {
  df$Gene <- factor(df$Gene, levels = rev(df$Gene))
  ggplot(df, aes(x = Gene, y = Mean)) +
    geom_col(fill = fill_col, width = 0.6) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45") +
    labs(title = title_txt, x = "Gene", y = expression("Mean " * Delta * "AUC")) +
    theme_classic(base_family = "Times New Roman") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
}

fill_cols <- c("#87CEEB", "#E020C0", "#9932CC")
pi_plots <- lapply(seq_along(names(perm_results)), function(i) {
  model_name <- names(perm_results)[i]
  make_pi_plot(perm_results[[model_name]]$summary, model_name, fill_cols[i])
})

ordered_index <- order(factor(gene_data$label, levels = c("FA", "FVPTC")))
gene_matrix_scaled <- scale(gene_data[ordered_index, gene_list, drop = FALSE])

annotation_row <- data.frame(Group = factor(gene_data$label[ordered_index], levels = c("FA", "FVPTC")))
rownames(annotation_row) <- rownames(gene_matrix_scaled)

ph <- pheatmap(
  gene_matrix_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(256),
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = annotation_row,
  annotation_colors = list(Group = c(FA = "#4E79A7", FVPTC = "#E15759")),
  gaps_row = sum(annotation_row$Group == "FA"),
  fontsize = 12,
  fontfamily = "Times New Roman",
  silent = TRUE
)

score_cols_for_plot <- grep("_score$", names(consensus_table), value = TRUE)
consensus_table_for_plot <- consensus_table %>%
  dplyr::filter(Total_score > 0) %>%
  dplyr::select(Rank = Consensus_rank, Gene, dplyr::all_of(score_cols_for_plot), Total = Total_score, `Wilcoxon p` = Wilcoxon_p) %>%
  dplyr::mutate(`Wilcoxon p` = signif(`Wilcoxon p`, 3))

consensus_table_grob <- tableGrob(
  consensus_table_for_plot,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontfamily = "Times New Roman", fontsize = 10)),
    colhead = list(
      bg_params = list(fill = "#E8E8E8", col = "black", lwd = 1),
      fg_params = list(fontfamily = "Times New Roman", fontsize = 11, fontface = "bold")
    )
  )
)

fig8_abc <- (
  wrap_plots(pi_plots, nrow = 1) /
    (wrap_elements(full = ph$gtable) | wrap_elements(consensus_table_grob)) +
    plot_layout(widths = c(1.25, 1))
) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20, family = "Times New Roman"))

save_figure_outputs(
  plot = fig8_abc,
  filename_base = "Figure_8ABC_Heatmap_PI_Consensus_Top3_Final_Models",
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "NUMERIC REPORT FOR FIGURE 8ABC",
    "==============================",
    "Permutation importance summaries:",
    capture.output(print(lapply(perm_results, function(x) x$summary))),
    "",
    "Consensus table:",
    capture.output(print(consensus_table)),
    "",
    "Selected panels:",
    capture.output(print(selected_panels))
  ),
  con = out_file("Report_Figure_8ABC_PI_Heatmap_Consensus.txt")
)

# -----------------------------------------------------
# 11. FIGURE 8D: DYNAMIC ELBOW-SUPPORTED PANEL LOLLIPOP
# -----------------------------------------------------

elbow_df <- consensus_table %>%
  dplyr::arrange(Consensus_rank) %>%
  dplyr::filter(Total_score > 0) %>%
  dplyr::mutate(
    Gene_label = paste0(Consensus_rank, ". ", Gene),
    Selected_status = dplyr::case_when(
      Gene %in% primary_panel ~
        "Primary Dynamic Elbow-Supported Panel",
      
      Total_score > 0 & !(Gene %in% primary_panel) ~
        "Secondary ML-Supported Candidate",
      
      TRUE ~
        "Additional Candidate"
    )
  )

dynamic_cutoff_x <- cut_rank

p_consensus_score <- ggplot(elbow_df, aes(x = Consensus_rank, y = Total_score)) +
  geom_line(linewidth = 0.8, color = "gray45") +
  geom_segment(
    aes(x = Consensus_rank, xend = Consensus_rank, y = 0, yend = Total_score),
    linewidth = 0.7,
    color = "gray60"
  ) +
  geom_point(aes(fill = Selected_status), size = 6, shape = 21, color = "black", stroke = 1) +
  geom_text(aes(label = Gene), vjust = -1.5, size = 4.5, family = "Times New Roman", fontface = "bold") +
  geom_vline(
    xintercept = dynamic_cutoff_x + 0.5,
    linetype = "dashed",
    linewidth = 0.8,
    color = "red"
  ) +
  annotate(
    "text",
    x = dynamic_cutoff_x + 0.6,
    y = max(elbow_df$Total_score, na.rm = TRUE),
    label = paste0("Dynamic Elbow Cutoff\nAfter Rank ", cut_rank),
    hjust = 0,
    vjust = 1,
    size = 4.5,
    family = "Times New Roman",
    color = "red"
  ) +
  scale_x_continuous(breaks = elbow_df$Consensus_rank, labels = elbow_df$Gene_label) +
  scale_fill_manual(
    values = c(
      "Primary Dynamic Elbow-Supported Panel" = "#5DADE2",
      "Secondary ML-Supported Candidate" = "#F39C12",
      "Additional Candidate" = "white"
    )
  ) +
  labs(
    title = "Consensus Leaderboard and Biomarker Selection",
    subtitle = "Dynamic elbow-supported panel defined by consensus ML score drop",
    x = "Consensus rank and gene",
    y = "Total consensus score",
    fill = "Biomarker classification"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray30"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  ) +
  expand_limits(y = max(elbow_df$Total_score, na.rm = TRUE) + 2)

save_figure_outputs(
  plot = p_consensus_score,
  filename_base = "Figure_8D_Lollipop_Dynamic_Elbow_Supported_Panel",
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "NUMERIC REPORT FOR FIGURE 8D",
    "===========================",
    "Consensus leaderboard:",
    capture.output(print(consensus_table)),
    "",
    paste0("Dynamic elbow cut rank: ", cut_rank),
    paste0("Primary dynamic elbow-supported panel: ", paste(primary_panel, collapse = ", ")),
    paste0("Secondary ML-supported candidates: ", paste(secondary_candidates, collapse = ", ")),
    paste0(
      "Expression-confirmed but not ML-supported genes: ",
      paste(expression_confirmed_not_ml_supported, collapse = ", ")
    )
  ),
  con = out_file("Report_Figure_8D_Lollipop_Dynamic_Elbow.txt")
)


# ---------------------------
# 12. FIGURE 8E AND SUPPLEMENTARY BOXPLOTS
# ---------------------------

plot_expression_boxplots <- function(genes, title_text) {
  plot_list <- list()
  
  for (g in genes) {
    p <- ggplot(gene_data, aes(x = label, y = .data[[g]], fill = label)) +
      geom_boxplot(outlier.colour = "red", alpha = 0.7, width = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.4) +
      scale_fill_manual(values = c("FA" = "#4E79A7", "FVPTC" = "#E15759")) +
      labs(title = g, x = NULL, y = "Expression Level (log2)") +
      theme_classic(base_family = "Times New Roman") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text = element_text(size = 11, color = "black")
      ) +
      ggpubr::stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("FA", "FVPTC")),
        label = "p.format",
        size = 5
      )
    
    plot_list[[g]] <- p
  }
  
  wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 18,
          family = "Times New Roman"
        )
      )
    )
}

wilcox_expression_report <- consensus_table %>%
  dplyr::select(
    Consensus_rank,
    Gene,
    Total_score,
    Wilcoxon_p,
    Expression_confirmed
  )

write.csv(
  wilcox_expression_report,
  out_file("Table_Expression_Boxplot_Wilcoxon_Report.csv"),
  row.names = FALSE
)

# Figure 8E: expression boxplots for the primary dynamic elbow-supported panel
if (length(primary_panel) > 0) {
  
  fig8e_expression <- plot_expression_boxplots(
    primary_panel,
    "Expression Distribution and P-values: Primary Dynamic Elbow-Supported Panel"
  )
  
  save_figure_outputs(
    plot = fig8e_expression,
    filename_base = "Figure_8E_Boxplots_Primary_Dynamic_Elbow_Panel",
    width = max(10, 4 * length(primary_panel)),
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  writeLines(
    c(
      "NUMERIC REPORT FOR FIGURE 8E",
      "===========================",
      paste0(
        "Primary dynamic elbow-supported panel: ",
        paste(primary_panel, collapse = ", ")
      ),
      "",
      "Wilcoxon/expression report:",
      capture.output(
        print(
          wilcox_expression_report %>%
            dplyr::filter(Gene %in% primary_panel)
        )
      )
    ),
    con = out_file("Report_Figure_8E_Dynamic_Elbow_Boxplots.txt")
  )
  
} else {
  warning("No primary dynamic elbow-supported genes were selected; Figure 8E was not generated.")
}

# Supplementary Figure S1: expression boxplots for secondary ML-supported candidates
if (length(supplementary_genes) > 0) {
  
  fig_supp_expression <- plot_expression_boxplots(
    supplementary_genes,
    "Supplementary: Expression Distribution of Secondary ML-Supported Candidates"
  )
  
  save_figure_outputs(
    plot = fig_supp_expression,
    filename_base = "Supplementary_Figure_S1_Secondary_Candidate_Boxplots",
    width = max(6, 4 * length(supplementary_genes)),
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  writeLines(
    c(
      "NUMERIC REPORT FOR SUPPLEMENTARY FIGURE S1",
      "=========================================",
      paste0(
        "Secondary genes: ",
        paste(supplementary_genes, collapse = ", ")
      ),
      "",
      "Wilcoxon/expression report:",
      capture.output(
        print(
          wilcox_expression_report %>%
            dplyr::filter(Gene %in% supplementary_genes)
        )
      )
    ),
    con = out_file("Report_Supplementary_Figure_S1_Secondary_Boxplots.txt")
  )
}

























# ==============================================================================
# MASTER REVIEW TABLE: ML CONSENSUS + UNIVARIATE EXPRESSION EVIDENCE
# ==============================================================================

cat("\n=======================================================\n")
cat("--- GENERATING MASTER BIOMARKER REVIEW TABLE ---\n")
cat("=======================================================\n")

# ---------------------------
# 1. Univariate expression + single-gene ROC report
# ---------------------------

univariate_gene_report <- lapply(gene_list, function(g) {
  
  x_fa <- gene_data[gene_data$label == "FA", g]
  x_fv <- gene_data[gene_data$label == "FVPTC", g]
  
  wilcox_res <- wilcox.test(x_fa, x_fv)
  
  median_fa <- median(x_fa, na.rm = TRUE)
  median_fv <- median(x_fv, na.rm = TRUE)
  mean_fa <- mean(x_fa, na.rm = TRUE)
  mean_fv <- mean(x_fv, na.rm = TRUE)
  
  roc_obj <- pROC::roc(
    response = gene_data$label,
    predictor = gene_data[[g]],
    levels = c("FA", "FVPTC"),
    quiet = TRUE
  )
  
  auc_raw <- as.numeric(pROC::auc(roc_obj))
  
  data.frame(
    Gene = g,
    Median_FA = median_fa,
    Median_FVPTC = median_fv,
    Median_log2_difference_FVPTC_minus_FA = median_fv - median_fa,
    Mean_FA = mean_fa,
    Mean_FVPTC = mean_fv,
    Mean_log2_difference_FVPTC_minus_FA = mean_fv - mean_fa,
    Wilcoxon_p_univariate = wilcox_res$p.value,
    Single_gene_AUC_raw = auc_raw,
    Single_gene_AUC_separation = max(auc_raw, 1 - auc_raw),
    Expression_direction = ifelse(
      median_fv > median_fa,
      "Higher in FVPTC",
      "Higher in FA"
    ),
    stringsAsFactors = FALSE
  )
})

univariate_gene_report <- do.call(rbind, univariate_gene_report)

univariate_gene_report$Wilcoxon_FDR_univariate <- p.adjust(
  univariate_gene_report$Wilcoxon_p_univariate,
  method = "BH"
)

univariate_gene_report <- univariate_gene_report %>%
  dplyr::arrange(
    dplyr::desc(Single_gene_AUC_separation),
    Wilcoxon_FDR_univariate
  )

write.csv(
  univariate_gene_report,
  out_file("Table_Univariate_Gene_Expression_and_Single_Gene_ROC.csv"),
  row.names = FALSE
)

# ---------------------------
# 2. Panel classification for each gene
# ---------------------------

panel_classification <- data.frame(
  Gene = gene_list,
  Panel_class = dplyr::case_when(
    gene_list %in% primary_panel ~
      "Primary dynamic elbow-supported candidate",
    
    gene_list %in% secondary_candidates ~
      "Secondary ML-supported candidate",
    
    gene_list %in% expression_confirmed_not_ml_supported ~
      "Expression-confirmed but not ML-supported",
    
    TRUE ~
      "Not selected by ML consensus"
  ),
  In_primary_dynamic_elbow_panel = gene_list %in% primary_panel,
  In_secondary_candidates = gene_list %in% secondary_candidates,
  Expression_confirmed_not_ML_supported =
    gene_list %in% expression_confirmed_not_ml_supported,
  stringsAsFactors = FALSE
)

# ---------------------------
# 3. Merge ML consensus + univariate evidence
# ---------------------------

master_biomarker_review <- consensus_table %>%
  dplyr::left_join(univariate_gene_report, by = "Gene") %>%
  dplyr::left_join(panel_classification, by = "Gene")

# Put the most important columns first, but keep all model-specific score/rank columns
priority_cols <- c(
  "Gene",
  "Panel_class",
  "Consensus_rank",
  "Total_score",
  "Score_drop",
  "Positive_models",
  "Mean_positive_delta_AUC",
  "Expression_confirmed",
  "Wilcoxon_p",
  "Wilcoxon_p_univariate",
  "Wilcoxon_FDR_univariate",
  "Single_gene_AUC_raw",
  "Single_gene_AUC_separation",
  "Expression_direction",
  "Median_FA",
  "Median_FVPTC",
  "Median_log2_difference_FVPTC_minus_FA",
  "Mean_FA",
  "Mean_FVPTC",
  "Mean_log2_difference_FVPTC_minus_FA",
  "In_primary_dynamic_elbow_panel",
  "In_secondary_candidates",
  "Expression_confirmed_not_ML_supported"
)

priority_cols <- priority_cols[priority_cols %in% colnames(master_biomarker_review)]
other_cols <- setdiff(colnames(master_biomarker_review), priority_cols)

master_biomarker_review <- master_biomarker_review %>%
  dplyr::select(dplyr::all_of(priority_cols), dplyr::all_of(other_cols)) %>%
  dplyr::arrange(Consensus_rank)

# ---------------------------
# 4. Add interpretation flags
# ---------------------------
master_biomarker_review <- master_biomarker_review %>%
  dplyr::mutate(
    Review_interpretation = dplyr::case_when(
      In_primary_dynamic_elbow_panel == TRUE ~
        "Primary dynamic elbow-supported ML-expression candidate",
      
      In_secondary_candidates == TRUE &
        Single_gene_AUC_separation >= 0.85 ~
        "Secondary ML-supported and strong univariate marker",
      
      In_secondary_candidates == TRUE &
        Expression_confirmed == TRUE ~
        "Secondary ML-supported expression-confirmed candidate",
      
      Expression_confirmed_not_ML_supported == TRUE ~
        "Expression-confirmed but not ML-supported",
      
      Single_gene_AUC_separation >= 0.85 &
        In_secondary_candidates == FALSE &
        In_primary_dynamic_elbow_panel == FALSE ~
        "Strong univariate marker but lower/no ML support",
      
      TRUE ~
        "Lower-priority candidate"
    )
  )





write.csv(
  master_biomarker_review,
  out_file("Table_MASTER_Biomarker_Review_ML_and_Univariate_Evidence.csv"),
  row.names = FALSE
)

cat("\nMaster biomarker review table saved:\n")
cat(out_file("Table_MASTER_Biomarker_Review_ML_and_Univariate_Evidence.csv"), "\n")

print(master_biomarker_review)


# ==============================================================================
# 12.5. GENERATE MASTER COMPREHENSIVE TEXT REPORT
# ==============================================================================

cat("\n=======================================================\n")
cat("--- GENERATING MASTER COMPREHENSIVE TEXT REPORT ---\n")
cat("=======================================================\n")

master_report_path <- out_file("MASTER_COMPREHENSIVE_REPORT.txt")

# Start redirecting console output to the text file
sink(master_report_path)

cat("================================================================================\n")
cat("                MASTER COMPREHENSIVE REPORT: THYROID ML PIPELINE                \n")
cat("================================================================================\n\n")

cat("DATE GENERATED: ", as.character(Sys.time()), "\n\n")

cat("--------------------------------------------------------------------------------\n")
cat("1. FINAL MODELS PERFORMANCE (SUPPORTING FIGURE 7)\n")
cat("--------------------------------------------------------------------------------\n")
print(performance_table)
cat("\n>>> Confusion Matrices:\n")
print(confusion_long)
cat("\n>>> Nadeau-Bengio Corrected Tests:\n")
print(nadeau_table)
cat("\n")

cat("--------------------------------------------------------------------------------\n")
cat("2. PERMUTATION IMPORTANCE SUMMARY (SUPPORTING FIGURE 8A-C)\n")
cat("--------------------------------------------------------------------------------\n")
for(model_name in names(perm_results)) {
  cat("\n--> Model: ", model_name, " PI Summary:\n")
  print(perm_results[[model_name]]$summary)
}
cat("\n")

cat("--------------------------------------------------------------------------------\n")
cat("3. CONSENSUS SCORING AND RANKING (SUPPORTING FIGURE 8D)\n")
cat("--------------------------------------------------------------------------------\n")
print(
  consensus_table %>% 
    dplyr::select(Consensus_rank, Gene, Total_score, Score_drop, Positive_models, Expression_confirmed)
)
cat("\n")

cat("--------------------------------------------------------------------------------\n")
cat("4. SELECTED BIOMARKER PANELS (SUPPORTING FIGURES 8D & 8E)\n")
cat("--------------------------------------------------------------------------------\n")
cat("Primary Dynamic Elbow-Supported Panel:     ", paste(primary_panel, collapse = ", "), "\n")
cat("Secondary ML Candidates:                   ", paste(secondary_candidates, collapse = ", "), "\n")
cat("Expression-Confirmed but not ML-Supported: ", paste(expression_confirmed_not_ml_supported, collapse = ", "), "\n\n")

cat("--------------------------------------------------------------------------------\n")
cat("5. EXPRESSION VALIDATION (SUPPORTING FIGURES 8E & S1)\n")
cat("--------------------------------------------------------------------------------\n")
print(wilcox_expression_report)
cat("\n")

cat("--------------------------------------------------------------------------------\n")
cat("6. MASTER BIOMARKER INTERPRETATION\n")
cat("--------------------------------------------------------------------------------\n")
print(
  master_biomarker_review %>% 
    dplyr::select(Gene, Panel_class, Consensus_rank, Single_gene_AUC_raw, Review_interpretation)
)
cat("\n")

cat("================================================================================\n")
cat("                              END OF REPORT                                     \n")
cat("================================================================================\n")

# Stop redirecting and close the file
sink()

cat("Master report saved to:\n", master_report_path, "\n")







# ---------------------------
# 13. FINAL SUMMARY AND SESSION INFO
# ---------------------------

capture.output(sessionInfo(), file = out_file("SessionInfo_Stage2_Final_Modeling.txt"))

writeLines(
  c(
    "FINAL SUMMARY",
    "=============",
    paste0("Final models: ", paste(final_model_names, collapse = ", ")),
    paste0("Secondary ML-supported candidates: ", paste(secondary_candidates, collapse = ", ")),
    paste0("Primary dynamic elbow-supported panel: ", paste(primary_panel, collapse = ", ")),
    paste0("Expression-confirmed but not ML-supported: ", paste(expression_confirmed_not_ml_supported, collapse = ", ")),
    "",
    "Key output files:",
    "- Table_Final_Model_Performance_Aggregated.csv",
    "- Table_Nadeau_Bengio_Corrected_Tests.csv",
    "- Table_Consensus_Scores_With_Expression_Confirmation.csv",
    "- Table_Selected_Panels.csv",
    "- Figure_7_CM_ROC_Top3_Final_Models.png/.tiff",
    "- Figure_8ABC_Heatmap_PI_Consensus_Top3_Final_Models.png/.tiff",
    "- Figure_8D_Lollipop_Dynamic_Elbow_Supported_Panel.png/.tiff",
    "- Figure_8E_Boxplots_Dynamic_Elbow_Supported_Panel.png/.tiff",
    "- Supplementary_Figure_S1_Secondary_Candidate_Boxplots.png/.tiff"
  ),
  con = out_file("Report_Final_Summary.txt")
)

cat("\n=======================================================\n")
cat(">>> STAGE 2 FINAL MODELING COMPLETED SUCCESSFULLY <<<\n")
cat("=======================================================\n")
cat("Final models:", paste(final_model_names, collapse = ", "), "\n")
cat("Secondary ML-supported candidates:", paste(secondary_candidates, collapse = ", "), "\n")
cat("Primary dynamic elbow-supported panel:", paste(primary_panel, collapse = ", "), "\n")
cat("Expression-confirmed but not ML-supported:", paste(expression_confirmed_not_ml_supported, collapse = ", "), "\n")
cat("Outputs saved in:\n", outdir, "\n")
cat("=======================================================\n")

# ==============================================================================
# 14. ONE-FILE COMPLETE REVIEW REPORT: XLSX + CSV + TEXT
# ==============================================================================

cat("\n=======================================================\n")
cat("--- GENERATING ONE-FILE COMPLETE REVIEW REPORT ---\n")
cat("=======================================================\n")

# This block keeps the modeling structure unchanged and only packages all
# generated numeric evidence into one coherent review file.

# Try to make an Excel workbook. If openxlsx is not installed and cannot be
# installed, the script still saves the master CSV and text report.
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("Optional package 'openxlsx' is not installed. Trying to install it from CRAN...")
  try(
    install.packages("openxlsx", repos = "https://cloud.r-project.org"),
    silent = TRUE
  )
}

make_sheet_name <- function(x, used_names = character()) {
  x <- gsub("[^A-Za-z0-9_ ]+", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- substr(x, 1, 31)
  if (nchar(x) == 0) x <- "Sheet"
  
  base <- x
  counter <- 1
  while (x %in% used_names) {
    suffix <- paste0("_", counter)
    x <- substr(base, 1, 31 - nchar(suffix))
    x <- paste0(x, suffix)
    counter <- counter + 1
  }
  x
}

to_data_frame <- function(x) {
  if (is.null(x)) {
    return(data.frame(Message = "Object not available", stringsAsFactors = FALSE))
  }
  
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  
  if (is.matrix(x) || is.table(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  
  if (is.vector(x) && !is.list(x)) {
    return(data.frame(Value = as.character(x), stringsAsFactors = FALSE))
  }
  
  data.frame(Text = capture.output(print(x)), stringsAsFactors = FALSE)
}

add_report_sheet <- function(sheet_list, sheet_name, object_value) {
  safe_sheet <- make_sheet_name(sheet_name, names(sheet_list))
  sheet_list[[safe_sheet]] <- to_data_frame(object_value)
  sheet_list
}

# ---------------------------
# Executive and interpretation sheets
# ---------------------------

executive_summary <- data.frame(
  Field = c(
    "Date_generated",
    "Project_directory",
    "Stage1_output_directory",
    "Stage2_output_directory",
    "Random_seed",
    "RNG_kind",
    "Total_samples",
    "FA_samples",
    "FVPTC_samples",
    "Genes_analyzed",
    "Final_models",
    "Primary_dynamic_elbow_supported_panel",
    "Secondary_ML_supported_candidates",
    "Dynamic_elbow_supported_genes",
    "Dynamic_elbow_cut_rank"
  ),
  Value = c(
    as.character(Sys.time()),
    project_dir,
    stage1_dir,
    outdir,
    as.character(SEED_GLOBAL),
    paste(RNGkind(), collapse = " | "),
    as.character(nrow(gene_data)),
    as.character(sum(gene_data$label == "FA")),
    as.character(sum(gene_data$label == "FVPTC")),
    paste(gene_list, collapse = ", "),
    paste(final_model_names, collapse = ", "),
    paste(primary_panel, collapse = ", "),
    paste(secondary_candidates, collapse = ", "),
    paste(elbow_supported_genes, collapse = ", "),
    as.character(cut_rank)
  ),
  stringsAsFactors = FALSE
)

interpretation_guide <- data.frame(
  Item = c(
    "MeanROC / AUC",
    "MCC",
    "Aggregated prediction",
    "Permutation importance Mean_delta_AUC",
    "Positive_models",
    "Total_score",
    "Wilcoxon_p / FDR",
    "Single_gene_AUC_separation",
    "Primary dynamic elbow-supported panel",
    "Secondary ML-supported candidates",
    "APOE interpretation rule"
  ),
  Interpretation = c(
    "Threshold-independent discrimination of the model or single gene.",
    "Balanced classification metric; useful when classes are not perfectly balanced.",
    "Repeated cross-validated FVPTC probabilities are averaged per sample before applying the fixed 0.50 threshold.",
    "Average decrease in AUC after permuting a gene; higher positive values indicate greater incremental multivariable model contribution.",
    "Number of final models in which the gene had positive permutation importance.",
    "Dynamic unweighted Borda-style consensus score across final models.",
    "Univariate expression difference between FA and FVPTC; FDR is BH-adjusted across the 8 genes.",
    "AUC corrected for direction, so values near 1 indicate stronger single-gene separation regardless of higher-in-FA or higher-in-FVPTC direction.",
    "Primary candidate panel selected by the dynamic elbow rule applied to cross-model consensus ML scores.",
    "Genes with ML support that are not in the primary dynamic elbow-supported panel.",
    "If APOE has high single-gene AUC/expression difference but lower consensus ML score, report it as an expression-supported secondary candidate rather than forcing it into the primary panel."
  ),
  stringsAsFactors = FALSE
)

selected_panels_table <- data.frame(
  Panel = c(
    "Final_models",
    "Primary_dynamic_elbow_supported_panel",
    "Secondary_ML_supported_candidates",
    "Expression_confirmed_not_ML_supported"
  ),
  Members = c(
    paste(final_model_names, collapse = ", "),
    paste(primary_panel, collapse = ", "),
    paste(secondary_candidates, collapse = ", "),
    paste(expression_confirmed_not_ml_supported, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

figure_file_index <- data.frame(
  Figure_or_report = c(
    "Supplementary Figure S0",
    "Figure 7",
    "Figure 8ABC",
    "Figure 8D",
    "Figure 8E",
    "Supplementary Figure S1",
    "Master comprehensive text report",
    "One-file Excel review report",
    "Master biomarker CSV"
  ),
  File = c(
    file.path(stage1_dir, "Supplementary_Figure_S0_Initial_Screening.png"),
    out_file("Figure_7_CM_ROC_Top3_Final_Models.png"),
    out_file("Figure_8ABC_Heatmap_PI_Consensus_Top3_Final_Models.png"),
    out_file("Figure_8D_Lollipop_Dynamic_Elbow_Supported_Panel.png"),
    out_file("Figure_8E_Boxplots_Primary_Dynamic_Elbow_Panel.png"),
    out_file("Supplementary_Figure_S1_Secondary_Candidate_Boxplots.png"),
    out_file("MASTER_COMPREHENSIVE_REPORT.txt"),
    out_file("ONE_FILE_All_Results_Review_Report.xlsx"),
    out_file("ONE_FILE_Master_Biomarker_Review.csv")
  ),
  Purpose = c(
    "Initial 7-model screening ROC distribution with mean line.",
    "Final top-3 confusion matrices and ROC curves.",
    "Permutation importance, expression heatmap, and consensus table.",
    "Consensus leaderboard and dynamic elbow-supported panel cutoff.",
    "Expression boxplots for the primary dynamic elbow-supported panel.",
    "Expression boxplots for secondary ML-supported candidates.",
    "Human-readable full text report.",
    "All key numeric outputs packaged as one workbook.",
    "Main single-table biomarker review for external checking."
  ),
  stringsAsFactors = FALSE
)

all_generated_files <- data.frame(
  Directory = c(rep(stage1_dir, length(list.files(stage1_dir, all.files = FALSE))),
                rep(outdir, length(list.files(outdir, all.files = FALSE)))),
  File = c(list.files(stage1_dir, all.files = FALSE),
           list.files(outdir, all.files = FALSE)),
  stringsAsFactors = FALSE
)

# ---------------------------
# Collect sheets
# ---------------------------

report_sheets <- list()

report_sheets <- add_report_sheet(report_sheets, "Executive_summary", executive_summary)
report_sheets <- add_report_sheet(report_sheets, "Interpretation_guide", interpretation_guide)
report_sheets <- add_report_sheet(report_sheets, "Selected_panels", selected_panels_table)

if (exists("summary_results")) {
  report_sheets <- add_report_sheet(report_sheets, "Stage1_screening_summary", summary_results)
}
if (exists("top3_candidates")) {
  report_sheets <- add_report_sheet(report_sheets, "Stage1_top3_candidates", top3_candidates)
}
if (exists("all_results")) {
  report_sheets <- add_report_sheet(report_sheets, "Stage1_run_level_ROC", all_results)
}
if (exists("all_resamples")) {
  report_sheets <- add_report_sheet(report_sheets, "Stage1_fold_resamples", all_resamples)
}

report_sheets <- add_report_sheet(report_sheets, "Final_selected_models", data.frame(Final_model_rank = 1:3, Model = final_model_names))
report_sheets <- add_report_sheet(report_sheets, "Best_tuning_parameters", best_tune_table)
report_sheets <- add_report_sheet(report_sheets, "Final_model_performance", performance_table)
report_sheets <- add_report_sheet(report_sheets, "Confusion_matrices", confusion_long)
report_sheets <- add_report_sheet(report_sheets, "Nadeau_Bengio_tests", nadeau_table)

report_sheets <- add_report_sheet(report_sheets, "Consensus_table", consensus_table)
report_sheets <- add_report_sheet(report_sheets, "Positive_PI_table", positive_pi_table)
report_sheets <- add_report_sheet(report_sheets, "Borda_model_scores", score_long)
report_sheets <- add_report_sheet(report_sheets, "Univariate_gene_ROC", univariate_gene_report)
report_sheets <- add_report_sheet(report_sheets, "Expression_Wilcoxon", wilcox_expression_report)
report_sheets <- add_report_sheet(report_sheets, "Master_biomarker_review", master_biomarker_review)

if ("APOE" %in% master_biomarker_review$Gene) {
  report_sheets <- add_report_sheet(
    report_sheets,
    "APOE_quick_check",
    master_biomarker_review %>% dplyr::filter(Gene == "APOE")
  )
}

for (model_name in names(perm_results)) {
  model_id <- safe_name(model_name)
  report_sheets <- add_report_sheet(
    report_sheets,
    paste0("PI_summary_", model_id),
    perm_results[[model_name]]$summary
  )
  report_sheets <- add_report_sheet(
    report_sheets,
    paste0("PI_raw_", model_id),
    perm_results[[model_name]]$raw
  )
}

for (model_name in names(agg_preds_list)) {
  model_id <- safe_name(model_name)
  report_sheets <- add_report_sheet(
    report_sheets,
    paste0("Aggregated_preds_", model_id),
    agg_preds_list[[model_name]]
  )
}

report_sheets <- add_report_sheet(report_sheets, "Figure_file_index", figure_file_index)
report_sheets <- add_report_sheet(report_sheets, "All_generated_files", all_generated_files)
report_sheets <- add_report_sheet(report_sheets, "Session_info", data.frame(Line = capture.output(sessionInfo()), stringsAsFactors = FALSE))

# ---------------------------
# Save the key one-file reports
# ---------------------------

write.csv(
  master_biomarker_review,
  out_file("ONE_FILE_Master_Biomarker_Review.csv"),
  row.names = FALSE
)

one_file_text_report <- out_file("ONE_FILE_Executive_Summary.txt")
writeLines(
  c(
    "ONE-FILE EXECUTIVE SUMMARY",
    "==========================",
    paste0("Date generated: ", as.character(Sys.time())),
    paste0("Final models: ", paste(final_model_names, collapse = ", ")),
    paste0("Secondary ML-supported candidates: ", paste(secondary_candidates, collapse = ", ")),
    paste0("Primary dynamic elbow-supported panel: ", paste(primary_panel, collapse = ", ")),
    paste0("Expression-confirmed but not ML-supported: ", paste(expression_confirmed_not_ml_supported, collapse = ", ")),
    "",
    "Final performance:",
    capture.output(print(performance_table)),
    "",
    "Master biomarker review:",
    capture.output(print(master_biomarker_review)),
    "",
    "Key files:",
    paste0("- ", out_file("ONE_FILE_All_Results_Review_Report.xlsx")),
    paste0("- ", out_file("ONE_FILE_Master_Biomarker_Review.csv")),
    paste0("- ", out_file("MASTER_COMPREHENSIVE_REPORT.txt"))
  ),
  con = one_file_text_report
)

xlsx_path <- out_file("ONE_FILE_All_Results_Review_Report.xlsx")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  
  for (sheet_name in names(report_sheets)) {
    openxlsx::addWorksheet(wb, sheetName = sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = report_sheets[[sheet_name]])
    
    if (ncol(report_sheets[[sheet_name]]) > 0) {
      openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = seq_len(ncol(report_sheets[[sheet_name]])),
        widths = "auto"
      )
    }
    
    openxlsx::freezePane(wb, sheet = sheet_name, firstActiveRow = 2)
  }
  
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("One-file Excel review report saved to:\n", xlsx_path, "\n")
} else {
  warning(
    "Package 'openxlsx' is not available. Excel workbook was not created. ",
    "Install it using install.packages('openxlsx') and rerun this script if XLSX is required. ",
    "The master CSV and text reports were still created."
  )
}

cat("\n=======================================================\n")
cat(">>> ONE-FILE REVIEW REPORT GENERATION COMPLETED <<<\n")
cat("=======================================================\n")
cat("Main Excel report:\n", xlsx_path, "\n")
cat("Master biomarker CSV:\n", out_file("ONE_FILE_Master_Biomarker_Review.csv"), "\n")
cat("Executive text summary:\n", one_file_text_report, "\n")
cat("All outputs directory:\n", outdir, "\n")
cat("=======================================================\n")


















