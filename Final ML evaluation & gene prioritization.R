# ==============================================================================
# 3B_METHOD_CORRECTIONS_PROFESSOR.R
# FA versus source-annotated FVPTC
#
# Purpose
# -------
# This script implements the methodological corrections requested by the
# research-methodology review, without overwriting the original 3.R outputs.
#
# Corrections implemented:
# 1. Nested repeated cross-validation for the final three classifiers.
# 2. Permutation importance calculated only on outer held-out samples.
# 3. Twenty permutations per gene, with dispersion and positive-frequency output.
# 4. Rank-stability summaries across outer resamples.
# 5. Model-normalized Borda scores so that every classifier contributes equal total weight.
# 6. Explicit zero score for mean Delta-AUC <= 0 within a classifier.
# 7. A 70% positive-resample stability filter reported only as a sensitivity analysis.
# 8. Spearman correlation analysis for all eight genes.
# 9. Post hoc panel-size sensitivity analysis for top 1 through top 8 genes.
#
# Important interpretation
# ------------------------
# The eight-gene feature set and the three model families were identified using
# the full discovery cohort in the preceding pipeline. Therefore, even this
# corrected analysis remains an internal exploratory prioritization analysis,
# not external validation.
#
# This script does NOT clear the R workspace and does NOT overwrite the original
# stage2_top3_final_reproducible folder.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Reproducibility
# ------------------------------------------------------------------------------

SEED_GLOBAL <- 123

RNGkind(
  "Mersenne-Twister",
  "Inversion",
  "Rejection"
)

set.seed(SEED_GLOBAL)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()
}

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "caret",
  "pROC",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "xgboost",
  "kernlab",
  "naivebayes"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running this script:\n",
    paste(missing_packages, collapse = ", ")
  )
}

library(caret)
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(xgboost)
library(kernlab)
library(naivebayes)

# ------------------------------------------------------------------------------
# 2. Paths and analysis settings
# ------------------------------------------------------------------------------

project_dir <- getwd()
original_stage2_dir <- file.path(
  project_dir,
  "ml-output",
  "stage2_top3_final_reproducible"
)

corrected_outdir <- file.path(
  project_dir,
  "ml-output",
  "stage2_methodologically_corrected"
)

dir.create(
  corrected_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(corrected_outdir)) {
  stop(
    "Could not create the corrected-analysis output directory:\n",
    corrected_outdir
  )
}

out_file <- function(filename) {
  file.path(
    corrected_outdir,
    filename
  )
}

safe_name <- function(x) {
  gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )
}

journal_font <- "Times New Roman"

outer_folds <- 5
outer_repeats <- 20
inner_folds <- 4
n_permutations <- 20

# This threshold is used only for a post hoc sensitivity analysis.
# It is NOT used as the sole primary selection rule.
stability_threshold <- 0.70

# Panel-size analysis is computationally intensive because it refits models.
run_panel_size_sensitivity <- TRUE

# ------------------------------------------------------------------------------
# 3. Read the processed expression matrix and exact labels
# ------------------------------------------------------------------------------

expression_file <- file.path(
  project_dir,
  "final_data_clean.rds"
)

label_file <- file.path(
  original_stage2_dir,
  "Report_Sample_Label_Order_Used_By_ML.csv"
)

selected_model_file <- file.path(
  original_stage2_dir,
  "Table_Final_Selected_Models_From_Stage1.csv"
)

required_input_files <- c(
  expression_file,
  label_file,
  selected_model_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Required input files are missing:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    )
  )
}

final_data <- readRDS(
  expression_file
)

sample_labels <- read.csv(
  label_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  arrange(
    Sample_order
  )

if (!identical(
  as.character(
    sample_labels$Sample_name
  ),
  colnames(
    final_data
  )
)) {
  stop(
    "Sample order mismatch between final_data_clean.rds and ",
    "Report_Sample_Label_Order_Used_By_ML.csv."
  )
}

sample_labels$Label <- factor(
  sample_labels$Label,
  levels = c(
    "FA",
    "FVPTC"
  )
)

if (any(is.na(sample_labels$Label))) {
  stop(
    "At least one sample label is not FA or FVPTC."
  )
}

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

missing_genes <- setdiff(
  gene_list,
  rownames(
    final_data
  )
)

if (length(missing_genes) > 0) {
  stop(
    "These genes are missing from final_data_clean.rds:\n",
    paste(
      missing_genes,
      collapse = ", "
    )
  )
}

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

for (gene_name in gene_list) {
  gene_data[[gene_name]] <- as.numeric(
    gene_data[[gene_name]]
  )
}

gene_data$label <- sample_labels$Label
rownames(gene_data) <- sample_labels$Sample_name

selected_models_table <- read.csv(
  selected_model_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  arrange(
    Final_model_rank
  )

final_model_names <- as.character(
  selected_models_table$Model
)

expected_models <- c(
  "XGBoost",
  "Naive Bayes",
  "KNN"
)

if (!identical(
  final_model_names,
  expected_models
)) {
  stop(
    "This corrected script expects the final models to be exactly:\n",
    paste(
      expected_models,
      collapse = ", "
    ),
    "\nFound:\n",
    paste(
      final_model_names,
      collapse = ", "
    )
  )
}

write.csv(
  data.frame(
    Total_samples = nrow(gene_data),
    FA = sum(gene_data$label == "FA"),
    FVPTC = sum(gene_data$label == "FVPTC"),
    Number_genes = length(gene_list),
    Outer_folds = outer_folds,
    Outer_repeats = outer_repeats,
    Inner_folds = inner_folds,
    Permutations_per_gene_per_outer_resample = n_permutations,
    stringsAsFactors = FALSE
  ),
  out_file(
    "Report_Corrected_Analysis_Design.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

safe_auc <- function(
    observed,
    probability
) {

  observed <- factor(
    observed,
    levels = c(
      "FA",
      "FVPTC"
    )
  )

  valid <- is.finite(
    probability
  ) & !is.na(
    observed
  )

  observed <- observed[
    valid
  ]

  probability <- probability[
    valid
  ]

  if (
    length(
      unique(
        observed
      )
    ) < 2
  ) {
    return(
      NA_real_
    )
  }

  roc_object <- pROC::roc(
    observed,
    probability,
    levels = c(
      "FA",
      "FVPTC"
    ),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(
    pROC::auc(
      roc_object
    )
  )
}

make_inner_control <- function(
    train_labels,
    seed
) {

  smallest_class <- min(
    table(
      train_labels
    )
  )

  k_inner <- min(
    inner_folds,
    as.integer(
      smallest_class
    )
  )

  if (k_inner < 2) {
    stop(
      "Too few samples in at least one class for inner cross-validation."
    )
  }

  set.seed(
    seed
  )

  inner_index <- caret::createFolds(
    train_labels,
    k = k_inner,
    returnTrain = TRUE
  )

  inner_index_out <- lapply(
    inner_index,
    function(train_index) {
      setdiff(
        seq_along(
          train_labels
        ),
        train_index
      )
    }
  )

  caret::trainControl(
    method = "cv",
    index = inner_index,
    indexOut = inner_index_out,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE,
    verboseIter = FALSE
  )
}

train_nested_model <- function(
    model_name,
    training_data,
    predictors,
    inner_control,
    seed
) {

  model_data <- training_data[
    ,
    c(
      predictors,
      "label"
    ),
    drop = FALSE
  ]

  set.seed(
    seed
  )

  if (model_name == "XGBoost") {

    tuning_grid <- expand.grid(
      nrounds = c(
        50,
        100,
        150
      ),
      max_depth = c(
        2,
        3
      ),
      eta = c(
        0.01,
        0.05,
        0.10
      ),
      gamma = 0,
      colsample_bytree = c(
        0.8,
        1.0
      ),
      min_child_weight = 1,
      subsample = c(
        0.8,
        1.0
      )
    )

    fitted_model <- caret::train(
      label ~ .,
      data = model_data,
      method = "xgbTree",
      metric = "ROC",
      trControl = inner_control,
      tuneGrid = tuning_grid,
      nthread = 1
    )

  } else if (model_name == "Naive Bayes") {

    tuning_grid <- expand.grid(
      laplace = c(
        0,
        1
      ),
      usekernel = c(
        TRUE,
        FALSE
      ),
      adjust = c(
        1,
        2
      )
    )

    fitted_model <- caret::train(
      label ~ .,
      data = model_data,
      method = "naive_bayes",
      metric = "ROC",
      trControl = inner_control,
      tuneGrid = tuning_grid
    )

  } else if (model_name == "KNN") {

    tuning_grid <- expand.grid(
      k = seq(
        3,
        15,
        by = 2
      )
    )

    fitted_model <- caret::train(
      label ~ .,
      data = model_data,
      method = "knn",
      metric = "ROC",
      trControl = inner_control,
      tuneGrid = tuning_grid,
      preProcess = c(
        "center",
        "scale"
      )
    )

  } else {

    stop(
      "Unsupported final model: ",
      model_name
    )
  }

  fitted_model
}

extract_tuning_row <- function(
    model_object,
    model_name,
    resample_name
) {

  tuning_row <- data.frame(
    Model = model_name,
    Resample = resample_name,
    stringsAsFactors = FALSE
  )

  cbind(
    tuning_row,
    model_object$bestTune
  )
}

extract_tune_grid <- function(
    tuning_row,
    model_name
) {

  if (model_name == "XGBoost") {

    tuning_row %>%
      dplyr::select(
        nrounds,
        max_depth,
        eta,
        gamma,
        colsample_bytree,
        min_child_weight,
        subsample
      )

  } else if (model_name == "Naive Bayes") {

    tuning_row %>%
      dplyr::select(
        laplace,
        usekernel,
        adjust
      )

  } else if (model_name == "KNN") {

    tuning_row %>%
      dplyr::select(
        k
      )

  } else {

    stop(
      "Unsupported model in extract_tune_grid(): ",
      model_name
    )
  }
}

fit_fixed_model <- function(
    model_name,
    training_data,
    predictors,
    tuning_grid,
    seed
) {

  model_data <- training_data[
    ,
    c(
      predictors,
      "label"
    ),
    drop = FALSE
  ]

  fixed_control <- caret::trainControl(
    method = "none",
    classProbs = TRUE,
    allowParallel = FALSE
  )

  set.seed(
    seed
  )

  if (model_name == "XGBoost") {

    caret::train(
      label ~ .,
      data = model_data,
      method = "xgbTree",
      trControl = fixed_control,
      tuneGrid = tuning_grid,
      nthread = 1
    )

  } else if (model_name == "Naive Bayes") {

    caret::train(
      label ~ .,
      data = model_data,
      method = "naive_bayes",
      trControl = fixed_control,
      tuneGrid = tuning_grid
    )

  } else if (model_name == "KNN") {

    caret::train(
      label ~ .,
      data = model_data,
      method = "knn",
      trControl = fixed_control,
      tuneGrid = tuning_grid,
      preProcess = c(
        "center",
        "scale"
      )
    )

  } else {

    stop(
      "Unsupported model in fit_fixed_model(): ",
      model_name
    )
  }
}

calculate_binary_metrics <- function(
    observed,
    probability,
    model_name
) {

  observed <- factor(
    observed,
    levels = c(
      "FA",
      "FVPTC"
    )
  )

  predicted <- factor(
    ifelse(
      probability >= 0.50,
      "FVPTC",
      "FA"
    ),
    levels = c(
      "FA",
      "FVPTC"
    )
  )

  confusion_matrix <- table(
    Predicted = predicted,
    Actual = observed
  )

  true_positive <- as.numeric(
    confusion_matrix[
      "FVPTC",
      "FVPTC"
    ]
  )

  true_negative <- as.numeric(
    confusion_matrix[
      "FA",
      "FA"
    ]
  )

  false_positive <- as.numeric(
    confusion_matrix[
      "FVPTC",
      "FA"
    ]
  )

  false_negative <- as.numeric(
    confusion_matrix[
      "FA",
      "FVPTC"
    ]
  )

  accuracy <- (
    true_positive +
      true_negative
  ) / sum(
    confusion_matrix
  )

  sensitivity <- ifelse(
    true_positive +
      false_negative ==
      0,
    NA_real_,
    true_positive /
      (
        true_positive +
          false_negative
      )
  )

  specificity <- ifelse(
    true_negative +
      false_positive ==
      0,
    NA_real_,
    true_negative /
      (
        true_negative +
          false_positive
      )
  )

  precision <- ifelse(
    true_positive +
      false_positive ==
      0,
    NA_real_,
    true_positive /
      (
        true_positive +
          false_positive
      )
  )

  negative_predictive_value <- ifelse(
    true_negative +
      false_negative ==
      0,
    NA_real_,
    true_negative /
      (
        true_negative +
          false_negative
      )
  )

  f1_score <- ifelse(
    is.na(
      precision
    ) ||
      is.na(
        sensitivity
      ) ||
      precision +
        sensitivity ==
        0,
    NA_real_,
    2 *
      precision *
      sensitivity /
      (
        precision +
          sensitivity
      )
  )

  balanced_accuracy <- mean(
    c(
      sensitivity,
      specificity
    ),
    na.rm = TRUE
  )

  mcc_denominator <- sqrt(
    (
      true_positive +
        false_positive
    ) *
      (
        true_positive +
          false_negative
      ) *
      (
        true_negative +
          false_positive
      ) *
      (
        true_negative +
          false_negative
      )
  )

  matthews_correlation_coefficient <- ifelse(
    mcc_denominator == 0,
    0,
    (
      true_positive *
        true_negative -
        false_positive *
        false_negative
    ) /
      mcc_denominator
  )

  data.frame(
    Model = model_name,
    AUC = safe_auc(
      observed,
      probability
    ),
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision_PPV = precision,
    NPV = negative_predictive_value,
    F1 = f1_score,
    Balanced_Accuracy = balanced_accuracy,
    MCC = matthews_correlation_coefficient,
    TP = true_positive,
    TN = true_negative,
    FP = false_positive,
    FN = false_negative,
    stringsAsFactors = FALSE
  )
}

nadeau_bengio_from_outer_auc <- function(
    auc_table,
    model_one,
    model_two,
    k_folds
) {

  paired_table <- auc_table %>%
    dplyr::filter(
      Model %in% c(
        model_one,
        model_two
      )
    ) %>%
    dplyr::select(
      Resample,
      Model,
      AUC
    ) %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = AUC
    )

  differences <- paired_table[[model_one]] -
    paired_table[[model_two]]

  differences <- differences[
    is.finite(
      differences
    )
  ]

  if (length(differences) < 2) {
    return(
      data.frame(
        Comparison = paste(
          model_one,
          "vs",
          model_two
        ),
        Mean_AUC_difference = NA_real_,
        Corrected_t = NA_real_,
        Corrected_p_value = NA_real_,
        N_paired_resamples = length(
          differences
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  mean_difference <- mean(
    differences
  )

  variance_difference <- stats::var(
    differences
  )

  number_resamples <- length(
    differences
  )

  correction <- (
    1 /
      number_resamples
  ) +
    (
      1 /
        (
          k_folds -
            1
        )
    )

  corrected_standard_error <- sqrt(
    variance_difference *
      correction
  )

  if (
    !is.finite(
      corrected_standard_error
    ) ||
      corrected_standard_error <=
        .Machine$double.eps
  ) {

    corrected_t <- ifelse(
      abs(
        mean_difference
      ) <=
        .Machine$double.eps,
      0,
      sign(
        mean_difference
      ) *
        Inf
    )

    corrected_p <- ifelse(
      is.finite(
        corrected_t
      ),
      1,
      0
    )

  } else {

    corrected_t <- mean_difference /
      corrected_standard_error

    corrected_p <- 2 *
      stats::pt(
        -abs(
          corrected_t
        ),
        df = number_resamples -
          1
      )
  }

  data.frame(
    Comparison = paste(
      model_one,
      "vs",
      model_two
    ),
    Mean_AUC_difference = mean_difference,
    Corrected_t = corrected_t,
    Corrected_p_value = corrected_p,
    N_paired_resamples = number_resamples,
    stringsAsFactors = FALSE
  )
}

save_plot_all_formats <- function(
    plot,
    filename_base,
    width,
    height,
    dpi = 600
) {

  output_base <- file.path(
    corrected_outdir,
    filename_base
  )

  ggsave(
    filename = paste0(
      output_base,
      ".png"
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  ggsave(
    filename = paste0(
      output_base,
      ".tiff"
    ),
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )

  pdf_saved <- tryCatch(
    {
      ggsave(
        filename = paste0(
          output_base,
          ".pdf"
        ),
        plot = plot,
        device = grDevices::cairo_pdf,
        width = width,
        height = height,
        units = "in",
        bg = "white"
      )
      TRUE
    },
    error = function(error_object) {
      warning(
        "PDF export failed for ",
        filename_base,
        ": ",
        conditionMessage(
          error_object
        )
      )
      FALSE
    }
  )

  invisible(
    pdf_saved
  )
}

# ------------------------------------------------------------------------------
# 5. Shared outer repeated cross-validation folds
# ------------------------------------------------------------------------------

set.seed(
  SEED_GLOBAL +
    5000
)

outer_index <- caret::createMultiFolds(
  gene_data$label,
  k = outer_folds,
  times = outer_repeats
)

outer_index_out <- lapply(
  outer_index,
  function(training_index) {
    setdiff(
      seq_len(
        nrow(
          gene_data
        )
      ),
      training_index
    )
  }
)

outer_predictions_list <- list()
outer_auc_list <- list()
outer_tuning_list <- list()
permutation_raw_list <- list()

prediction_counter <- 1
auc_counter <- 1
tuning_counter <- 1
permutation_counter <- 1

cat(
  "\nStarting nested repeated cross-validation and held-out permutation importance.\n"
)

for (
  outer_number in seq_along(
    outer_index
  )
) {

  resample_name <- names(
    outer_index
  )[
    outer_number
  ]

  cat(
    "\nOuter resample ",
    outer_number,
    " of ",
    length(
      outer_index
    ),
    ": ",
    resample_name,
    "\n",
    sep = ""
  )

  training_index <- outer_index[[
    outer_number
  ]]

  testing_index <- outer_index_out[[
    outer_number
  ]]

  outer_training_data <- gene_data[
    training_index,
    ,
    drop = FALSE
  ]

  outer_testing_data <- gene_data[
    testing_index,
    ,
    drop = FALSE
  ]

  # The same inner folds are used for all three model families in this
  # outer resample.
  inner_control <- make_inner_control(
    train_labels = outer_training_data$label,
    seed = SEED_GLOBAL +
      10000 +
      outer_number
  )

  for (
    model_number in seq_along(
      final_model_names
    )
  ) {

    model_name <- final_model_names[
      model_number
    ]

    cat(
      "  Model: ",
      model_name,
      "\n",
      sep = ""
    )

    nested_model <- train_nested_model(
      model_name = model_name,
      training_data = outer_training_data,
      predictors = gene_list,
      inner_control = inner_control,
      seed = SEED_GLOBAL +
        20000 +
        outer_number *
        100 +
        model_number
    )

    outer_tuning_list[[
      tuning_counter
    ]] <- extract_tuning_row(
      model_object = nested_model,
      model_name = model_name,
      resample_name = resample_name
    )

    tuning_counter <- tuning_counter +
      1

    base_probability <- predict(
      nested_model,
      newdata = outer_testing_data[
        ,
        gene_list,
        drop = FALSE
      ],
      type = "prob"
    )[
      ,
      "FVPTC"
    ]

    base_auc <- safe_auc(
      outer_testing_data$label,
      base_probability
    )

    outer_predictions_list[[
      prediction_counter
    ]] <- data.frame(
      Model = model_name,
      Resample = resample_name,
      rowIndex = testing_index,
      Sample = rownames(
        gene_data
      )[
        testing_index
      ],
      obs = outer_testing_data$label,
      Probability_FVPTC = as.numeric(
        base_probability
      ),
      stringsAsFactors = FALSE
    )

    prediction_counter <- prediction_counter +
      1

    outer_auc_list[[
      auc_counter
    ]] <- data.frame(
      Model = model_name,
      Resample = resample_name,
      AUC = base_auc,
      N_test = length(
        testing_index
      ),
      stringsAsFactors = FALSE
    )

    auc_counter <- auc_counter +
      1

    for (
      permutation_number in seq_len(
        n_permutations
      )
    ) {

      for (
        gene_number in seq_along(
          gene_list
        )
      ) {

        gene_name <- gene_list[
          gene_number
        ]

        permuted_test_data <- outer_testing_data[
          ,
          gene_list,
          drop = FALSE
        ]

        # This deterministic seed produces the same permutation order for
        # every classifier within the same outer resample, gene, and repeat.
        set.seed(
          SEED_GLOBAL +
            100000 +
            outer_number *
            10000 +
            permutation_number *
            100 +
            gene_number
        )

        permuted_test_data[[
          gene_name
        ]] <- sample(
          permuted_test_data[[
            gene_name
          ]],
          size = nrow(
            permuted_test_data
          ),
          replace = FALSE
        )

        permuted_probability <- predict(
          nested_model,
          newdata = permuted_test_data,
          type = "prob"
        )[
          ,
          "FVPTC"
        ]

        permuted_auc <- safe_auc(
          outer_testing_data$label,
          permuted_probability
        )

        permutation_raw_list[[
          permutation_counter
        ]] <- data.frame(
          Model = model_name,
          Resample = resample_name,
          Permutation_run = permutation_number,
          Gene = gene_name,
          Base_AUC = base_auc,
          Permuted_AUC = permuted_auc,
          delta_auc = as.numeric(
            base_auc -
              permuted_auc
          ),
          stringsAsFactors = FALSE
        )

        permutation_counter <- permutation_counter +
          1
      }
    }
  }
}

outer_predictions <- dplyr::bind_rows(
  outer_predictions_list
)

outer_auc <- dplyr::bind_rows(
  outer_auc_list
)

outer_tuning <- dplyr::bind_rows(
  outer_tuning_list
)

permutation_raw <- dplyr::bind_rows(
  permutation_raw_list
)

write.csv(
  outer_predictions,
  out_file(
    "Table_NestedCV_Outer_HeldOut_Predictions_Raw.csv"
  ),
  row.names = FALSE
)

write.csv(
  outer_auc,
  out_file(
    "Table_NestedCV_Outer_Resample_AUC.csv"
  ),
  row.names = FALSE
)

write.csv(
  outer_tuning,
  out_file(
    "Table_NestedCV_Outer_Resample_Best_Tuning.csv"
  ),
  row.names = FALSE
)

write.csv(
  permutation_raw,
  out_file(
    "Table_CrossFitted_PI_Raw_All_Models.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Nested-CV aggregated predictions and performance
# ------------------------------------------------------------------------------

aggregated_predictions <- outer_predictions %>%
  dplyr::group_by(
    Model,
    rowIndex,
    Sample
  ) %>%
  dplyr::summarise(
    obs = dplyr::first(
      obs
    ),
    Mean_Probability_FVPTC = mean(
      Probability_FVPTC,
      na.rm = TRUE
    ),
    Number_held_out_predictions = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Model,
    rowIndex
  )

write.csv(
  aggregated_predictions,
  out_file(
    "Table_NestedCV_Aggregated_Predictions.csv"
  ),
  row.names = FALSE
)

nested_performance <- dplyr::bind_rows(
  lapply(
    final_model_names,
    function(model_name) {

      model_predictions <- aggregated_predictions %>%
        dplyr::filter(
          Model ==
            model_name
        )

      calculate_binary_metrics(
        observed = model_predictions$obs,
        probability = model_predictions$Mean_Probability_FVPTC,
        model_name = model_name
      )
    }
  )
)

write.csv(
  nested_performance,
  out_file(
    "Table_NestedCV_Final_Model_Performance_Aggregated.csv"
  ),
  row.names = FALSE
)

nested_pairwise_comparisons <- combn(
  final_model_names,
  2,
  simplify = FALSE
)

nested_nadeau_bengio <- dplyr::bind_rows(
  lapply(
    nested_pairwise_comparisons,
    function(model_pair) {

      nadeau_bengio_from_outer_auc(
        auc_table = outer_auc,
        model_one = model_pair[
          1
        ],
        model_two = model_pair[
          2
        ],
        k_folds = outer_folds
      )
    }
  )
)

write.csv(
  nested_nadeau_bengio,
  out_file(
    "Table_NestedCV_Nadeau_Bengio_Corrected_Tests.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Held-out permutation-importance summaries
# ------------------------------------------------------------------------------

permutation_fold_means <- permutation_raw %>%
  dplyr::group_by(
    Model,
    Resample,
    Gene
  ) %>%
  dplyr::summarise(
    Fold_mean_delta_AUC = mean(
      delta_auc,
      na.rm = TRUE
    ),
    Fold_SD_delta_AUC = stats::sd(
      delta_auc,
      na.rm = TRUE
    ),
    N_permutations = sum(
      is.finite(
        delta_auc
      )
    ),
    Positive_permutations = sum(
      delta_auc >
        0,
      na.rm = TRUE
    ),
    Positive_permutation_fraction = mean(
      delta_auc >
        0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

permutation_raw_frequency <- permutation_raw %>%
  dplyr::group_by(
    Model,
    Gene
  ) %>%
  dplyr::summarise(
    Positive_permutations = sum(
      delta_auc >
        0,
      na.rm = TRUE
    ),
    N_permutations_total = sum(
      is.finite(
        delta_auc
      )
    ),
    Positive_permutation_fraction = mean(
      delta_auc >
        0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

permutation_summary <- permutation_fold_means %>%
  dplyr::group_by(
    Model,
    Gene
  ) %>%
  dplyr::summarise(
    Mean = mean(
      Fold_mean_delta_AUC,
      na.rm = TRUE
    ),
    SD = stats::sd(
      Fold_mean_delta_AUC,
      na.rm = TRUE
    ),
    Median = stats::median(
      Fold_mean_delta_AUC,
      na.rm = TRUE
    ),
    Min = min(
      Fold_mean_delta_AUC,
      na.rm = TRUE
    ),
    Max = max(
      Fold_mean_delta_AUC,
      na.rm = TRUE
    ),
    Positive_runs = sum(
      Fold_mean_delta_AUC >
        0,
      na.rm = TRUE
    ),
    N_runs = sum(
      is.finite(
        Fold_mean_delta_AUC
      )
    ),
    Positive_fraction = mean(
      Fold_mean_delta_AUC >
        0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    permutation_raw_frequency,
    by = c(
      "Model",
      "Gene"
    )
  ) %>%
  dplyr::arrange(
    Model,
    dplyr::desc(
      Mean
    )
  )

write.csv(
  permutation_fold_means,
  out_file(
    "Table_CrossFitted_PI_Outer_Resample_Means.csv"
  ),
  row.names = FALSE
)

write.csv(
  permutation_summary,
  out_file(
    "Table_CrossFitted_PI_Summary_All_Models.csv"
  ),
  row.names = FALSE
)

for (
  model_name in final_model_names
) {

  model_raw <- permutation_raw %>%
    dplyr::filter(
      Model ==
        model_name
    )

  model_summary <- permutation_summary %>%
    dplyr::filter(
      Model ==
        model_name
    ) %>%
    dplyr::select(
      -Model
    )

  write.csv(
    model_raw,
    out_file(
      paste0(
        "Table_PI_Raw_",
        safe_name(
          model_name
        ),
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    model_summary,
    out_file(
      paste0(
        "Table_PI_Summary_",
        safe_name(
          model_name
        ),
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 8. Rank stability across outer resamples
# ------------------------------------------------------------------------------

fold_rank_table <- permutation_fold_means %>%
  dplyr::group_by(
    Model,
    Resample
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Fold_mean_delta_AUC
    ),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    Positive_in_resample = Fold_mean_delta_AUC >
      0,
    Fold_rank = ifelse(
      Positive_in_resample,
      cumsum(
        Positive_in_resample
      ),
      NA_integer_
    ),
    In_top_four = !is.na(
      Fold_rank
    ) &
      Fold_rank <=
        4
  ) %>%
  dplyr::ungroup()

rank_stability <- fold_rank_table %>%
  dplyr::group_by(
    Model,
    Gene
  ) %>%
  dplyr::summarise(
    Positive_resample_frequency = mean(
      Positive_in_resample,
      na.rm = TRUE
    ),
    Top4_resample_frequency = mean(
      In_top_four,
      na.rm = TRUE
    ),
    Median_positive_rank = ifelse(
      all(
        is.na(
          Fold_rank
        )
      ),
      NA_real_,
      stats::median(
        Fold_rank,
        na.rm = TRUE
      )
    ),
    Rank_IQR = ifelse(
      all(
        is.na(
          Fold_rank
        )
      ),
      NA_real_,
      stats::IQR(
        Fold_rank,
        na.rm = TRUE
      )
    ),
    Number_outer_resamples = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Model,
    dplyr::desc(
      Top4_resample_frequency
    ),
    dplyr::desc(
      Positive_resample_frequency
    )
  )

write.csv(
  fold_rank_table,
  out_file(
    "Table_CrossFitted_PI_Rank_By_Outer_Resample.csv"
  ),
  row.names = FALSE
)

write.csv(
  rank_stability,
  out_file(
    "Table_CrossFitted_PI_Rank_Stability.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Model-normalized Borda scoring
# ------------------------------------------------------------------------------

model_normalized_scores <- permutation_summary %>%
  dplyr::group_by(
    Model
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Mean
    ),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    Eligible_mean_positive = is.finite(
      Mean
    ) &
      Mean >
        0,
    PI_rank = ifelse(
      Eligible_mean_positive,
      cumsum(
        Eligible_mean_positive
      ),
      NA_integer_
    ),
    Number_eligible_genes = sum(
      Eligible_mean_positive
    ),
    Raw_Borda_score = ifelse(
      Eligible_mean_positive,
      Number_eligible_genes -
        PI_rank +
        1,
      0
    ),
    Model_total_raw_Borda = sum(
      Raw_Borda_score
    ),
    Normalized_Borda_score = ifelse(
      Model_total_raw_Borda >
        0,
      Raw_Borda_score /
        Model_total_raw_Borda,
      0
    ),
    Eligible_stability_70pct = Eligible_mean_positive &
      Positive_fraction >=
        stability_threshold
  ) %>%
  dplyr::ungroup()

# Stability-filtered scoring is reported only as a sensitivity analysis.
stability_filtered_scores <- model_normalized_scores %>%
  dplyr::group_by(
    Model
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Mean
    ),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    Stability_rank = ifelse(
      Eligible_stability_70pct,
      cumsum(
        Eligible_stability_70pct
      ),
      NA_integer_
    ),
    Number_stability_eligible = sum(
      Eligible_stability_70pct
    ),
    Stability_raw_Borda = ifelse(
      Eligible_stability_70pct,
      Number_stability_eligible -
        Stability_rank +
        1,
      0
    ),
    Stability_model_total = sum(
      Stability_raw_Borda
    ),
    Stability_normalized_Borda = ifelse(
      Stability_model_total >
        0,
      Stability_raw_Borda /
        Stability_model_total,
      0
    )
  ) %>%
  dplyr::ungroup()

write.csv(
  model_normalized_scores,
  out_file(
    "Table_Model_Normalized_Borda_Scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  stability_filtered_scores,
  out_file(
    "Table_Model_Normalized_Borda_70pct_Stability_Sensitivity.csv"
  ),
  row.names = FALSE
)

normalized_score_wide <- model_normalized_scores %>%
  dplyr::mutate(
    Model_ID = safe_name(
      Model
    )
  ) %>%
  dplyr::select(
    Gene,
    Model_ID,
    Normalized_Borda_score
  ) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = Normalized_Borda_score,
    values_fill = 0,
    names_glue = "{Model_ID}_normalized_score"
  )

mean_delta_wide <- model_normalized_scores %>%
  dplyr::mutate(
    Model_ID = safe_name(
      Model
    )
  ) %>%
  dplyr::select(
    Gene,
    Model_ID,
    Mean
  ) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = Mean,
    values_fill = 0,
    names_glue = "{Model_ID}_mean_delta_AUC"
  )

positive_fraction_wide <- model_normalized_scores %>%
  dplyr::mutate(
    Model_ID = safe_name(
      Model
    )
  ) %>%
  dplyr::select(
    Gene,
    Model_ID,
    Positive_fraction
  ) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = Positive_fraction,
    values_fill = 0,
    names_glue = "{Model_ID}_positive_fraction"
  )

normalized_score_columns <- grep(
  "_normalized_score$",
  names(
    normalized_score_wide
  ),
  value = TRUE
)

corrected_consensus <- data.frame(
  Gene = gene_list,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    normalized_score_wide,
    by = "Gene"
  ) %>%
  dplyr::left_join(
    mean_delta_wide,
    by = "Gene"
  ) %>%
  dplyr::left_join(
    positive_fraction_wide,
    by = "Gene"
  )

for (
  score_column in normalized_score_columns
) {

  corrected_consensus[[
    score_column
  ]] <- tidyr::replace_na(
    corrected_consensus[[
      score_column
    ]],
    0
  )
}

mean_delta_columns <- grep(
  "_mean_delta_AUC$",
  names(
    corrected_consensus
  ),
  value = TRUE
)

positive_fraction_columns <- grep(
  "_positive_fraction$",
  names(
    corrected_consensus
  ),
  value = TRUE
)

# Univariate expression evidence is retained as a complementary descriptive
# analysis within the same 43-sample cohort. It is not independent validation.
wilcoxon_expression_table <- dplyr::bind_rows(
  lapply(
    gene_list,
    function(gene_name) {

      fa_values <- gene_data[
        gene_data$label ==
          "FA",
        gene_name
      ]

      fvptc_values <- gene_data[
        gene_data$label ==
          "FVPTC",
        gene_name
      ]

      wilcoxon_result <- suppressWarnings(
        stats::wilcox.test(
          fvptc_values,
          fa_values,
          alternative = "two.sided",
          exact = FALSE
        )
      )

      data.frame(
        Gene = gene_name,
        Median_FA = stats::median(
          fa_values,
          na.rm = TRUE
        ),
        Median_FVPTC = stats::median(
          fvptc_values,
          na.rm = TRUE
        ),
        Median_difference_FVPTC_minus_FA = stats::median(
          fvptc_values,
          na.rm = TRUE
        ) -
          stats::median(
            fa_values,
            na.rm = TRUE
          ),
        Expression_direction = ifelse(
          stats::median(
            fvptc_values,
            na.rm = TRUE
          ) >
            stats::median(
              fa_values,
              na.rm = TRUE
            ),
          "Higher in FVPTC",
          "Higher in FA"
        ),
        Wilcoxon_p = wilcoxon_result$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    Wilcoxon_BH_FDR = stats::p.adjust(
      Wilcoxon_p,
      method = "BH"
    ),
    Expression_confirmed = Wilcoxon_BH_FDR <
      0.05
  )

write.csv(
  wilcoxon_expression_table,
  out_file(
    "Table_Corrected_Univariate_Expression_Wilcoxon.csv"
  ),
  row.names = FALSE
)

corrected_consensus <- corrected_consensus %>%
  dplyr::left_join(
    wilcoxon_expression_table,
    by = "Gene"
  )

corrected_consensus <- corrected_consensus %>%
  dplyr::mutate(
    Normalized_consensus_score = rowMeans(
      dplyr::across(
        dplyr::all_of(
          normalized_score_columns
        )
      ),
      na.rm = TRUE
    ),
    Total_score = Normalized_consensus_score,
    Positive_models = rowSums(
      dplyr::across(
        dplyr::all_of(
          mean_delta_columns
        ),
        ~ .x >
          0
      ),
      na.rm = TRUE
    ),
    Mean_delta_AUC_across_models = rowMeans(
      dplyr::across(
        dplyr::all_of(
          mean_delta_columns
        )
      ),
      na.rm = TRUE
    ),
    Mean_positive_frequency_across_models = rowMeans(
      dplyr::across(
        dplyr::all_of(
          positive_fraction_columns
        )
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Normalized_consensus_score
    ),
    dplyr::desc(
      Positive_models
    ),
    dplyr::desc(
      Mean_delta_AUC_across_models
    )
  ) %>%
  dplyr::mutate(
    Consensus_rank = dplyr::row_number(),
    Score_drop = Normalized_consensus_score -
      dplyr::lead(
        Normalized_consensus_score
      )
  ) %>%
  dplyr::select(
    Consensus_rank,
    Gene,
    Normalized_consensus_score,
    Total_score,
    Score_drop,
    Positive_models,
    Mean_delta_AUC_across_models,
    Mean_positive_frequency_across_models,
    dplyr::everything()
  )

# Dynamic elbow is retained as a transparent heuristic, not as a gold standard.
minimum_elbow_size <- 3

eligible_score_drops <- corrected_consensus %>%
  dplyr::filter(
    Consensus_rank >=
      minimum_elbow_size,
    Normalized_consensus_score >
      0,
    !is.na(
      Score_drop
    )
  )

if (
  nrow(
    eligible_score_drops
  ) >
    0 &&
    max(
      eligible_score_drops$Score_drop,
      na.rm = TRUE
    ) >
    0
) {

  corrected_cut_rank <- eligible_score_drops$Consensus_rank[
    which.max(
      eligible_score_drops$Score_drop
    )
  ]

} else {

  corrected_cut_rank <- minimum_elbow_size
}

corrected_primary_panel <- corrected_consensus %>%
  dplyr::filter(
    Consensus_rank <=
      corrected_cut_rank,
    Normalized_consensus_score >
      0
  ) %>%
  dplyr::pull(
    Gene
  )

corrected_secondary_candidates <- corrected_consensus %>%
  dplyr::filter(
    Normalized_consensus_score >
      0,
    !Gene %in%
      corrected_primary_panel
  ) %>%
  dplyr::pull(
    Gene
  )

corrected_zero_score_genes <- corrected_consensus %>%
  dplyr::filter(
    Normalized_consensus_score ==
      0
  ) %>%
  dplyr::pull(
    Gene
  )

corrected_expression_confirmed_not_ml_supported <- corrected_consensus %>%
  dplyr::filter(
    Expression_confirmed,
    Normalized_consensus_score ==
      0
  ) %>%
  dplyr::pull(
    Gene
  )

write.csv(
  corrected_consensus,
  out_file(
    "Table_Corrected_Consensus_Normalized_Borda.csv"
  ),
  row.names = FALSE
)

corrected_selected_panels <- list(
  final_models = final_model_names,
  primary_dynamic_elbow_panel = corrected_primary_panel,
  secondary_candidates = corrected_secondary_candidates,
  expression_confirmed_not_ml_supported =
    corrected_expression_confirmed_not_ml_supported,
  zero_normalized_consensus_genes = corrected_zero_score_genes,
  dynamic_elbow_cut_rank = corrected_cut_rank,
  scoring_method = paste0(
    "Equal-model-weight normalized Borda; mean Delta-AUC <= 0 receives zero within that model; ",
    "70% positive-resample threshold is sensitivity analysis only."
  )
)

saveRDS(
  corrected_selected_panels,
  out_file(
    "Selected_Panels_Corrected_Normalized_Borda.rds"
  )
)

write.csv(
  data.frame(
    Panel = c(
      "Corrected_primary_dynamic_elbow_panel",
      "Corrected_secondary_candidates",
      "Expression_confirmed_not_ML_supported",
      "Zero_normalized_consensus_genes"
    ),
    Genes = c(
      paste(
        corrected_primary_panel,
        collapse = ", "
      ),
      paste(
        corrected_secondary_candidates,
        collapse = ", "
      ),
      paste(
        corrected_expression_confirmed_not_ml_supported,
        collapse = ", "
      ),
      paste(
        corrected_zero_score_genes,
        collapse = ", "
      )
    ),
    stringsAsFactors = FALSE
  ),
  out_file(
    "Table_Selected_Panels_Corrected.csv"
  ),
  row.names = FALSE
)

# Compare the corrected ranking with the original ranking when available.
original_consensus_file <- file.path(
  original_stage2_dir,
  "Table_Consensus_Scores_With_Expression_Confirmation.csv"
)

if (file.exists(
  original_consensus_file
)) {

  original_consensus <- read.csv(
    original_consensus_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ranking_comparison <- corrected_consensus %>%
    dplyr::select(
      Gene,
      Corrected_rank = Consensus_rank,
      Corrected_normalized_score = Normalized_consensus_score
    ) %>%
    dplyr::left_join(
      original_consensus %>%
        dplyr::select(
          Gene,
          Original_rank = Consensus_rank,
          Original_total_score = Total_score
        ),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      Rank_change_corrected_minus_original = Corrected_rank -
        Original_rank
    ) %>%
    dplyr::arrange(
      Corrected_rank
    )

  write.csv(
    ranking_comparison,
    out_file(
      "Table_Original_vs_Corrected_Consensus_Ranking.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 10. Stability-filtered consensus sensitivity analysis
# ------------------------------------------------------------------------------

stability_score_wide <- stability_filtered_scores %>%
  dplyr::mutate(
    Model_ID = safe_name(
      Model
    )
  ) %>%
  dplyr::select(
    Gene,
    Model_ID,
    Stability_normalized_Borda
  ) %>%
  tidyr::pivot_wider(
    names_from = Model_ID,
    values_from = Stability_normalized_Borda,
    values_fill = 0,
    names_glue = "{Model_ID}_stability_score"
  )

stability_score_columns <- grep(
  "_stability_score$",
  names(
    stability_score_wide
  ),
  value = TRUE
)

stability_consensus <- data.frame(
  Gene = gene_list,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    stability_score_wide,
    by = "Gene"
  )

for (
  score_column in stability_score_columns
) {

  stability_consensus[[
    score_column
  ]] <- tidyr::replace_na(
    stability_consensus[[
      score_column
    ]],
    0
  )
}

stability_consensus <- stability_consensus %>%
  dplyr::mutate(
    Stability_consensus_score = rowMeans(
      dplyr::across(
        dplyr::all_of(
          stability_score_columns
        )
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Stability_consensus_score
    )
  ) %>%
  dplyr::mutate(
    Stability_consensus_rank = dplyr::row_number()
  )

write.csv(
  stability_consensus,
  out_file(
    "Table_Stability_Filtered_70pct_Consensus_Sensitivity.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Spearman correlation analysis for all eight genes
# ------------------------------------------------------------------------------

spearman_matrix <- stats::cor(
  gene_data[
    ,
    gene_list,
    drop = FALSE
  ],
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  spearman_matrix,
  out_file(
    "Table_All8_Genes_Spearman_Correlation_Matrix.csv"
  ),
  row.names = TRUE
)

correlation_pairs <- list()
correlation_counter <- 1

for (
  first_gene_number in seq_len(
    length(
      gene_list
    ) -
      1
  )
) {

  for (
    second_gene_number in seq(
      first_gene_number +
        1,
      length(
        gene_list
      )
    )
  ) {

    first_gene <- gene_list[
      first_gene_number
    ]

    second_gene <- gene_list[
      second_gene_number
    ]

    correlation_test <- suppressWarnings(
      stats::cor.test(
        gene_data[[
          first_gene
        ]],
        gene_data[[
          second_gene
        ]],
        method = "spearman",
        exact = FALSE
      )
    )

    correlation_pairs[[
      correlation_counter
    ]] <- data.frame(
      Gene_1 = first_gene,
      Gene_2 = second_gene,
      Spearman_rho = unname(
        correlation_test$estimate
      ),
      P_value = correlation_test$p.value,
      Absolute_rho = abs(
        unname(
          correlation_test$estimate
        )
      ),
      High_correlation_abs_rho_ge_0_80 = abs(
        unname(
          correlation_test$estimate
        )
      ) >=
        0.80,
      stringsAsFactors = FALSE
    )

    correlation_counter <- correlation_counter +
      1
  }
}

correlation_pair_table <- dplyr::bind_rows(
  correlation_pairs
) %>%
  dplyr::mutate(
    BH_FDR = stats::p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Absolute_rho
    )
  )

write.csv(
  correlation_pair_table,
  out_file(
    "Table_All8_Genes_Spearman_Correlation_Pairs.csv"
  ),
  row.names = FALSE
)

correlation_plot_data <- as.data.frame(
  as.table(
    spearman_matrix
  ),
  stringsAsFactors = FALSE
)

colnames(
  correlation_plot_data
) <- c(
  "Gene_1",
  "Gene_2",
  "Spearman_rho"
)

correlation_plot_data$Gene_1 <- factor(
  correlation_plot_data$Gene_1,
  levels = gene_list
)

correlation_plot_data$Gene_2 <- factor(
  correlation_plot_data$Gene_2,
  levels = rev(
    gene_list
  )
)

correlation_plot <- ggplot(
  correlation_plot_data,
  aes(
    x = Gene_1,
    y = Gene_2,
    fill = Spearman_rho
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman_rho
      )
    ),
    family = journal_font,
    size = 3.1
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  labs(
    title = "Spearman correlation among all eight candidate genes",
    subtitle = paste0(
      "All ",
      nrow(
        gene_data
      ),
      " integrated samples; descriptive pooled-cohort analysis"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 11,
    base_family = journal_font
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

save_plot_all_formats(
  plot = correlation_plot,
  filename_base = "Supplementary_Correction_All8_Spearman_Correlation",
  width = 8.5,
  height = 7.5
)

# ------------------------------------------------------------------------------
# 12. Corrected consensus and stability plots
# ------------------------------------------------------------------------------

consensus_plot_data <- corrected_consensus %>%
  dplyr::mutate(
    Gene = factor(
      Gene,
      levels = rev(
        Gene
      )
    )
  )

corrected_consensus_plot <- ggplot(
  consensus_plot_data,
  aes(
    x = Normalized_consensus_score,
    y = Gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Normalized_consensus_score,
      yend = Gene
    ),
    linewidth = 0.7,
    color = "gray45"
  ) +
  geom_point(
    size = 3
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  labs(
    title = "Equal-model-weight normalized Borda consensus",
    subtitle = paste0(
      "Mean held-out permutation importance; dynamic-elbow cut after rank ",
      corrected_cut_rank
    ),
    x = "Normalized consensus score",
    y = NULL
  ) +
  theme_classic(
    base_size = 11,
    base_family = journal_font
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

save_plot_all_formats(
  plot = corrected_consensus_plot,
  filename_base = "Correction_Figure_Normalized_Borda_Consensus",
  width = 8.5,
  height = 5.5
)

rank_stability_plot <- rank_stability %>%
  dplyr::mutate(
    Gene = factor(
      Gene,
      levels = rev(
        gene_list
      )
    )
  ) %>%
  ggplot(
    aes(
      x = Top4_resample_frequency,
      y = Gene
    )
  ) +
  geom_point(
    size = 2.6
  ) +
  facet_wrap(
    ~ Model,
    nrow = 1
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    title = "Cross-fitted rank stability",
    subtitle = "Frequency of ranking among the top four positive-importance genes",
    x = "Top-four frequency across outer resamples",
    y = NULL
  ) +
  theme_classic(
    base_size = 10.5,
    base_family = journal_font
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

save_plot_all_formats(
  plot = rank_stability_plot,
  filename_base = "Supplementary_Correction_PI_Rank_Stability",
  width = 13,
  height = 5.5
)

# ------------------------------------------------------------------------------
# 13. Post hoc panel-size sensitivity analysis
# ------------------------------------------------------------------------------

panel_size_performance <- NULL

if (isTRUE(
  run_panel_size_sensitivity
)) {

  cat(
    "\nStarting post hoc top-k panel-size sensitivity analysis.\n"
  )

  ranked_genes <- as.character(
    corrected_consensus$Gene
  )

  panel_prediction_list <- list()
  panel_prediction_counter <- 1

  for (
    number_genes in seq_along(
      ranked_genes
    )
  ) {

    selected_genes <- ranked_genes[
      seq_len(
        number_genes
      )
    ]

    cat(
      "  Panel size ",
      number_genes,
      ": ",
      paste(
        selected_genes,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )

    for (
      outer_number in seq_along(
        outer_index
      )
    ) {

      resample_name <- names(
        outer_index
      )[
        outer_number
      ]

      training_index <- outer_index[[
        outer_number
      ]]

      testing_index <- outer_index_out[[
        outer_number
      ]]

      outer_training_data <- gene_data[
        training_index,
        ,
        drop = FALSE
      ]

      outer_testing_data <- gene_data[
        testing_index,
        ,
        drop = FALSE
      ]

      for (
        model_number in seq_along(
          final_model_names
        )
      ) {

        model_name <- final_model_names[
          model_number
        ]

        tuning_row <- outer_tuning %>%
          dplyr::filter(
            Model ==
              model_name,
            Resample ==
              resample_name
          )

        if (nrow(
          tuning_row
        ) !=
          1) {
          stop(
            "Expected exactly one tuning row for ",
            model_name,
            " / ",
            resample_name,
            "."
          )
        }

        fixed_tune_grid <- extract_tune_grid(
          tuning_row = tuning_row,
          model_name = model_name
        )

        fixed_model <- fit_fixed_model(
          model_name = model_name,
          training_data = outer_training_data,
          predictors = selected_genes,
          tuning_grid = fixed_tune_grid,
          seed = SEED_GLOBAL +
            300000 +
            number_genes *
            10000 +
            outer_number *
            100 +
            model_number
        )

        panel_probability <- predict(
          fixed_model,
          newdata = outer_testing_data[
            ,
            selected_genes,
            drop = FALSE
          ],
          type = "prob"
        )[
          ,
          "FVPTC"
        ]

        panel_prediction_list[[
          panel_prediction_counter
        ]] <- data.frame(
          Number_genes = number_genes,
          Genes = paste(
            selected_genes,
            collapse = ", "
          ),
          Model = model_name,
          Resample = resample_name,
          rowIndex = testing_index,
          Sample = rownames(
            gene_data
          )[
            testing_index
          ],
          obs = outer_testing_data$label,
          Probability_FVPTC = as.numeric(
            panel_probability
          ),
          stringsAsFactors = FALSE
        )

        panel_prediction_counter <- panel_prediction_counter +
          1
      }
    }
  }

  panel_predictions_raw <- dplyr::bind_rows(
    panel_prediction_list
  )

  write.csv(
    panel_predictions_raw,
    out_file(
      "Table_Panel_Size_Sensitivity_HeldOut_Predictions_Raw.csv"
    ),
    row.names = FALSE
  )

  panel_predictions_aggregated <- panel_predictions_raw %>%
    dplyr::group_by(
      Number_genes,
      Genes,
      Model,
      rowIndex,
      Sample
    ) %>%
    dplyr::summarise(
      obs = dplyr::first(
        obs
      ),
      Mean_Probability_FVPTC = mean(
        Probability_FVPTC,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    panel_predictions_aggregated,
    out_file(
      "Table_Panel_Size_Sensitivity_Aggregated_Predictions.csv"
    ),
    row.names = FALSE
  )

  panel_size_performance <- panel_predictions_aggregated %>%
    dplyr::group_by(
      Number_genes,
      Genes,
      Model
    ) %>%
    dplyr::summarise(
      AUC = safe_auc(
        obs,
        Mean_Probability_FVPTC
      ),
      .groups = "drop"
    ) %>%
    dplyr::group_by(
      Model
    ) %>%
    dplyr::mutate(
      Model_maximum_AUC = max(
        AUC,
        na.rm = TRUE
      ),
      Difference_from_model_maximum = Model_maximum_AUC -
        AUC,
      Within_0_01_of_model_maximum = Difference_from_model_maximum <=
        0.01
    ) %>%
    dplyr::ungroup()

  panel_size_mean <- panel_size_performance %>%
    dplyr::group_by(
      Number_genes,
      Genes
    ) %>%
    dplyr::summarise(
      Model = "Mean across models",
      AUC = mean(
        AUC,
        na.rm = TRUE
      ),
      Model_maximum_AUC = NA_real_,
      Difference_from_model_maximum = NA_real_,
      Within_0_01_of_model_maximum = NA,
      .groups = "drop"
    )

  panel_size_performance_with_mean <- dplyr::bind_rows(
    panel_size_performance,
    panel_size_mean
  )

  write.csv(
    panel_size_performance_with_mean,
    out_file(
      "Table_Panel_Size_Sensitivity_Performance.csv"
    ),
    row.names = FALSE
  )

  panel_size_plot <- ggplot(
    panel_size_performance_with_mean,
    aes(
      x = Number_genes,
      y = AUC,
      group = Model,
      linetype = Model,
      shape = Model
    )
  ) +
    geom_line(
      linewidth = 0.75
    ) +
    geom_point(
      size = 2.4
    ) +
    scale_x_continuous(
      breaks = seq_along(
        ranked_genes
      )
    ) +
    scale_y_continuous(
      limits = c(
        0.50,
        1.00
      )
    ) +
    labs(
      title = "Post hoc panel-size sensitivity",
      subtitle = "Top-k order fixed by the corrected cross-fitted normalized consensus ranking",
      x = "Number of top-ranked genes",
      y = "Aggregated held-out ROC-AUC",
      linetype = NULL,
      shape = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family = journal_font
    ) +
    theme(
      axis.text = element_text(
        color = "black"
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.position = "bottom"
    )

  save_plot_all_formats(
    plot = panel_size_plot,
    filename_base = "Supplementary_Correction_Panel_Size_Sensitivity",
    width = 10,
    height = 6.2
  )
}

# ------------------------------------------------------------------------------
# 14. Human-readable methodological report
# ------------------------------------------------------------------------------

high_correlation_pairs <- correlation_pair_table %>%
  dplyr::filter(
    High_correlation_abs_rho_ge_0_80
  )

report_lines <- c(
  "METHODOLOGICALLY CORRECTED ML ANALYSIS",
  "=====================================",
  "",
  paste0(
    "Generated: ",
    as.character(
      Sys.time()
    )
  ),
  "",
  "DESIGN",
  "------",
  paste0(
    "Outer CV: stratified ",
    outer_folds,
    "-fold CV repeated ",
    outer_repeats,
    " times."
  ),
  paste0(
    "Inner tuning CV: ",
    inner_folds,
    " folds within each outer-training set."
  ),
  paste0(
    "Permutation repeats: ",
    n_permutations,
    " per gene per outer resample."
  ),
  "Permutation importance was evaluated only on outer held-out samples.",
  "",
  "IMPORTANT REMAINING LIMITATION",
  "------------------------------",
  paste0(
    "The eight-gene feature set and the three model families were selected using ",
    "the full discovery cohort before this corrected stage. Therefore, the analysis ",
    "remains an internal exploratory prioritization analysis and is not external validation."
  ),
  "",
  "NESTED-CV AGGREGATED PERFORMANCE",
  "--------------------------------",
  capture.output(
    print(
      nested_performance
    )
  ),
  "",
  "NADEAU-BENGIO CORRECTED COMPARISONS",
  "-----------------------------------",
  capture.output(
    print(
      nested_nadeau_bengio
    )
  ),
  "",
  "CORRECTED NORMALIZED CONSENSUS",
  "------------------------------",
  capture.output(
    print(
      corrected_consensus %>%
        dplyr::select(
          Consensus_rank,
          Gene,
          Normalized_consensus_score,
          Positive_models,
          Mean_delta_AUC_across_models,
          Mean_positive_frequency_across_models
        )
    )
  ),
  "",
  paste0(
    "Corrected dynamic-elbow cut rank: ",
    corrected_cut_rank
  ),
  paste0(
    "Corrected primary panel: ",
    paste(
      corrected_primary_panel,
      collapse = ", "
    )
  ),
  paste0(
    "Corrected secondary candidates: ",
    paste(
      corrected_secondary_candidates,
      collapse = ", "
    )
  ),
  paste0(
    "Expression-confirmed but not ML-supported: ",
    paste(
      corrected_expression_confirmed_not_ml_supported,
      collapse = ", "
    )
  ),
  "",
  "HIGH CORRELATION PAIRS (ABS SPEARMAN RHO >= 0.80)",
  "-------------------------------------------------",
  if (
    nrow(
      high_correlation_pairs
    ) ==
      0
  ) {
    "None detected."
  } else {
    capture.output(
      print(
        high_correlation_pairs
      )
    )
  },
  "",
  "INTERPRETATION RULES",
  "--------------------",
  paste0(
    "- Mean held-out Delta-AUC <= 0 receives zero Borda score within that model."
  ),
  paste0(
    "- Each model's raw Borda scores are divided by their within-model sum, ",
    "so every model contributes equal total weight."
  ),
  paste0(
    "- Positive-resample frequency and top-four frequency are reported as stability measures."
  ),
  paste0(
    "- The 70% stability threshold is a sensitivity analysis, not a preregistered primary rule."
  ),
  paste0(
    "- Panel-size curves are post hoc internal sensitivity analyses and do not constitute ",
    "independent validation or an unbiased final panel-size selection procedure."
  )
)

writeLines(
  report_lines,
  con = out_file(
    "MASTER_Methodological_Corrections_Report.txt"
  )
)

capture.output(
  sessionInfo(),
  file = out_file(
    "SessionInfo_Methodologically_Corrected_Analysis.txt"
  )
)

cat(
  "\n=======================================================\n"
)

cat(
  "Corrected analysis completed.\n"
)

cat(
  "Outputs were saved in:\n",
  corrected_outdir,
  "\n",
  sep = ""
)

cat(
  "Corrected primary panel:\n",
  paste(
    corrected_primary_panel,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Do not update the manuscript or Figure 4 until the corrected output tables ",
  "have been reviewed.\n",
  sep = ""
)

cat(
  "=======================================================\n"
)
