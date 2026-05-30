# =============================================================================
# MSS project: MRS construction and analysis
# File: code/05_mrs_construction_and_analysis.R
#
# Purpose:
#   - Construct a metabolic risk score (MRS) associated with sleep health
#   - Evaluate biomarker importance using LightGBM + SHAP
#   - Estimate biomarker weights using ridge regression
#   - Calculate participant-level MRS
#   - Examine association of MRS with MSS and health outcomes
#   - Optionally perform mediation analysis
#
# Required input:
#   - data/derived/final_imp.csv
#
# Notes:
#   - This repository does not include participant-level biomarker data.
#   - Users must ensure that the required biomarker variables are available in
#     final_imp.csv or merge them from authorized source data locally.
# =============================================================================

source("code/00_setup.R")

# ---- 1. Load analysis-ready dataset ------------------------------------------

file_final_imp <- file.path(dir_derived, "final_imp.csv")
final_imp <- safe_read_csv(file_final_imp)

message("Loaded final_imp with N = ", nrow(final_imp))

if (!"mss_std" %in% names(final_imp) && "mss" %in% names(final_imp)) {
  final_imp$mss_std <- scale_numeric(final_imp$mss)
}

# ---- 2. Biomarker setup ------------------------------------------------------

# Biomarkers used in the MRS workflow.
biomarkers <- c(
  "glucose", "triglycerides", "cholesterol", "hdl_cholesterol",
  "ldl", "alanine_aminotransferase", "albumin", "creatinine",
  "cystatin_c", "urea", "calcium", "phosphate", "c_reactive_protein"
)

biomarkers_present <- intersect(biomarkers, names(final_imp))

if (length(biomarkers_present) < 5) {
  stop(
    "Too few biomarker variables found in final_imp.csv. Found: ",
    paste(biomarkers_present, collapse = ", "),
    call. = FALSE
  )
}

message("Biomarkers available: ", paste(biomarkers_present, collapse = ", "))

# ---- 3. Log-scale standardization -------------------------------------------

for (bm in biomarkers_present) {
  final_imp[[paste0(bm, "_std")]] <- log_scale(final_imp[[bm]])
}

std_biomarker_vars <- paste0(biomarkers_present, "_std")

# Optional: rename selected standardized biomarkers to concise labels
rename_map <- c(
  "glucose_std" = "Glucose",
  "triglycerides_std" = "Triglycerides",
  "cholesterol_std" = "Cholesterol",
  "hdl_cholesterol_std" = "HDL",
  "ldl_std" = "LDL",
  "alanine_aminotransferase_std" = "ALT",
  "albumin_std" = "Albumin",
  "creatinine_std" = "Creatinine",
  "cystatin_c_std" = "CYC",
  "urea_std" = "Urea",
  "calcium_std" = "Calcium",
  "phosphate_std" = "Phosphate",
  "c_reactive_protein_std" = "CRP"
)

rename_map <- rename_map[names(rename_map) %in% names(final_imp)]

if (length(rename_map) > 0) {
  final_imp <- final_imp %>%
    rename(!!!rename_map)
}

# Update standardized biomarker variable names after renaming
std_biomarker_vars <- intersect(unname(rename_map), names(final_imp))

safe_write_csv(
  final_imp,
  file.path(dir_derived, "final_imp_with_standardized_biomarkers.csv")
)

# ---- 4. Train/validation split ----------------------------------------------

set.seed(250503)

model_data <- final_imp %>%
  select(any_of(c("mss", std_biomarker_vars))) %>%
  na.omit()

if (!"mss" %in% names(model_data)) {
  stop("Variable 'mss' is required for MRS construction.", call. = FALSE)
}

n_total <- nrow(model_data)
train_index <- sample(seq_len(n_total), size = floor(0.7 * n_total))

sm_train <- model_data[train_index, , drop = FALSE]
sm_validate <- model_data[-train_index, , drop = FALSE]

X_train <- as.matrix(sm_train[, std_biomarker_vars, drop = FALSE])
y_train <- sm_train$mss

X_validate <- as.matrix(sm_validate[, std_biomarker_vars, drop = FALSE])
y_validate <- sm_validate$mss

message("Training set N = ", nrow(sm_train))
message("Validation set N = ", nrow(sm_validate))

# ---- 5. LightGBM + SHAP feature importance ----------------------------------

shap_importance_df <- NULL

if (requireNamespace("lightgbm", quietly = TRUE) &&
    requireNamespace("shapviz", quietly = TRUE)) {
  
  dtrain <- lightgbm::lgb.Dataset(data = X_train, label = y_train)
  dvalidate <- lightgbm::lgb.Dataset(data = X_validate, label = y_validate)
  
  params <- list(
    objective = "regression",
    metric = "rmse",
    boosting_type = "gbdt",
    num_leaves = 31,
    learning_rate = 0.05,
    verbose = -1
  )
  
  gb_model <- lightgbm::lgb.train(
    params = params,
    data = dtrain,
    nrounds = 300,
    valids = list(val = dvalidate),
    early_stopping_rounds = 10
  )
  
  sv <- shapviz::shapviz(gb_model, X_pred = X_train, X = X_train)
  shap_matrix <- as.matrix(shapviz::get_shap_values(sv))
  
  shap_importance_df <- data.frame(
    Feature = colnames(shap_matrix),
    Importance = colMeans(abs(shap_matrix), na.rm = TRUE),
    Direction = colMeans(shap_matrix, na.rm = TRUE),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Importance))
  
  shap_importance_df$Cumulative_Contribution <-
    cumsum(shap_importance_df$Importance) / sum(shap_importance_df$Importance)
  
  safe_write_csv(
    shap_importance_df,
    file.path(dir_tables, "mrs_shap_feature_importance.csv")
  )
  
  # validation performance
  y_pred_train <- predict(gb_model, X_train)
  y_pred_validate <- predict(gb_model, X_validate)
  
  perf_df <- data.frame(
    Dataset = c("Train", "Validation"),
    RMSE = c(
      sqrt(mean((y_pred_train - y_train)^2, na.rm = TRUE)),
      sqrt(mean((y_pred_validate - y_validate)^2, na.rm = TRUE))
    ),
    R2 = c(
      1 - sum((y_train - y_pred_train)^2) / sum((y_train - mean(y_train))^2),
      1 - sum((y_validate - y_pred_validate)^2) / sum((y_validate - mean(y_validate))^2)
    )
  )
  
  safe_write_csv(
    perf_df,
    file.path(dir_tables, "mrs_lightgbm_performance.csv")
  )
} else {
  message("Skipping LightGBM/SHAP analysis because required packages are not installed.")
}

# ---- 6. Select biomarkers for ridge regression -------------------------------

# If SHAP results are available, keep features up to 75% cumulative contribution.
# Otherwise, use all available standardized biomarkers.
if (!is.null(shap_importance_df) && nrow(shap_importance_df) > 0) {
  selected_features <- shap_importance_df %>%
    filter(Cumulative_Contribution <= 0.75) %>%
    pull(Feature)
  
  if (length(selected_features) < 3) {
    selected_features <- shap_importance_df %>%
      slice_head(n = min(5, nrow(shap_importance_df))) %>%
      pull(Feature)
  }
} else {
  selected_features <- std_biomarker_vars
}

selected_features <- intersect(selected_features, names(final_imp))

message("Selected biomarkers for ridge regression: ", paste(selected_features, collapse = ", "))

# ---- 7. Ridge regression to derive weights ----------------------------------

ridge_input <- final_imp %>%
  select(any_of(c("mss", selected_features))) %>%
  na.omit()

calculate_ridge_betas <- function(data, selected_features, outcome_var = "mss") {
  X <- as.matrix(data[, selected_features, drop = FALSE])
  y <- data[[outcome_var]]
  
  set.seed(123)
  cv_model <- glmnet::cv.glmnet(X, y, alpha = 0)
  best_lambda <- cv_model$lambda.min
  
  final_model <- glmnet::glmnet(X, y, alpha = 0, lambda = best_lambda)
  beta_mat <- as.matrix(coef(final_model))
  
  out <- data.frame(
    Feature = selected_features,
    Beta = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(out))) {
    feature_i <- out$Feature[i]
    if (feature_i %in% rownames(beta_mat)) {
      out$Beta[i] <- beta_mat[feature_i, 1]
    } else {
      out$Beta[i] <- 0
    }
  }
  
  out %>%
    arrange(desc(abs(Beta)))
}

beta_coefficients <- calculate_ridge_betas(
  data = ridge_input,
  selected_features = selected_features,
  outcome_var = "mss"
)

safe_write_csv(
  beta_coefficients,
  file.path(dir_tables, "mrs_ridge_coefficients.csv")
)

# ---- 8. Calculate participant-level MRS -------------------------------------

calculate_metabolite_score <- function(data, beta_df) {
  data$MRS <- 0
  
  for (i in seq_len(nrow(beta_df))) {
    feature_i <- beta_df$Feature[i]
    beta_i <- beta_df$Beta[i]
    
    if (feature_i %in% names(data)) {
      data$MRS <- data$MRS + data[[feature_i]] * beta_i
    }
  }
  
  data
}

final_mrs <- calculate_metabolite_score(final_imp, beta_coefficients)

final_mrs$MRS_std <- scale_numeric(final_mrs$MRS)

safe_write_csv(
  final_mrs,
  file.path(dir_derived, "final_mrs.csv")
)

# ---- 9. Correlation between MRS and MSS -------------------------------------

cor_results <- list()

if (all(c("MRS", "mss") %in% names(final_mrs))) {
  cor_train <- calculate_metabolite_score(
    final_imp[train_index, , drop = FALSE],
    beta_coefficients
  )
  
  cor_validate <- calculate_metabolite_score(
    final_imp[-train_index, , drop = FALSE],
    beta_coefficients
  )
  
  cor_results$train <- data.frame(
    Dataset = "Train",
    Correlation = tryCatch(
      cor.test(cor_train$MRS, cor_train$mss, method = "pearson")$estimate,
      error = function(e) NA_real_
    ),
    P_Value = tryCatch(
      cor.test(cor_train$MRS, cor_train$mss, method = "pearson")$p.value,
      error = function(e) NA_real_
    )
  )
  
  cor_results$validation <- data.frame(
    Dataset = "Validation",
    Correlation = tryCatch(
      cor.test(cor_validate$MRS, cor_validate$mss, method = "pearson")$estimate,
      error = function(e) NA_real_
    ),
    P_Value = tryCatch(
      cor.test(cor_validate$MRS, cor_validate$mss, method = "pearson")$p.value,
      error = function(e) NA_real_
    )
  )
  
  safe_write_csv(
    bind_rows(cor_results),
    file.path(dir_tables, "mrs_mss_correlation.csv")
  )
}

# ---- 10. MRS and disease outcomes -------------------------------------------

disease_params <- list(
  t2d = list(time = "t_t2d", event = "t2d", fh = "fh_dm"),
  obese = list(time = "t_obese", event = "obese", fh = "fh_cmd"),
  hp = list(time = "t_hp", event = "hp", fh = "fh_hp"),
  lipid = list(time = "t_lipid", event = "lipid", fh = "fh_cmd"),
  ihd = list(time = "t_ihd", event = "ihd", fh = "fh_heart"),
  af = list(time = "t_af", event = "af", fh = "fh_heart"),
  conduction = list(time = "t_conduction", event = "conduction", fh = "fh_heart"),
  hf = list(time = "t_hf", event = "hf", fh = "fh_heart"),
  is = list(time = "t_is", event = "is", fh = "fh_stroke"),
  dementia_ad = list(time = "t_dementia_ad", event = "dementia_ad", fh = "fh_ad"),
  depressive = list(time = "t_depressive", event = "depressive", fh = "fh_depression"),
  anxiety = list(time = "t_anxiety", event = "anxiety", fh = "fh_depression"),
  cancer = list(time = "t_cancer", event = "cancer", fh = "fh_cancer"),
  death = list(time = "t_death", event = "death", fh = "fh_cmd")
)

categorical_covariates <- c(
  "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "smoke", "drink", "mvpa3", "fh_var"
)

model_adjustments <- list(
  model2 = c("age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
             "ms", "smoke", "drink", "mvpa3", "fh_var"),
  model3 = c("age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
             "ms", "smoke", "drink", "mvpa3", "fh_var", "bmi", "sys")
)

build_disease_dataset <- function(data, disease_name, params) {
  out <- data %>%
    filter(!( .data[[params$event]] == 1 & .data[[params$time]] < 0 ))
  
  if (disease_name == "t2d" && "glucose" %in% names(out)) {
    out <- out %>% filter(glucose < 11 | is.na(glucose))
  }
  
  if (disease_name == "obese" && "bmi" %in% names(out)) {
    out <- out %>% filter(bmi < 30 | is.na(bmi))
  }
  
  out
}

analysis_datasets <- lapply(names(disease_params), function(dz) {
  build_disease_dataset(final_mrs, dz, disease_params[[dz]])
})
names(analysis_datasets) <- names(disease_params)

run_continuous_analysis <- function(datasets, disease_params, predictors = c("mss_std", "MRS_std")) {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    for (pred in predictors) {
      if (!pred %in% names(data)) next
      
      data$neg_pred <- -data[[pred]]
      
      for (model_name in names(model_adjustments)) {
        adj_vars <- model_adjustments[[model_name]]
        adj_vars <- sub("fh_var", current$fh, adj_vars)
        adj_vars <- intersect(adj_vars, names(data))
        
        rhs <- build_formula(
          predictor = "neg_pred",
          adjustments = adj_vars,
          categorical_covariates = categorical_covariates,
          is_categorical_predictor = FALSE
        )
        
        fit <- tryCatch(
          fit_cox_model(data, current$time, current$event, rhs),
          error = function(e) NULL
        )
        
        if (is.null(fit)) next
        
        est <- extract_cox_term(fit, "neg_pred")
        
        results[[paste(disease, pred, model_name, sep = "_")]] <- data.frame(
          Disease = disease,
          Predictor = pred,
          Model = model_name,
          Term = "per SD decrement",
          HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
          P_Value = format_pvalue(est$p_value),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  bind_rows(results)
}

mrs_disease_results <- run_continuous_analysis(
  datasets = analysis_datasets,
  disease_params = disease_params,
  predictors = intersect(c("mss_std", "MRS_std"), names(final_mrs))
)

safe_write_csv(
  mrs_disease_results,
  file.path(dir_tables, "mrs_disease_association_results.csv")
)

# ---- 11. Optional mediation analysis ----------------------------------------

if (requireNamespace("CMAverse", quietly = TRUE)) {
  message("CMAverse detected: mediation analysis can be performed if desired.")
  
  # This section is intentionally conservative for repository use.
  # Users may adapt the example below to their locally approved workflow.
  
  run_mediation <- function(data, time_var, event_var, mediator, base_covars) {
    tryCatch({
      CMAverse::cmest(
        data = data,
        model = "rb",
        outcome = time_var,
        event = event_var,
        exposure = "mss",
        mediator = mediator,
        basec = base_covars,
        EMint = FALSE,
        mreg = list("linear"),
        yreg = "coxph",
        astar = 0,
        a = 1,
        mval = list(1),
        estimation = "imputation",
        inference = "bootstrap",
        nboot = 200
      )
    }, error = function(e) NULL)
  }
  
  mediation_results <- list()
  
  for (disease in names(analysis_datasets)) {
    data <- analysis_datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c("mss", "MRS_std", current$time, current$event) %in% names(data))) next
    
    base_covars <- intersect(
      c("age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
        "ms", "smoke", "drink", "mvpa3", current$fh),
      names(data)
    )
    
    data_med <- data %>%
      select(any_of(c("mss", "MRS_std", current$time, current$event, base_covars))) %>%
      na.omit()
    
    if (nrow(data_med) < 100) next
    
    fit_med <- run_mediation(
      data = data_med,
      time_var = current$time,
      event_var = current$event,
      mediator = "MRS_std",
      base_covars = base_covars
    )
    
    if (is.null(fit_med)) next
    
    mediation_results[[disease]] <- data.frame(
      Disease = disease,
      Mediator = "MRS_std",
      Rte = fit_med$effect.pe["Rte"],
      Rpnie = fit_med$effect.pe["Rpnie"],
      pm = fit_med$effect.pe["pm"],
      stringsAsFactors = FALSE
    )
  }
  
  if (length(mediation_results) > 0) {
    safe_write_csv(
      bind_rows(mediation_results),
      file.path(dir_tables, "mrs_mediation_summary.csv")
    )
  }
} else {
  message("CMAverse not installed; skipping mediation analysis.")
}

message("MRS construction and analysis complete.")
