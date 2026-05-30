# =============================================================================
# MSS project: sensitivity analyses
# File: code/04_sensitivity_analyses.R
#
# Purpose:
#   - Conduct sensitivity analyses for the association between sleep health and
#     disease outcomes
#   - Perform lag analysis
#   - Adjust for interval between baseline assessment and accelerometer date
#   - Evaluate additional sleep metrics
#   - Examine consistency/discordance between subjective and objective sleep
#
# Required input:
#   - data/derived/final_imp.csv
#
# Optional input:
#   - baseline assessment date file, if available and permitted
# =============================================================================

source("code/00_setup.R")

# ---- 1. Load analysis-ready dataset ------------------------------------------

file_final_imp <- file.path(dir_derived, "final_imp.csv")
final_imp <- safe_read_csv(file_final_imp)

message("Loaded final_imp with N = ", nrow(final_imp))

if (!"mss_std" %in% names(final_imp) && "mss" %in% names(final_imp)) {
  final_imp$mss_std <- scale_numeric(final_imp$mss)
}

# ---- 2. Disease map and shared settings --------------------------------------

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

model3_adjustments <- c(
  "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "ms", "smoke", "drink", "mvpa3", "fh_var", "bmi", "ldl", "glucose", "sys"
)

build_disease_dataset <- function(data, disease_name, params, lag_days = 0) {
  out <- data %>%
    filter(!( .data[[params$event]] == 1 & .data[[params$time]] < lag_days ))
  
  if (disease_name == "t2d" && "glucose" %in% names(out)) {
    out <- out %>% filter(glucose < 11 | is.na(glucose))
  }
  
  if (disease_name == "obese" && "bmi" %in% names(out)) {
    out <- out %>% filter(bmi < 30 | is.na(bmi))
  }
  
  out
}

run_continuous_mss_model3 <- function(datasets, disease_params, predictor = "mss") {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c(predictor, current$time, current$event) %in% names(data))) next
    
    data$neg_pred <- -data[[predictor]]
    
    adj_vars <- sub("fh_var", current$fh, model3_adjustments)
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
    
    results[[disease]] <- data.frame(
      Disease = disease,
      Predictor = predictor,
      Model = "model3",
      Term = ifelse(grepl("_std$", predictor), "per SD decrement", "per point decrement"),
      N = nrow(data),
      Cases = sum(data[[current$event]], na.rm = TRUE),
      HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
      P_Value = format_pvalue(est$p_value),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(results)
}

run_categorical_mss3_model3 <- function(datasets, disease_params, predictor = "mss3", trend_var = "mss") {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c(predictor, trend_var, current$time, current$event) %in% names(data))) next
    
    data[[predictor]] <- factor(data[[predictor]])
    if ("2" %in% levels(data[[predictor]])) {
      data[[predictor]] <- relevel(data[[predictor]], ref = "2")
    }
    
    medians <- tapply(data[[trend_var]], data[[predictor]], median, na.rm = TRUE)
    data$trend_median <- medians[as.character(data[[predictor]])]
    
    adj_vars <- sub("fh_var", current$fh, model3_adjustments)
    adj_vars <- intersect(adj_vars, names(data))
    
    rhs_cat <- build_formula(
      predictor = predictor,
      adjustments = adj_vars,
      categorical_covariates = categorical_covariates,
      is_categorical_predictor = TRUE
    )
    
    rhs_trend <- build_formula(
      predictor = "trend_median",
      adjustments = adj_vars,
      categorical_covariates = categorical_covariates,
      is_categorical_predictor = FALSE
    )
    
    fit_cat <- tryCatch(
      fit_cox_model(data, current$time, current$event, rhs_cat),
      error = function(e) NULL
    )
    
    fit_trend <- tryCatch(
      fit_cox_model(data, current$time, current$event, rhs_trend),
      error = function(e) NULL
    )
    
    if (is.null(fit_cat)) next
    
    sm <- summary(fit_cat)
    coef_names <- rownames(sm$coefficients)
    pred_terms <- grep(paste0("^factor\\(", predictor, "\\)"), coef_names, value = TRUE)
    
    row_out <- data.frame(
      Disease = disease,
      Predictor = predictor,
      Model = "model3",
      Reference = "1.00 (ref)",
      Level1_HR_95CI = NA_character_,
      Level1_P = NA_character_,
      Level2_HR_95CI = NA_character_,
      Level2_P = NA_character_,
      P_trend = if (!is.null(fit_trend)) summary(fit_trend)$coefficients[1, "Pr(>|z|)"] else NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(pred_terms)) {
      term <- pred_terms[i]
      est <- extract_cox_term(fit_cat, term)
      
      if (i == 1) {
        row_out$Level1_HR_95CI <- format_hr_ci(est$hr, est$ci_low, est$ci_high)
        row_out$Level1_P <- format_pvalue(est$p_value)
      } else if (i == 2) {
        row_out$Level2_HR_95CI <- format_hr_ci(est$hr, est$ci_low, est$ci_high)
        row_out$Level2_P <- format_pvalue(est$p_value)
      }
    }
    
    results[[disease]] <- row_out
  }
  
  out <- bind_rows(results)
  if (nrow(out) > 0) {
    out$P_trend_FDR <- p.adjust(out$P_trend, method = "fdr")
  }
  out
}

# ---- 3. Sensitivity analysis 1: lag analysis --------------------------------

lag_days <- 730.5

lag_datasets <- lapply(names(disease_params), function(dz) {
  build_disease_dataset(final_imp, dz, disease_params[[dz]], lag_days = lag_days)
})
names(lag_datasets) <- names(disease_params)

lag_cat_results <- run_categorical_mss3_model3(
  datasets = lag_datasets,
  disease_params = disease_params,
  predictor = "mss3",
  trend_var = "mss"
)

lag_cont_results <- run_continuous_mss_model3(
  datasets = lag_datasets,
  disease_params = disease_params,
  predictor = "mss"
)

safe_write_csv(
  lag_cat_results,
  file.path(dir_tables, "sensitivity_lag_mss3_model3.csv")
)

safe_write_csv(
  lag_cont_results,
  file.path(dir_tables, "sensitivity_lag_mss_continuous_model3.csv")
)

# incidence summaries for lag analysis
lag_incidence <- bind_rows(
  lapply(names(lag_datasets), function(disease) {
    data <- lag_datasets[[disease]]
    current <- disease_params[[disease]]
    if (!all(c("mss3", current$time, current$event) %in% names(data))) return(NULL)
    
    calc_incidence_summary(
      data = data,
      group_var = "mss3",
      time_var = current$time,
      event_var = current$event
    ) %>%
      mutate(Disease = disease)
  })
)

safe_write_csv(
  lag_incidence,
  file.path(dir_tables, "sensitivity_lag_incidence_by_mss3.csv")
)

# ---- 4. Sensitivity analysis 2: adjust for interval --------------------------

# Optional restricted input:
# If available, load baseline assessment date and calculate interval between
# baseline assessment and accelerometer measurement date.

file_baseline <- file.path(dir_raw, "ukb_baseline_assessment_date.csv")

if (file.exists(file_baseline)) {
  baseline_df <- safe_read_csv(file_baseline) %>%
    select(f.eid, baseline = f.53.0.0) %>%
    mutate(baseline = as.Date(baseline))
  
  data_interval <- final_imp %>%
    inner_join(baseline_df, by = "f.eid") %>%
    mutate(
      calendar_date = as.Date(calendar_date),
      interval = as.numeric(calendar_date - baseline)
    )
  
  message("Interval-adjusted dataset N = ", nrow(data_interval))
  
  model3_adjustments_interval <- c(
    "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
    "ms", "smoke", "drink", "mvpa3", "fh_var", "bmi", "ldl", "glucose", "sys", "interval"
  )
  
  run_continuous_with_interval <- function(data, disease_name, params) {
    if (!all(c("mss", params$time, params$event, "interval") %in% names(data))) return(NULL)
    
    data <- build_disease_dataset(data, disease_name, params, lag_days = 0)
    data$neg_mss <- -data$mss
    
    adj_vars <- sub("fh_var", params$fh, model3_adjustments_interval)
    adj_vars <- intersect(adj_vars, names(data))
    
    rhs <- build_formula(
      predictor = "neg_mss",
      adjustments = adj_vars,
      categorical_covariates = categorical_covariates,
      is_categorical_predictor = FALSE
    )
    
    fit <- tryCatch(
      fit_cox_model(data, params$time, params$event, rhs),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NULL)
    
    est <- extract_cox_term(fit, "neg_mss")
    
    data.frame(
      Disease = disease_name,
      Predictor = "mss",
      Model = "model3 + interval",
      Term = "per point decrement",
      N = nrow(data),
      Cases = sum(data[[params$event]], na.rm = TRUE),
      HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
      P_Value = format_pvalue(est$p_value),
      stringsAsFactors = FALSE
    )
  }
  
  interval_results <- bind_rows(
    lapply(names(disease_params), function(disease) {
      run_continuous_with_interval(data_interval, disease, disease_params[[disease]])
    })
  )
  
  safe_write_csv(
    interval_results,
    file.path(dir_tables, "sensitivity_interval_adjusted_mss.csv")
  )
} else {
  message("Baseline assessment date file not found; skipping interval-adjusted analysis.")
}

# ---- 5. Sensitivity analysis 3: additional sleep metrics ---------------------

additional_sleep_vars <- c(
  "insomnia", "doze", "snore", "daytime_nap", "regularity",
  "epijetlag", "preference", "socialjetlag", "onset",
  "efficiency", "duration", "duration_report"
)

additional_sleep_vars <- intersect(additional_sleep_vars, names(final_imp))

run_additional_sleep_analysis <- function(datasets, disease_params, predictors) {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    for (pred in predictors) {
      if (!pred %in% names(data)) next
      
      # Treat these variables as binary healthy/unhealthy indicators
      data$neg_pred <- ifelse(is.na(data[[pred]]), NA, 1 - data[[pred]])
      
      adj_vars <- sub("fh_var", current$fh, model3_adjustments)
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
      
      results[[paste(disease, pred, sep = "_")]] <- data.frame(
        Disease = disease,
        Predictor = pred,
        Model = "model3",
        Term = "decrement (healthy to unhealthy)",
        HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
        P_Value = format_pvalue(est$p_value),
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(results)
}

additional_sleep_results <- run_additional_sleep_analysis(
  datasets = lapply(names(disease_params), function(dz) {
    build_disease_dataset(final_imp, dz, disease_params[[dz]], lag_days = 0)
  }) %>% setNames(names(disease_params)),
  disease_params = disease_params,
  predictors = additional_sleep_vars
)

safe_write_csv(
  additional_sleep_results,
  file.path(dir_tables, "sensitivity_additional_sleep_metrics.csv")
)

# ---- 6. Sensitivity analysis 4: consistency / discordance --------------------

# These indicators compare subjective and objective sleep domains where available.
data_discordance <- final_imp

required_discordance_vars <- c(
  "duration", "duration_report", "insomnia", "onset", "doze", "daytime_nap",
  "mss", "hss"
)

if (all(required_discordance_vars %in% names(data_discordance))) {
  data_discordance <- data_discordance %>%
    mutate(
      dur_consistent = ifelse(duration == duration_report, 1, 0),
      ins_consistent = ifelse(insomnia == onset, 1, 0),
      doze_consistent = ifelse(doze == daytime_nap, 1, 0)
    )
  
  if (all(c("regularity", "socialjetlag", "efficiency", "epijetlag") %in% names(data_discordance))) {
    data_discordance <- data_discordance %>%
      mutate(
        mss_obj = regularity + socialjetlag + onset + efficiency + duration + daytime_nap + epijetlag,
        mss_sub = insomnia + doze + snore
      )
    
    mss_z <- scale_numeric(data_discordance$mss)
    hss_z <- scale_numeric(data_discordance$hss)
    mss_obj_z <- scale_numeric(data_discordance$mss_obj)
    mss_sub_z <- scale_numeric(data_discordance$mss_sub)
    
    data_discordance <- data_discordance %>%
      mutate(
        mss_hss_consistent_z = ifelse(abs(mss_z - hss_z) <= 1, 1, 0),
        obj_sub_consistent_z = ifelse(abs(mss_obj_z - mss_sub_z) <= 1, 1, 0)
      )
  }
  
  discordance_predictors <- intersect(
    c(
      "dur_consistent", "ins_consistent", "doze_consistent",
      "mss_hss_consistent_z", "obj_sub_consistent_z"
    ),
    names(data_discordance)
  )
  
  run_discordance_analysis <- function(data, disease_params, predictors) {
    results <- list()
    
    for (disease in names(disease_params)) {
      current <- disease_params[[disease]]
      data_dz <- build_disease_dataset(data, disease, current, lag_days = 0)
      
      for (pred in predictors) {
        if (!pred %in% names(data_dz)) next
        
        data_dz[[pred]] <- factor(data_dz[[pred]], levels = c(0, 1))
        
        adj_vars <- sub("fh_var", current$fh, model3_adjustments)
        adj_vars <- intersect(adj_vars, names(data_dz))
        
        rhs <- build_formula(
          predictor = pred,
          adjustments = adj_vars,
          categorical_covariates = categorical_covariates,
          is_categorical_predictor = TRUE
        )
        
        fit <- tryCatch(
          fit_cox_model(data_dz, current$time, current$event, rhs),
          error = function(e) NULL
        )
        
        if (is.null(fit)) next
        
        sm <- summary(fit)
        rn <- rownames(sm$coefficients)
        term <- grep(paste0("^factor\\(", pred, "\\)1$"), rn, value = TRUE)
        
        if (length(term) == 0) next
        
        est <- extract_cox_term(fit, term[1])
        
        results[[paste(disease, pred, sep = "_")]] <- data.frame(
          Disease = disease,
          Predictor = pred,
          Model = "model3",
          Reference = "0",
          Level1_HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
          Level1_P = format_pvalue(est$p_value),
          stringsAsFactors = FALSE
        )
      }
    }
    
    bind_rows(results)
  }
  
  discordance_results <- run_discordance_analysis(
    data = data_discordance,
    disease_params = disease_params,
    predictors = discordance_predictors
  )
  
  safe_write_csv(
    discordance_results,
    file.path(dir_tables, "sensitivity_discordance_results.csv")
  )
} else {
  message("Required discordance variables not available; skipping consistency/discordance analysis.")
}

message("Sensitivity analyses complete.")