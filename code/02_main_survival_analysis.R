# =============================================================================
# MSS project: main survival analysis
# File: code/02_main_survival_analysis.R
#
# Purpose:
#   - Load the analysis-ready dataset
#   - Build disease-specific analytic datasets
#   - Summarize baseline characteristics
#   - Run main Cox regression analyses for MSS and HSS
#   - Calculate incidence rates
#   - Compare MSS and HSS model performance
#
# Required input:
#   - data/derived/final_imp.csv
#
# Outputs:
#   - output/tables/table1_baseline_mss3.csv
#   - output/tables/cox_results_mss.csv
#   - output/tables/cox_results_hss.csv
#   - output/tables/incidence_by_mss3.csv
#   - output/tables/model_comparison_mss_vs_hss.csv
# =============================================================================

source("code/00_setup.R")

# ---- 1. Load analysis-ready dataset ------------------------------------------

file_final_imp <- file.path(dir_derived, "final_imp.csv")
final_imp <- safe_read_csv(file_final_imp)

message("Loaded final_imp with N = ", nrow(final_imp))

# ---- 2. Derived variables ----------------------------------------------------

if ("mss" %in% names(final_imp) && !"mss_std" %in% names(final_imp)) {
  final_imp$mss_std <- scale_numeric(final_imp$mss)
}

if ("hss" %in% names(final_imp) && !"hss_std" %in% names(final_imp)) {
  final_imp$hss_std <- scale_numeric(final_imp$hss)
}

# ---- 3. Disease map and disease-specific datasets ----------------------------

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

build_disease_dataset <- function(data, disease_name, params) {
  check_required_columns(data, c(params$time, params$event))
  
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
  build_disease_dataset(final_imp, dz, disease_params[[dz]])
})
names(analysis_datasets) <- names(disease_params)

invisible(lapply(names(analysis_datasets), function(nm) {
  message("Dataset ", nm, ": N = ", nrow(analysis_datasets[[nm]]))
}))

# ---- 4. Baseline characteristics (Table 1 style) ----------------------------

vars_to_describe <- c(
  "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3", "light3", "noise3",
  "ms", "smoke", "drink", "coffee3", "tea3", "mvpa3",
  "insomnia", "doze", "snore", "regularity", "epijetlag", "socialjetlag",
  "onset", "efficiency", "daytime_nap", "duration", "duration_report",
  "mss", "hss", "bmi", "sys", "ldl", "glucose"
)

vars_to_describe <- intersect(vars_to_describe, names(final_imp))

table1_data <- final_imp %>%
  mutate(
    mss3 = as.factor(mss3)
  ) %>%
  to_factor_if_present(
    vars = c(
      "sex", "ethnic", "urban", "edu", "employ", "tdi3", "light3", "noise3",
      "smoke", "drink", "coffee3", "tea3", "mvpa3",
      "insomnia", "doze", "snore", "regularity", "epijetlag", "socialjetlag",
      "onset", "efficiency", "daytime_nap", "duration", "duration_report"
    )
  )

summarize_numeric <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

summarize_factor <- function(x) {
  tb <- table(x, useNA = "no")
  pct <- prop.table(tb) * 100
  paste0(names(tb), "=", as.integer(tb), " (", sprintf("%.1f", pct), "%)", collapse = " | ")
}

make_baseline_table <- function(data, group_var, vars) {
  grouped <- split(data, data[[group_var]])
  
  res_list <- lapply(vars, function(v) {
    x <- data[[v]]
    
    if (is.numeric(x)) {
      vals <- sapply(grouped, function(df) summarize_numeric(df[[v]]))
      pooled <- summarize_numeric(x)
      out <- data.frame(
        Variable = v,
        Subcategory = NA_character_,
        Pooled = pooled,
        stringsAsFactors = FALSE
      )
      for (nm in names(vals)) out[[paste0("Group_", nm)]] <- vals[[nm]]
      out
    } else {
      lvls <- levels(as.factor(x))
      out <- lapply(lvls, function(lv) {
        pooled_n <- sum(x == lv, na.rm = TRUE)
        pooled_pct <- mean(x == lv, na.rm = TRUE) * 100
        
        row <- data.frame(
          Variable = v,
          Subcategory = lv,
          Pooled = paste0(pooled_n, " (", sprintf("%.1f", pooled_pct), "%)"),
          stringsAsFactors = FALSE
        )
        
        for (nm in names(grouped)) {
          vec <- grouped[[nm]][[v]]
          n_lv <- sum(vec == lv, na.rm = TRUE)
          pct_lv <- mean(vec == lv, na.rm = TRUE) * 100
          row[[paste0("Group_", nm)]] <- paste0(n_lv, " (", sprintf("%.1f", pct_lv), "%)")
        }
        
        row
      })
      bind_rows(out)
    }
  })
  
  bind_rows(res_list)
}

table1_out <- make_baseline_table(table1_data, "mss3", vars_to_describe)

safe_write_csv(table1_out, file.path(dir_tables, "table1_baseline_mss3.csv"))

# ---- 5. Cox model settings ---------------------------------------------------

categorical_covariates <- c(
  "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "smoke", "drink", "mvpa3", "fh_var"
)

base_adjustments <- c("age", "sex", "ethnic")

model2_common <- c(
  base_adjustments,
  "urban", "edu", "employ", "tdi3",
  "ms", "smoke", "drink", "mvpa3", "fh_var"
)

model3_common <- c(
  model2_common,
  "bmi", "ldl", "glucose", "sys"
)

model_definitions <- list(
  model1 = list(adjustments = base_adjustments),
  model2 = list(adjustments = model2_common),
  model3 = list(adjustments = model3_common)
)

# ---- 6. Main MSS analysis ----------------------------------------------------

run_main_cox_analysis <- function(datasets, disease_params, predictor_cat, predictor_cont, predictor_std) {
  results_cat <- list()
  results_cont <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    # categorical predictor
    if (predictor_cat %in% names(data)) {
      data[[predictor_cat]] <- factor(data[[predictor_cat]])
      if ("2" %in% levels(data[[predictor_cat]])) {
        data[[predictor_cat]] <- relevel(data[[predictor_cat]], ref = "2")
      }
      
      for (model_name in names(model_definitions)) {
        adj_vars <- model_definitions[[model_name]]$adjustments
        if ("fh_var" %in% adj_vars) {
          adj_vars <- sub("fh_var", current$fh, adj_vars)
        }
        
        rhs_cat <- build_formula(
          predictor = predictor_cat,
          adjustments = adj_vars,
          categorical_covariates = categorical_covariates,
          is_categorical_predictor = TRUE
        )
        
        rhs_trend <- build_formula(
          predictor = paste0("as.numeric(", predictor_cat, ")"),
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
        
        if (!is.null(fit_cat)) {
          sm <- summary(fit_cat)
          coef_names <- rownames(sm$coefficients)
          pred_terms <- grep(paste0("^factor\\(", predictor_cat, "\\)"), coef_names, value = TRUE)
          
          row_out <- data.frame(
            Disease = disease,
            Predictor = predictor_cat,
            Model = model_name,
            Reference = "1.00 (ref)",
            Level1_HR_95CI = NA_character_,
            Level1_P = NA_character_,
            Level2_HR_95CI = NA_character_,
            Level2_P = NA_character_,
            P_trend = NA_real_,
            stringsAsFactors = FALSE
          )
          
          for (i in seq_along(pred_terms)) {
            term <- pred_terms[i]
            est <- extract_cox_term(fit_cat, term)
            hr_ci <- format_hr_ci(est$hr, est$ci_low, est$ci_high)
            pval <- format_pvalue(est$p_value)
            
            if (i == 1) {
              row_out$Level1_HR_95CI <- hr_ci
              row_out$Level1_P <- pval
            } else if (i == 2) {
              row_out$Level2_HR_95CI <- hr_ci
              row_out$Level2_P <- pval
            }
          }
          
          if (!is.null(fit_trend)) {
            trend_sm <- summary(fit_trend)
            row_out$P_trend <- trend_sm$coefficients[1, "Pr(>|z|)"]
          }
          
          results_cat[[paste(disease, predictor_cat, model_name, sep = "_")]] <- row_out
        }
      }
    }
    
    # continuous predictors
    for (pred in c(predictor_cont, predictor_std)) {
      if (!pred %in% names(data)) next
      
      for (model_name in names(model_definitions)) {
        adj_vars <- model_definitions[[model_name]]$adjustments
        if ("fh_var" %in% adj_vars) {
          adj_vars <- sub("fh_var", current$fh, adj_vars)
        }
        
        data2 <- make_decrement_var(data, pred)$data
        neg_var <- make_decrement_var(data, pred)$new_var
        
        rhs_cont <- build_formula(
          predictor = neg_var,
          adjustments = adj_vars,
          categorical_covariates = categorical_covariates,
          is_categorical_predictor = FALSE
        )
        
        fit_cont <- tryCatch(
          fit_cox_model(data2, current$time, current$event, rhs_cont),
          error = function(e) NULL
        )
        
        if (!is.null(fit_cont)) {
          est <- extract_cox_term(fit_cont, neg_var)
          
          results_cont[[paste(disease, pred, model_name, sep = "_")]] <- data.frame(
            Disease = disease,
            Predictor = pred,
            Model = model_name,
            Term = ifelse(grepl("_std$", pred), "per SD decrement", "per point decrement"),
            HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
            P_Value = format_pvalue(est$p_value),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  list(
    categorical = bind_rows(results_cat),
    continuous = bind_rows(results_cont)
  )
}

mss_results <- run_main_cox_analysis(
  datasets = analysis_datasets,
  disease_params = disease_params,
  predictor_cat = "mss3",
  predictor_cont = "mss",
  predictor_std = "mss_std"
)

if (nrow(mss_results$categorical) > 0) {
  mss_results$categorical$P_trend_FDR <- p.adjust(mss_results$categorical$P_trend, method = "fdr")
}
safe_write_csv(mss_results$categorical, file.path(dir_tables, "cox_results_mss_categorical.csv"))
safe_write_csv(mss_results$continuous, file.path(dir_tables, "cox_results_mss_continuous.csv"))

# ---- 7. Main HSS analysis ----------------------------------------------------

hss_results <- run_main_cox_analysis(
  datasets = analysis_datasets,
  disease_params = disease_params,
  predictor_cat = "hss3",
  predictor_cont = "hss",
  predictor_std = "hss_std"
)

if (nrow(hss_results$categorical) > 0) {
  hss_results$categorical$P_trend_FDR <- p.adjust(hss_results$categorical$P_trend, method = "fdr")
}
safe_write_csv(hss_results$categorical, file.path(dir_tables, "cox_results_hss_categorical.csv"))
safe_write_csv(hss_results$continuous, file.path(dir_tables, "cox_results_hss_continuous.csv"))

# ---- 8. Incidence rate by MSS category --------------------------------------

incidence_results <- lapply(names(analysis_datasets), function(disease) {
  data <- analysis_datasets[[disease]]
  current <- disease_params[[disease]]
  
  if (!all(c("mss3", current$time, current$event) %in% names(data))) return(NULL)
  
  calc_incidence_summary(
    data = data,
    group_var = "mss3",
    time_var = current$time,
    event_var = current$event
  ) %>%
    mutate(
      Disease = disease
    ) %>%
    select(Disease, everything())
})

incidence_results <- bind_rows(incidence_results)
safe_write_csv(incidence_results, file.path(dir_tables, "incidence_by_mss3.csv"))

# ---- 9. Model performance: MSS vs HSS ---------------------------------------

get_cindex <- function(model, data, time_var, event_var) {
  pred <- tryCatch(predict(model), error = function(e) NULL)
  if (is.null(pred)) return(NA_real_)
  
  tryCatch({
    survcomp::concordance.index(
      pred,
      surv.time = data[[time_var]],
      surv.event = data[[event_var]],
      method = "noether"
    )$c.index
  }, error = function(e) NA_real_)
}

get_auc5 <- function(model, data, time_var, event_var) {
  pred <- tryCatch(predict(model, type = "lp"), error = function(e) NULL)
  if (is.null(pred)) return(NA_real_)
  
  tryCatch({
    roc_obj <- timeROC::timeROC(
      T = data[[time_var]],
      delta = data[[event_var]],
      marker = pred,
      cause = 1,
      weighting = "marginal",
      times = 1825,
      ROC = FALSE
    )
    roc_obj$AUC[2]
  }, error = function(e) NA_real_)
}

compare_mss_hss <- function(datasets, disease_params) {
  res <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c("mss_std", "hss_std") %in% names(data))) next
    
    covariates <- c(
      "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
      "ms", "smoke", "drink", "mvpa3", current$fh, "bmi", "ldl", "glucose", "sys"
    )
    covariates <- intersect(covariates, names(data))
    
    rhs_mss <- build_formula(
      predictor = "neg_mss_std",
      adjustments = covariates,
      categorical_covariates = categorical_covariates,
      is_categorical_predictor = FALSE
    )
    
    rhs_hss <- build_formula(
      predictor = "neg_hss_std",
      adjustments = covariates,
      categorical_covariates = categorical_covariates,
      is_categorical_predictor = FALSE
    )
    
    data$neg_mss_std <- -data$mss_std
    data$neg_hss_std <- -data$hss_std
    
    fit_mss <- tryCatch(
      fit_cox_model(data, current$time, current$event, rhs_mss),
      error = function(e) NULL
    )
    
    fit_hss <- tryCatch(
      fit_cox_model(data, current$time, current$event, rhs_hss),
      error = function(e) NULL
    )
    
    if (is.null(fit_mss) || is.null(fit_hss)) next
    
    est_mss <- extract_cox_term(fit_mss, "neg_mss_std")
    est_hss <- extract_cox_term(fit_hss, "neg_hss_std")
    
    res[[disease]] <- data.frame(
      Disease = disease,
      MSS_HR_95CI = format_hr_ci(est_mss$hr, est_mss$ci_low, est_mss$ci_high),
      MSS_P = format_pvalue(est_mss$p_value),
      HSS_HR_95CI = format_hr_ci(est_hss$hr, est_hss$ci_low, est_hss$ci_high),
      HSS_P = format_pvalue(est_hss$p_value),
      MSS_AIC = tryCatch(extractAIC(fit_mss)[2], error = function(e) NA_real_),
      HSS_AIC = tryCatch(extractAIC(fit_hss)[2], error = function(e) NA_real_),
      MSS_Cindex = get_cindex(fit_mss, data, current$time, current$event),
      HSS_Cindex = get_cindex(fit_hss, data, current$time, current$event),
      MSS_AUC5 = get_auc5(fit_mss, data, current$time, current$event),
      HSS_AUC5 = get_auc5(fit_hss, data, current$time, current$event),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(res)
}

# Optional packages used here
if (requireNamespace("survcomp", quietly = TRUE) &&
    requireNamespace("timeROC", quietly = TRUE)) {
  model_comparison <- compare_mss_hss(analysis_datasets, disease_params)
  safe_write_csv(model_comparison, file.path(dir_tables, "model_comparison_mss_vs_hss.csv"))
} else {
  message("Skipping model comparison: packages 'survcomp' and/or 'timeROC' are not installed.")
}

message("Main survival analysis complete.")
