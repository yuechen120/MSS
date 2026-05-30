# =============================================================================
# MSS project: sex-stratified and interaction analyses
# File: code/03_sex_stratified_and_interaction.R
#
# Purpose:
#   - Conduct sex-stratified Cox regression analyses
#   - Test interaction by sex for MSS and disease outcomes
#   - Export sex-disaggregated sample numbers and model results
#
# Notes:
#   - Sex is treated as a biological variable as recorded in the source dataset.
#   - This script assumes that final_imp.csv has already been created.
# =============================================================================

source("code/00_setup.R")

# ---- 1. Load analysis-ready dataset ------------------------------------------

file_final_imp <- file.path(dir_derived, "final_imp.csv")
final_imp <- safe_read_csv(file_final_imp)

message("Loaded final_imp with N = ", nrow(final_imp))

# ---- 2. Basic checks and derived variables -----------------------------------

required_core_vars <- c("sex", "mss", "mss3")
check_required_columns(final_imp, required_core_vars, "final_imp")

if (!"mss_std" %in% names(final_imp) && "mss" %in% names(final_imp)) {
  final_imp$mss_std <- scale_numeric(final_imp$mss)
}

# Standardize sex coding if possible
# Assumption: sex is coded as 0/1 or 1/0; adjust labels below if needed.
# You should revise this block if your local coding differs.
if (is.numeric(final_imp$sex) || is.integer(final_imp$sex)) {
  sex_levels_present <- sort(unique(final_imp$sex[!is.na(final_imp$sex)]))
  message("Observed sex coding: ", paste(sex_levels_present, collapse = ", "))
}

# ---- 3. Disease map ----------------------------------------------------------

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

# ---- 4. Model settings -------------------------------------------------------

categorical_covariates <- c(
  "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "smoke", "drink", "mvpa3", "fh_var"
)

model3_adjustments <- c(
  "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "ms", "smoke", "drink", "mvpa3", "fh_var", "bmi", "ldl", "glucose", "sys"
)

# ---- 5. Sex-disaggregated counts --------------------------------------------

get_sex_counts <- function(data, disease_name, params) {
  data %>%
    group_by(sex) %>%
    summarise(
      N = n(),
      Cases = sum(.data[[params$event]], na.rm = TRUE),
      Person_Years = sum(.data[[params$time]], na.rm = TRUE) / 365.25,
      .groups = "drop"
    ) %>%
    mutate(
      Disease = disease_name
    ) %>%
    select(Disease, sex, N, Cases, Person_Years)
}

sex_counts <- bind_rows(
  lapply(names(analysis_datasets), function(disease) {
    get_sex_counts(analysis_datasets[[disease]], disease, disease_params[[disease]])
  })
)

safe_write_csv(
  sex_counts,
  file.path(dir_tables, "sex_disaggregated_sample_counts.csv")
)

# ---- 6. Sex-stratified analysis ---------------------------------------------

run_sex_stratified_analysis <- function(datasets, disease_params, predictor = "mss_std") {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c("sex", predictor, current$time, current$event) %in% names(data))) next
    
    sex_values <- sort(unique(data$sex[!is.na(data$sex)]))
    
    for (sx in sex_values) {
      data_sub <- data %>% filter(sex == sx)
      
      if (nrow(data_sub) == 0) next
      
      adj_vars <- model3_adjustments
      adj_vars <- sub("fh_var", current$fh, adj_vars)
      
      # remove sex from stratified model
      adj_vars <- setdiff(adj_vars, "sex")
      
      # create decrement variable
      data_sub$neg_mss_std <- -data_sub[[predictor]]
      
      rhs <- build_formula(
        predictor = "neg_mss_std",
        adjustments = adj_vars,
        categorical_covariates = categorical_covariates,
        is_categorical_predictor = FALSE
      )
      
      fit <- tryCatch(
        fit_cox_model(data_sub, current$time, current$event, rhs),
        error = function(e) NULL
      )
      
      if (is.null(fit)) next
      
      est <- extract_cox_term(fit, "neg_mss_std")
      
      results[[paste(disease, sx, sep = "_")]] <- data.frame(
        Disease = disease,
        Sex = sx,
        N = nrow(data_sub),
        Cases = sum(data_sub[[current$event]], na.rm = TRUE),
        Predictor = predictor,
        Model = "model3",
        Term = "per SD decrement",
        HR_95CI = format_hr_ci(est$hr, est$ci_low, est$ci_high),
        P_Value = format_pvalue(est$p_value),
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(results)
}

sex_stratified_results <- run_sex_stratified_analysis(
  datasets = analysis_datasets,
  disease_params = disease_params,
  predictor = "mss_std"
)

safe_write_csv(
  sex_stratified_results,
  file.path(dir_tables, "sex_stratified_mss_results.csv")
)

# ---- 7. Interaction by sex ---------------------------------------------------

run_interaction_analysis <- function(datasets, disease_params, predictor = "mss_std") {
  results <- list()
  
  for (disease in names(datasets)) {
    data <- datasets[[disease]]
    current <- disease_params[[disease]]
    
    if (!all(c("sex", predictor, current$time, current$event) %in% names(data))) next
    
    adj_vars <- model3_adjustments
    adj_vars <- sub("fh_var", current$fh, adj_vars)
    
    # remove sex main effect from adjustment list because it will be included explicitly
    adj_vars_no_sex <- setdiff(adj_vars, "sex")
    
    data$neg_mss_std <- -data[[predictor]]
    
    adj_terms <- vapply(adj_vars_no_sex, function(x) {
      if (x %in% categorical_covariates) {
        paste0("factor(", x, ")")
      } else {
        x
      }
    }, character(1))
    
    rhs <- paste(
      c("neg_mss_std * factor(sex)", adj_terms),
      collapse = " + "
    )
    
    fit <- tryCatch(
      fit_cox_model(data, current$time, current$event, rhs),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      results[[disease]] <- data.frame(
        Disease = disease,
        Interaction_Term = "neg_mss_std:factor(sex)",
        P_Value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    sm <- summary(fit)
    rn <- rownames(sm$coefficients)
    
    interaction_row <- grep("^neg_mss_std:factor\\(sex\\)", rn, value = TRUE)
    
    p_int <- if (length(interaction_row) > 0) {
      sm$coefficients[interaction_row[1], "Pr(>|z|)"]
    } else {
      NA_real_
    }
    
    results[[disease]] <- data.frame(
      Disease = disease,
      Interaction_Term = "neg_mss_std × sex",
      P_Value = p_int,
      stringsAsFactors = FALSE
    )
  }
  
  out <- bind_rows(results)
  if (nrow(out) > 0) {
    out$FDR_P_Value <- p.adjust(out$P_Value, method = "BH")
  }
  out
}

interaction_results <- run_interaction_analysis(
  datasets = analysis_datasets,
  disease_params = disease_params,
  predictor = "mss_std"
)

safe_write_csv(
  interaction_results,
  file.path(dir_tables, "interaction_results_by_sex.csv")
)

# ---- 8. Summary file for source data / supplement ----------------------------

sex_summary <- sex_stratified_results %>%
  left_join(
    interaction_results %>% select(Disease, Interaction_P = P_Value, Interaction_FDR = FDR_P_Value),
    by = "Disease"
  )

safe_write_csv(
  sex_summary,
  file.path(dir_tables, "sex_stratified_results_with_interaction.csv")
)

message("Sex-stratified and interaction analyses complete.")