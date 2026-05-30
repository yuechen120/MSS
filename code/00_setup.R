# =============================================================================
# MSS project: setup and shared utilities
# File: code/00_setup.R
#
# Purpose:
#   - Load required packages
#   - Define project directories
#   - Provide shared helper functions used across analysis scripts
#
# Notes:
#   - This repository does NOT include restricted participant-level data.
#   - Users must obtain access to the original datasets separately and update
#     file paths as needed.
#   - Analyses were developed in R 4.4.3.
# =============================================================================

rm(list = ls())

# ---- 1. Package management ---------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "janitor",
  "rlang",
  "survival",
  "mice",
  "ranger",
  "openxlsx",
  "ggplot2",
  "pROC",
  "glmnet",
  "car",
  "Hmisc",
  "ComplexHeatmap",
  "circlize",
  "lightgbm",
  "shapviz"
)

load_required_packages <- function(pkgs, install_if_missing = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    message("Missing packages: ", paste(missing_pkgs, collapse = ", "))
    if (isTRUE(install_if_missing)) {
      install.packages(missing_pkgs)
    } else {
      message(
        "Please install missing packages manually before running the scripts.\n",
        "If needed, rerun with install_if_missing = TRUE."
      )
    }
  }
  
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }))
}

load_required_packages(required_packages, install_if_missing = FALSE)

# ---- 2. Project directories --------------------------------------------------

# Assumes this script is sourced from within the repository.
# Example:
# source("code/00_setup.R")

project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)

dir_code   <- file.path(project_root, "code")
dir_data   <- file.path(project_root, "data")
dir_raw    <- file.path(dir_data, "raw")
dir_derived <- file.path(dir_data, "derived")
dir_output <- file.path(project_root, "output")
dir_docs   <- file.path(project_root, "docs")
dir_figures <- file.path(dir_output, "figures")
dir_tables  <- file.path(dir_output, "tables")

dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_docs, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_root)

# ---- 3. Data access note -----------------------------------------------------

message(
  "\nData access note:\n",
  "- Raw participant-level data are not included in this repository.\n",
  "- Update file paths in downstream scripts to point to your authorized local data sources.\n"
)

# ---- 4. Helper functions -----------------------------------------------------

format_pvalue <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_hr_ci <- function(hr, ci_low, ci_high, digits = 2) {
  if (any(is.na(c(hr, ci_low, ci_high)))) return(NA_character_)
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    hr, ci_low, ci_high
  )
}

format_beta_ci <- function(beta, ci_low, ci_high, digits = 3) {
  if (any(is.na(c(beta, ci_low, ci_high)))) return(NA_character_)
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    beta, ci_low, ci_high
  )
}

check_required_columns <- function(data, required_cols, object_name = "data") {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", object_name, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

deduplicate_by_id <- function(data, id_col = "f.eid") {
  check_required_columns(data, id_col, deparse(substitute(data)))
  data %>% distinct(.data[[id_col]], .keep_all = TRUE)
}

report_duplicates <- function(data, id_col = "f.eid", object_name = deparse(substitute(data))) {
  check_required_columns(data, id_col, object_name)
  n_dup <- data %>%
    count(.data[[id_col]]) %>%
    filter(n > 1) %>%
    nrow()
  
  if (n_dup > 0) {
    message(object_name, ": ", n_dup, " duplicated IDs detected.")
  } else {
    message(object_name, ": no duplicated IDs detected.")
  }
}

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, ...)
}

safe_write_csv <- function(data, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path, ...)
}

safe_write_xlsx <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(x, file = path, ...)
}

to_factor_if_present <- function(data, vars) {
  vars_present <- intersect(vars, names(data))
  data %>%
    mutate(across(all_of(vars_present), as.factor))
}

scale_numeric <- function(x) {
  as.numeric(scale(x))
}

log_scale <- function(x) {
  log_x <- ifelse(x > 0, log(x), NA_real_)
  as.numeric(scale(log_x))
}

build_formula <- function(predictor, adjustments = NULL,
                          categorical_covariates = NULL,
                          is_categorical_predictor = FALSE) {
  adj_terms <- character(0)
  
  if (!is.null(adjustments) && length(adjustments) > 0) {
    adj_terms <- vapply(adjustments, function(x) {
      if (!is.null(categorical_covariates) && x %in% categorical_covariates) {
        paste0("factor(", x, ")")
      } else {
        x
      }
    }, character(1))
  }
  
  pred_term <- if (isTRUE(is_categorical_predictor)) {
    paste0("factor(", predictor, ")")
  } else {
    predictor
  }
  
  paste(c(pred_term, adj_terms), collapse = " + ")
}

build_surv_formula <- function(time_var, event_var, rhs_terms) {
  as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs_terms)
  )
}

make_decrement_var <- function(data, varname, binary_vars = NULL) {
  stopifnot(varname %in% names(data))
  new_name <- paste0("neg_", varname)
  
  if (!is.null(binary_vars) && varname %in% binary_vars) {
    data[[new_name]] <- ifelse(is.na(data[[varname]]), NA, 1 - data[[varname]])
  } else {
    data[[new_name]] <- ifelse(is.na(data[[varname]]), NA, -1 * data[[varname]])
  }
  
  list(data = data, new_var = new_name)
}

fit_cox_model <- function(data, time_var, event_var, rhs_formula) {
  full_formula <- build_surv_formula(time_var, event_var, rhs_formula)
  survival::coxph(full_formula, data = data)
}

extract_cox_term <- function(model, term_name) {
  sm <- summary(model)
  
  if (!term_name %in% rownames(sm$coefficients)) {
    return(data.frame(
      term = term_name,
      beta = NA_real_,
      hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  data.frame(
    term = term_name,
    beta = sm$coefficients[term_name, "coef"],
    hr = sm$coefficients[term_name, "exp(coef)"],
    ci_low = sm$conf.int[term_name, "lower .95"],
    ci_high = sm$conf.int[term_name, "upper .95"],
    p_value = sm$coefficients[term_name, "Pr(>|z|)"]
  )
}

calc_incidence_summary <- function(data, group_var, time_var, event_var) {
  check_required_columns(data, c(group_var, time_var, event_var))
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      N = n(),
      Cases = sum(.data[[event_var]], na.rm = TRUE),
      Person_Years = sum(.data[[time_var]], na.rm = TRUE) / 365.25,
      Incidence_Rate = (Cases / Person_Years) * 1000,
      .groups = "drop"
    )
}

save_session_info <- function(path = file.path(dir_docs, "sessionInfo.txt")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), con = path)
}

# ---- 5. Common analysis settings --------------------------------------------

set.seed(20250326)

categorical_covariates_default <- c(
  "sex", "ethnic", "urban", "edu", "employ", "tdi3",
  "smoke", "drink", "mvpa3", "fh_var"
)

base_adjustments_default <- c("age", "sex", "ethnic")

model2_adjustments_default <- c(
  base_adjustments_default,
  "urban", "edu", "employ", "tdi3",
  "ms", "smoke", "drink", "mvpa3", "fh_var"
)

model3_adjustments_default <- c(
  model2_adjustments_default,
  "bmi", "ldl", "glucose", "sys"
)

# ---- 6. Save session info ----------------------------------------------------

save_session_info()

message("Setup complete.")
