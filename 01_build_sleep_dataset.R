# =============================================================================
# MSS project: build sleep dataset
# File: code/01_build_sleep_dataset.R
#
# Purpose:
#   - Construct accelerometer-derived and questionnaire-derived sleep variables
#   - Apply accelerometer QC
#   - Derive MSS and HSS
#   - Merge analytic covariate tables
#   - Generate analysis-ready datasets
#
# Important:
#   - This script requires authorized access to restricted source data.
#   - Replace placeholder file names and object-loading steps with your local setup.
# =============================================================================

source("code/00_setup.R")

# ---- 1. Input file placeholders ----------------------------------------------

# Replace these paths with your local authorized file locations.
file_acc_calibration <- file.path(dir_raw, "ukb_accelerometer_calibration.csv")
file_acc_weartime    <- file.path(dir_raw, "ukb_accelerometer_wear_time.csv")
file_acc_rawstats    <- file.path(dir_raw, "ukb_accelerometer_raw_statistics.csv")

file_part4_night     <- file.path(dir_raw, "part4_nightsummary_sleep_cleaned.csv")
file_part5_person    <- file.path(dir_raw, "part5_personsummary_sleep.csv")
file_sleep_question  <- file.path(dir_raw, "ukb_sleep_questionnaire.csv")

# Optional: downstream merged covariate tables can be loaded from pre-derived files
# or generated separately by authorized users.
# Examples:
# file_health         <- file.path(dir_derived, "health.csv")
# file_demo           <- file.path(dir_derived, "demo.csv")

# ---- 2. Load restricted source files -----------------------------------------

calib1 <- safe_read_csv(file_acc_calibration) %>%
  select(f.eid, f.90017.0.0, f.90016.0.0) %>%
  filter(!is.na(f.90016.0.0))

calib2 <- safe_read_csv(file_acc_weartime) %>%
  select(f.eid, f.90015.0.0) %>%
  filter(!is.na(f.90015.0.0))

calib3 <- safe_read_csv(file_acc_rawstats) %>%
  select(f.eid, f.90002.0.0, f.90180.0.0)

p4night <- safe_read_csv(file_part4_night)
p5person <- safe_read_csv(file_part5_person)
slq <- safe_read_csv(file_sleep_question)

# ---- 3. Accelerometer QC ----------------------------------------------------

calib_qc <- calib1 %>%
  inner_join(calib2, by = "f.eid") %>%
  inner_join(calib3, by = "f.eid") %>%
  filter(
    f.90016.0.0 >= 1,
    f.90017.0.0 >= 1,
    f.90015.0.0 >= 1,
    is.na(f.90002.0.0) | f.90002.0.0 == 0,
    is.na(f.90180.0.0) | f.90180.0.0 == 0
  )

message("Accelerometer QC retained N = ", nrow(calib_qc))

# ---- 4. Derive SRI and daytime nap ------------------------------------------

check_required_columns(
  p4night,
  c("ID", "SleepRegularityIndex", "duration_sib_wakinghours_atleast15min"),
  "p4night"
)

sri <- p4night %>%
  mutate(SleepRegularityIndex = abs(SleepRegularityIndex)) %>%
  group_by(ID) %>%
  summarise(
    daycount = sum(!is.na(SleepRegularityIndex)),
    sri = ifelse(daycount > 0, mean(SleepRegularityIndex, na.rm = TRUE), NA_real_),
    nap_daycount = sum(!is.na(duration_sib_wakinghours_atleast15min)),
    nap = ifelse(
      nap_daycount > 0,
      mean(duration_sib_wakinghours_atleast15min, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

# ---- 5. Derive accelerometer-based sleep timing metrics ----------------------

check_required_columns(
  p5person,
  c(
    "ID", "sleeponset_wei", "dur_spt_min_wei", "sleep_efficiency_wei",
    "dur_spt_min_WE", "dur_spt_min_WD", "sleeponset_WD", "sleeponset_WE",
    "M5_mean_peakLUX_wei", "calendar_date"
  ),
  "p5person"
)

p5part <- p5person %>%
  select(
    ID, sleeponset_wei, dur_spt_min_wei, sleep_efficiency_wei,
    M5_mean_peakLUX_wei, dur_spt_min_WE, dur_spt_min_WD,
    sleeponset_WD, sleeponset_WE, calendar_date
  ) %>%
  mutate(
    MSF = sleeponset_WE + (dur_spt_min_WE / 60) * 0.5,
    MSW = sleeponset_WD + (dur_spt_min_WD / 60) * 0.5,
    SJL = MSF - MSW,
    SJLr = abs(SJL),
    SD = dur_spt_min_wei / 60
  )

acc_sleep <- p5part %>%
  select(
    ID, SJLr, SD, sleeponset_wei, sleep_efficiency_wei,
    M5_mean_peakLUX_wei, calendar_date
  ) %>%
  inner_join(
    sri %>% select(ID, daycount, sri, nap),
    by = "ID"
  ) %>%
  mutate(
    sri2 = ifelse(sri >= quantile(sri, 0.25, na.rm = TRUE), 1, 0),
    sjl2 = ifelse(SJLr <= 1, 1, 0),
    so2  = ifelse(sleeponset_wei >= 22 & sleeponset_wei <= 24, 1, 0),
    se2  = ifelse(sleep_efficiency_wei >= 0.85, 1, 0),
    sd2  = ifelse(SD >= 7 & SD <= 9, 1, 0),
    nap2 = case_when(
      is.na(nap) ~ NA_real_,
      nap <= 1 ~ 1,
      nap > 1 ~ 0
    ),
    f.eid = stringr::str_extract(ID, "^[0-9]+(?=_)") %>% as.numeric()
  )

# ---- 6. Derive questionnaire-based sleep variables --------------------------

check_required_columns(
  slq,
  c("f.eid", "f.1220.0.0", "f.1210.0.0", "f.1200.0.0", "f.1180.0.0", "f.1160.0.0"),
  "slq"
)

slq_part <- slq %>%
  select(f.eid, f.1220.0.0, f.1210.0.0, f.1200.0.0, f.1180.0.0, f.1160.0.0) %>%
  mutate(
    across(c(f.1220.0.0, f.1210.0.0, f.1200.0.0, f.1180.0.0, f.1160.0.0),
           ~ replace(., . < 0, NA))
  ) %>%
  mutate(
    inso = case_when(
      f.1200.0.0 == 3 ~ 0,
      is.na(f.1200.0.0) ~ NA_real_,
      TRUE ~ 1
    ),
    snore = case_when(
      f.1210.0.0 == 1 ~ 0,
      f.1210.0.0 == 2 ~ 1,
      is.na(f.1210.0.0) ~ NA_real_
    ),
    doze = case_when(
      f.1220.0.0 == 2 ~ 0,
      f.1220.0.0 %in% c(0, 1) ~ 1,
      is.na(f.1220.0.0) ~ NA_real_
    ),
    prefer = case_when(
      f.1180.0.0 %in% c(2, 3) ~ 1,
      f.1180.0.0 %in% c(1, 4) ~ 0,
      is.na(f.1180.0.0) ~ NA_real_,
      TRUE ~ NA_real_
    ),
    sd2_sr = case_when(
      f.1160.0.0 >= 7 & f.1160.0.0 <= 9 ~ 1,
      is.na(f.1160.0.0) ~ NA_real_,
      TRUE ~ 0
    )
  ) %>%
  select(f.eid, inso, snore, doze, prefer, sd2_sr, f.1180.0.0)

# ---- 7. Derive epidemiological jetlag ---------------------------------------

midsleep <- p5part %>%
  select(ID, MSF, MSW) %>%
  mutate(
    f.eid = stringr::str_extract(ID, "^[0-9]+(?=_)") %>% as.numeric(),
    Midsleep = case_when(
      is.na(MSF) & is.na(MSW) ~ NA_real_,
      is.na(MSF) ~ MSW,
      is.na(MSW) ~ MSF,
      TRUE ~ MSF * 2 / 7 + MSW * 5 / 7
    )
  ) %>%
  select(f.eid, Midsleep)

epijl <- slq_part %>%
  select(f.eid, f.1180.0.0) %>%
  inner_join(midsleep, by = "f.eid")

epijl_model <- lm(Midsleep ~ f.1180.0.0, data = epijl, na.action = na.exclude)

epijl <- epijl %>%
  mutate(
    residuals_Midsleep = residuals(epijl_model),
    epijetlag = case_when(
      is.na(residuals_Midsleep) ~ NA_real_,
      abs(residuals_Midsleep) < 1 ~ 1,
      abs(residuals_Midsleep) >= 1 ~ 0
    )
  ) %>%
  select(f.eid, Midsleep, residuals_Midsleep, epijetlag) %>%
  distinct(f.eid, .keep_all = TRUE)

# ---- 8. Construct MSS and HSS ------------------------------------------------

mss_data <- acc_sleep %>%
  inner_join(slq_part, by = "f.eid") %>%
  left_join(epijl, by = "f.eid") %>%
  mutate(
    mss = inso + doze + snore + sri2 + sjl2 + so2 + se2 + sd2 + nap2 + epijetlag,
    hss = inso + doze + snore + sd2_sr + prefer,
    mss3 = case_when(
      mss < 4 ~ 0,
      mss >= 4 & mss < 8 ~ 1,
      mss >= 8 ~ 2,
      TRUE ~ NA_real_
    ),
    hss3 = case_when(
      hss < 2 ~ 0,
      hss >= 2 & hss < 4 ~ 1,
      hss >= 4 ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  inner_join(calib_qc, by = "f.eid") %>%
  select(
    -ID,
    -f.90002.0.0,
    -f.90180.0.0,
    -f.90017.0.0,
    -f.90016.0.0,
    -f.90015.0.0
  )

message("Rows after MSS/HSS construction and QC merge: ", nrow(mss_data))

# ---- 9. Merge with downstream covariate tables -------------------------------

# IMPORTANT:
# The objects below are placeholders for derived covariate tables generated from
# authorized source data. They are not included in this repository.
#
# Expected examples:
#   health, alcohol, anthro, demo, diet, edu, employ, ethnic, fh, light,
#   noise, PA, region, smoke, stimulus, bloodbiomarker, bp
#
# Load or create these objects before running the merge step below.

required_objects <- c(
  "health", "alcohol", "anthro", "demo", "diet", "edu", "employ",
  "ethnic", "fh", "light", "noise", "PA", "region", "smoke",
  "stimulus", "bloodbiomarker", "bp"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "The following required derived objects are missing: ",
    paste(missing_objects, collapse = ", "),
    "\nPlease load or generate them before running the merge step.",
    call. = FALSE
  )
}

df_list <- list(
  mss_data, health, alcohol, anthro, demo, diet, edu, employ,
  ethnic, fh, light, noise, PA, region, smoke, stimulus,
  bloodbiomarker, bp
)

invisible(lapply(seq_along(df_list), function(i) {
  report_duplicates(df_list[[i]], id_col = "f.eid", object_name = paste0("df_list[[", i, "]]"))
}))

df_list <- lapply(df_list, deduplicate_by_id, id_col = "f.eid")

final <- purrr::reduce(df_list, dplyr::inner_join, by = "f.eid")

message("Rows after merging all covariate tables: ", nrow(final))

# ---- 10. Date processing -----------------------------------------------------

if ("calendar_date" %in% names(final)) {
  final <- final %>% mutate(calendar_date = as.Date(calendar_date))
} else if ("calendar_date.x" %in% names(final)) {
  final <- final %>% mutate(calendar_date = as.Date(calendar_date.x))
} else if ("calendar_date.y" %in% names(final)) {
  final <- final %>% mutate(calendar_date = as.Date(calendar_date.y))
} else {
  stop("No calendar_date column found in final dataset.", call. = FALSE)
}

date_vars <- c(
  "ihd_date", "af_date", "is_date", "hf_date", "t2d_date", "cmd_date",
  "obese_date", "hp_date", "lipid_date", "dementia_ad_date", "death_date",
  "depressive_date", "anxiety_date", "conduction_date", "cancer_date"
)

date_vars_present <- intersect(date_vars, names(final))
final <- final %>%
  mutate(across(all_of(date_vars_present), as.Date))

# ---- 11. Follow-up time ------------------------------------------------------

make_followup_time <- function(event_date, baseline_date) {
  ifelse(
    !is.na(event_date),
    as.numeric(event_date - baseline_date),
    NA_real_
  )
}

time_map <- list(
  t_ihd = "ihd_date",
  t_af = "af_date",
  t_is = "is_date",
  t_hf = "hf_date",
  t_t2d = "t2d_date",
  t_cmd = "cmd_date",
  t_obese = "obese_date",
  t_hp = "hp_date",
  t_lipid = "lipid_date",
  t_dementia_ad = "dementia_ad_date",
  t_death = "death_date",
  t_depressive = "depressive_date",
  t_anxiety = "anxiety_date",
  t_conduction = "conduction_date",
  t_cancer = "cancer_date"
)

for (new_var in names(time_map)) {
  date_var <- time_map[[new_var]]
  if (date_var %in% names(final)) {
    final[[new_var]] <- make_followup_time(final[[date_var]], final$calendar_date)
  }
}

# ---- 12. Build analysis dataset ---------------------------------------------

model_vars <- c(
  "f.eid", "calendar_date",
  "t_ihd", "t_af", "t_is", "t_hf", "t_hp", "t_t2d", "t_obese", "t_lipid",
  "t_dementia_ad", "t_depressive", "t_anxiety", "t_cmd", "t_death",
  "t_conduction", "t_cancer",
  "cmd", "t2d", "af", "hf", "ihd", "is", "hp", "lipid", "obese",
  "dementia_ad", "depressive", "anxiety", "death", "conduction", "cancer_all",
  "SJLr", "SD", "sleeponset_wei", "sleep_efficiency_wei", "sri", "nap",
  "age", "sex", "ethnic", "urban", "edu", "employ", "tdi3", "light3",
  "noise3", "ms", "smoke", "drink", "coffee3", "tea3", "mvpa3", "bmi",
  "ldl", "glucose", "sys",
  "fh_cmd", "fh_dm", "fh_heart", "fh_stroke", "fh_depression", "fh_ad",
  "fh_hp", "fh_cancer",
  "inso", "doze", "snore", "prefer", "sri2", "sjl2", "so2", "se2", "sd2",
  "sd2_sr", "nap2", "epijetlag", "Midsleep", "residuals_Midsleep",
  "mss", "mss3", "hss", "hss3"
)

model_vars <- intersect(model_vars, names(final))
final_sub <- final[, model_vars]

safe_write_csv(final_sub, file.path(dir_derived, "final_sub.csv"))

# ---- 13. Complete-case dataset and imputation --------------------------------

selected_sleep_vars <- c(
  "inso", "doze", "snore", "prefer",
  "sri2", "sjl2", "so2", "se2", "sd2", "nap2", "epijetlag", "mss"
)
selected_sleep_vars <- intersect(selected_sleep_vars, names(final_sub))

final_complete <- final_sub %>%
  filter(complete.cases(across(all_of(selected_sleep_vars))))

if ("SD" %in% names(final_complete)) {
  final_complete <- final_complete %>%
    filter(SD >= 4 & SD <= 12)
}

message("Final complete-case sample size: ", nrow(final_complete))

safe_write_csv(final_complete, file.path(dir_derived, "final_complete.csv"))

# Multiple imputation
imputed_obj <- mice::mice(
  final_complete,
  m = 1,
  maxit = 10,
  method = "rf",
  seed = 20250326,
  printFlag = FALSE
)

final_imp <- mice::complete(imputed_obj, action = 1)

if ("cancer_all" %in% names(final_imp) && !"cancer" %in% names(final_imp)) {
  final_imp <- final_imp %>% rename(cancer = cancer_all)
}

final_imp <- final_imp %>%
  rename(
    insomnia = inso,
    preference = prefer,
    regularity = sri2,
    socialjetlag = sjl2,
    onset = so2,
    efficiency = se2,
    duration = sd2,
    duration_report = sd2_sr,
    daytime_nap = nap2
  )

safe_write_csv(final_imp, file.path(dir_derived, "final_imp.csv"))

message("Dataset construction complete.")