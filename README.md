# Code for Association Between Multidimensional Sleep Health, Metabolic Risk Score, and Incident Chronic Diseases (MSS)

This repository contains the analysis code used in the study of multidimensional sleep score (MSS), sleep-related metabolic risk, and chronic disease outcomes based on data from the UK Biobank, the Rugao Longitudinal Aging Study (RuLAS), and the PCAMCI study.

---

## Overview

This repository documents the analytical workflow for:

- Construction of accelerometer-derived and questionnaire-derived sleep variables
- Derivation of multidimensional sleep score (MSS) and related sleep summary measures
- Main survival analyses for associations between sleep health and multiple health outcomes
- Sex-stratified analyses and interaction testing by sex
- Sensitivity analyses including lag analyses and consistency/discordance analyses
- Construction of a sleep-related metabolic risk score (MRS) using biomarker data, LightGBM + SHAP feature importance, and ridge regression
- Generation of manuscript-style figures and supplementary plots

The code is provided to support transparency and reproducibility of the analytical workflow.

---

## Repository Structure
MSS/
├─ README.md
├─ .gitignore
├─ code/
│  ├─ 00_setup.R
│  ├─ 01_build_sleep_dataset.R
│  ├─ 02_main_survival_analysis.R
│  ├─ 03_sex_stratified_and_interaction.R
│  ├─ 04_sensitivity_analyses.R
│  ├─ 05_mrs_construction_and_analysis.R
│  └─ 06_figures_and_tables.R
├─ docs/
│  └─ sessionInfo.txt
└─ output/
---

## Script Descriptions

- **code/00_setup.R**  
  Loads required packages, defines project directories, and provides shared helper functions.

- **code/01_build_sleep_dataset.R**  
  Builds accelerometer-derived and questionnaire-derived sleep variables, applies accelerometer quality control, derives MSS and HSS, merges covariate tables, and creates analysis-ready datasets.

- **code/02_main_survival_analysis.R**  
  Runs the main Cox regression analyses, generates baseline summaries and incidence rates, and compares MSS and HSS in relation to disease outcomes.

- **code/03_sex_stratified_and_interaction.R**  
  Performs sex-stratified analyses and interaction testing by sex.

- **code/04_sensitivity_analyses.R**  
  Conducts sensitivity analyses including lag analyses, interval adjustment, additional sleep metric analyses, and consistency/discordance analyses.

- **code/05_mrs_construction_and_analysis.R**  
  Constructs metabolic risk score (MRS), performs biomarker importance analysis using LightGBM and SHAP, estimates biomarker weights via ridge regression, and evaluates MRS in relation to MSS and health outcomes.

- **code/06_figures_and_tables.R**  
  Generates manuscript-style figures and supplementary plots from processed analysis outputs.

---

## Data Availability

The underlying datasets are **not included** in this repository.

- **UK Biobank data**: Available to eligible researchers through the [UK Biobank Access Management System](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access).

- **RuLAS data**: Used in collaboration with the RuLAS research team; not publicly available due to privacy and ethical restrictions, but may be requested from principal investigators subject to approvals.

- **PCAMCI data**: Not publicly available for privacy/ethical reasons; available upon reasonable request subject to approvals.

Because of these restrictions, this repository contains **analysis scripts only** and does **not** include raw participant-level data.

---

## Software Environment

Analyses were developed using:

- R 4.4.3  
- RStudio 2025.05.1+513

A session information file is provided at: `docs/sessionInfo.txt`

---

## Reproducibility Notes

This repository is **not** a one-click fully reproducible pipeline due to data access restrictions and local data preprocessing steps. To use these scripts, you must:

- Obtain authorized access to source datasets
- Place raw/derived files in your local project structure
- Update placeholder file paths as needed
- Provide locally derived covariate tables not included here
- Adapt variable names and data structures based on your local extracts and preprocessing
- Replace placeholder filenames avoiding exposure of restricted directories

For questions regarding the analytical workflow, please contact the corresponding author.

---

Thank you for your interest!
