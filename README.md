# Code for Association between multidimensional sleep health, its metabolic risk score, and incident chronic diseases

# MSS

This repository contains the analysis code used in the study of multidimensional sleep score (MSS), sleep-related metabolic risk, and chronic disease outcomes using data from the UK Biobank, the Rugao Longitudinal Aging Study (RuLAS), and the PCAMCI study.

## Overview

The repository documents the analytical workflow used for:

- construction of accelerometer-derived and questionnaire-derived sleep variables;
- derivation of multidimensional sleep score (MSS) and related sleep summary measures;
- main survival analyses for associations between sleep health and multiple health outcomes;
- sex-stratified analyses and interaction testing by sex;
- sensitivity analyses, including lag analyses and consistency/discordance analyses;
- construction of a sleep-related metabolic risk score (MRS) using biomarker data, LightGBM + SHAP feature importance, and ridge regression;
- generation of manuscript-style figures and supplementary plots.

The code is provided to support transparency and reproducibility of the analytical workflow.

## Repository structure

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

#Script descriptions
code/00_setup.R
#Loads required packages, defines project directories, and provides shared helper functions.
code/01_build_sleep_dataset.R
#Builds accelerometer-derived and questionnaire-derived sleep variables, applies accelerometer quality control, derives MSS and HSS, merges covariate tables, and creates analysis-ready datasets.
code/02_main_survival_analysis.R
#Runs the main Cox regression analyses, generates baseline summaries, incidence rates, and compares MSS and HSS in relation to disease outcomes.
code/03_sex_stratified_and_interaction.R
#Performs sex-stratified analyses and interaction testing by sex.
code/04_sensitivity_analyses.R
#Conducts sensitivity analyses, including lag analyses, interval adjustment, additional sleep metric analyses, and consistency/discordance analyses.
code/05_mrs_construction_and_analysis.R
#Constructs the metabolic risk score (MRS), performs biomarker importance analysis using LightGBM and SHAP, estimates biomarker weights using #ridge regression, and evaluates MRS in relation to MSS and health outcomes.
code/06_figures_and_tables.R
#Generates manuscript-style figures and supplementary plots from the processed analysis outputs.
Data availability
#The underlying datasets are not included in this repository.

UK Biobank data are available to eligible researchers upon application through the UK Biobank Access Management System:
https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access
RuLAS data were used in collaboration with the RuLAS research team and are not publicly available owing to participant privacy and ethical restrictions, but may be available from the principal investigators upon reasonable request and subject to institutional and ethical approvals.
PCAMCI data are not publicly available owing to participant privacy and ethical restrictions, but may be available from the principal investigator upon reasonable request and subject to institutional and ethical approvals.
Because of these restrictions, this repository provides analysis scripts only and does not contain raw participant-level data.


#Software environment
Analyses were developed using:

R 4.4.3
RStudio 2025.05.1+513
A session information file is provided in docs/sessionInfo.txt.

#Reproducibility notes
This repository is not intended to be a one-click fully reproducible pipeline because the underlying data are restricted and some intermediate derived datasets must be generated locally under approved data access conditions. To use these scripts:

Obtain authorized access to the source datasets.
Place raw and/or derived files in your local project structure.
Update placeholder file paths where necessary.

Some scripts assume the presence of locally derived covariate tables that are not included in the public repository.
Variable names and data structures may need to be adapted depending on local data extracts and authorized preprocessing steps.
Placeholder file names are used in the public code to avoid exposing restricted directory structures.

For questions regarding the analytical workflow, please contact the corresponding author.
