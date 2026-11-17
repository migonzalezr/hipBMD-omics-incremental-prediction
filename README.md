# hipBMD-omics-incremental-prediction

Analysis code for:

> **Beyond association: polygenic and gut microbial contributions to hip bone mineral density in a multi-cohort evaluation.**
> Gonzalez-Ramirez M.I., *et al.* (manuscript under review).

The study asks whether polygenic and gut microbial profiles add predictive information on total hip bone
mineral density (BMD) beyond age, sex, weight and other clinical variables. Weights, thresholds and models are
fixed in one set of participants and evaluated in others, with relatives kept apart across every analysis
boundary (kinship firewall).

## Analysis sets

| Set | n | Role |
|---|---|---|
| GWAS sample | 4,379 | LOS participants without metagenomes; derives the LOS polygenic weights |
| Modeling discovery set | 1,430 | Fits PRS thresholds, the metagenomic risk score (MRS) and all models |
| Internal holdout | 351 | Evaluates the frozen models (same sequencing platform) |
| LC subset | 495 | LOS participants whose metagenomes were sequenced on a second platform |
| Guangzhou cohort | 492 | External cohort (Chinese women, different DXA scanner) |

External polygenic weights: UK Biobank total femur BMD summary statistics (N = 28,938; NGDC Genome Variation Map, GVP000082).

## Repository layout

```
config/paths.R          project root (run every script from the repository root)
plink/                  genotype steps (PLINK 2 / PLINK 1.9)
R/                      analysis scripts, numbered in run order
R/figures/              main-text figures 1-4
external/PLS-PBs/       principal-balance functions (not included; see below)
```

## Run order

| Step | Script | What it does | Paper |
|---|---|---|---|
| 1 | `R/01_cohort_census.R` | Assigns participants to the GWAS pool or the evaluation set | Fig. 1a |
| 2 | `plink/01_king.sh` | KING kinship table | Methods |
| 3 | `R/02_kinship_firewall.R` | Removes relatives across boundaries; family-aware 80:20 discovery/holdout split | Fig. 1a |
| 4 | `R/03_gwas_inputs.R`, `plink/02_gwas.sh` | Independent LOS GWAS | Methods |
| 5 | `R/04_sample_sizes.R` | Participant-flow cascade with consistency checks | Fig. 1a |
| 6 | `plink/03_prs_los.sh`, `plink/04_prs_external.sh` | Clumping and thresholding PRS (LOS and external weights) | Methods |
| 7 | `R/05a_prs_threshold_los.R`, `R/05b_prs_threshold_external.R` | Threshold selection in the discovery set (5 x 5-fold CV) | Fig. 2a |
| 8 | `R/06_prs_permutation.R` | Permutation test of the selected scores | Fig. 2c |
| 9 | `R/07_mrs_crossfit_frozen.R`, `R/08_mrs_threshold.R` | Compositional MRS (principal balances), cross-fitted and frozen | Methods |
| 10 | `plink/05_prs_guangzhou.sh` | Frozen PRS applied to the Guangzhou cohort | Methods |
| 11 | `R/09_models_internal_holdout.R` | Elastic Net and XGBoost, seven nested specifications; internal holdout | Table 2 |
| 12 | `R/10_bootstrap_ci_holdout.R`, `R/11_discovery_cv_increments.R` | Bootstrap CIs; discovery-set increments | Table 2 |
| 13 | `R/12_transfer_lc_subset.R`, `R/13_transfer_guangzhou.R` | Transfer of the frozen models | Table 2, Fig. 4 |
| 14 | `R/14_diversity_and_species.R`, `R/15_balance_covariate_decomposition.R` | Microbiome associations at three levels | Fig. 3 |
| 15 | `R/16_host_variables_microbiome.R`, `R/17_adjustment_sensitivity.R`, `R/18_zero_replacement_sensitivity.R` | Host variables and sensitivity analyses | Fig. S2, Tables S9-S10 |
| 16 | `R/19_model_interpretation.R` | Permutation importance, SHAP, coefficients | Fig. S6, Table S14 |
| 17 | `R/20_calibration_metrics.R` | Calibration of the holdout and the LC subset | Fig. 4 |
| 18 | `R/21_table1_cohorts.R`, `R/22_table2_results.R` | Tables 1 and 2 | Tables 1-2 |
| 19 | `R/figures/fig1_study_design.R` ... `fig4_calibration.R` | Main-text figures | Figs. 1-4 |

## Expected data layout (not included)

Individual-level data are not distributed with this repository. The scripts expect:

```
data/genotypes/    los_qc.{bed,bim,fam}, los_qc_pca.eigenvec, guangzhou_chr{1..22}.{bed,bim,fam}
data/phenotypes/   los_hip_bmd.xlsx, los_clinical_imputed.csv, los_lc_clinical_imputed.csv, guangzhou_clinical.csv
data/microbiome/   los_bgi_species_abundance.csv, los_lc_species_abundance.csv, guangzhou_species_abundance.csv,
                   metagenome_sample_ids.xlsx, metagenome_sample_filter.xlsx, metagenome_id_map.xlsx,
                   metagenome_id_corrections.csv (optional; columns old, new)
data/external/     ukb_total_femur_bmd_sumstats.tsv.gz
```

`derived/`, `results/` and `figures/` are created by the scripts.

Data availability: LOS gut metagenomes, SRA PRJNA1015234 and PRJNA1015228; Guangzhou whole-genome sequencing,
GSA-Human HRA004900, and metagenomes, SRA PRJNA986283 and PRJNA1011937; UK Biobank GWAS summary statistics,
NGDC GVM GVP000082. LOS whole-genome sequencing data are being deposited in the NIA Aging Research Biobank.

## External functions

Principal balances were derived with the partial-least-squares functions of
<https://github.com/NesrstovaV/PLS-PBs>. Clone that repository into `external/PLS-PBs/` before running the
MRS scripts.

## Software

R 4.6.0 with data.table, readxl, glmnet 5.0, xgboost 3.2.1.1, vegan 2.7.5, zCompositions 1.6.2, compositions,
pls, MASS, ggplot2 and patchwork; PLINK v1.9 and PLINK 2. Random seeds are fixed (2026).
