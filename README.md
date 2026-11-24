# hipBMD-omics-incremental-prediction

Analysis code for:

> **Beyond association: polygenic and gut microbial contributions to hip bone mineral density in a multi-cohort evaluation.**
> Gonzalez-Ramirez M.I., *et al.* (submitted to *Advanced Science*).

The study asks whether polygenic and gut microbial profiles add predictive information on total hip bone
mineral density (BMD) beyond clinical variables. Weights, thresholds and models are fixed in one set of
participants and evaluated in others, with relatives kept apart across every analysis boundary.

## Analysis sets

| Set | n | Role |
|---|---|---|
| GWAS sample | 4,379 | LOS participants without metagenomes; derives the LOS polygenic weights |
| Modeling discovery set | 1,430 | Fits PRS thresholds, the metagenomic risk score (MRS) and all models |
| Internal holdout | 351 | Evaluates the frozen models (same sequencing platform) |
| LC subset | 495 | LOS participants whose metagenomes were sequenced on a second platform |
| Guangzhou cohort | 492 | External cohort (Chinese women, different DXA scanner) |

External polygenic weights: UK Biobank total femur BMD summary statistics (N = 28,938; NGDC Genome Variation Map, GVP000082).

## Clinical base model

LOS sets: age, sex, weight, height, recorded race (African American, White and Asian indicators), sun exposure,
regular exercise, smoking, and consumption of alcohol, cheese, milk and yogurt, plus the first two genetic
principal components. Body mass index is not a predictor. Guangzhou: a reduced model trained in the modeling
discovery set with age, sex, weight, height, an Asian race indicator, alcohol consumption and regular exercise.
Body mass index appears only as a covariate in the microbiome adjustment analyses (scripts 15-17).

## Repository layout

```
config/paths.R          project root and sample-label helper (run every script from the repository root)
plink/                  genotype steps (PLINK 2 / PLINK 1.9)
R/                      analysis scripts, numbered in run order
R/figures/              figure scripts (main text and supporting information)
python/                 Bacteroides vulgatus lookup, power summary and three supporting figures
external/PLS-PBs/       principal-balance functions (not included; see below)
```

## Run order

| Step | Script(s) | What it does |
|---|---|---|
| 1 | `R/01_cohort_census.R` | Assigns participants to the GWAS pool or the evaluation set |
| 2 | `plink/01_king.sh`, `R/02_kinship_firewall.R` | KING kinship; removes relatives across boundaries; family-aware discovery/holdout split |
| 3 | `R/03_gwas_inputs.R`, `plink/02_gwas.sh` | Independent LOS GWAS |
| 4 | `R/04_sample_sizes.R` | Participant-flow cascade with consistency checks |
| 5 | `plink/03_prs_los.sh`, `plink/04_prs_external.sh` | Clumping-and-thresholding PRS (LOS and external weights) |
| 6 | `R/05a_prs_threshold_los.R`, `R/05b_prs_threshold_external.R` | Threshold selection in the discovery set (5 x 5-fold CV) |
| 7 | `R/06_prs_permutation.R` | Fixed-score permutation reference |
| 8 | `R/07_mrs_crossfit_frozen.R`, `R/08_mrs_threshold.R` | Compositional MRS (principal balances): cross-fitted, two supervision schemes, threshold curve, frozen score |
| 9 | `plink/05_prs_guangzhou.sh`, `plink/06_prs_restricted.sh` | Frozen PRS in Guangzhou; LOS scores restricted to transferable and Guangzhou variants |
| 10 | `R/09_models_internal_holdout.R`, `R/10_bootstrap_ci_holdout.R`, `R/11_discovery_cv_increments.R` | Elastic Net and XGBoost, seven nested specifications; holdout CIs; discovery contrasts |
| 11 | `R/12_transfer_lc_subset.R`, `R/13_transfer_guangzhou.R` | Transfer of the frozen models |
| 12 | `R/14_diversity_and_species.R` ... `R/18_zero_replacement_sensitivity.R` | Microbiome associations at three levels and sensitivity analyses |
| 13 | `R/19_model_interpretation.R`, `R/20_calibration_metrics.R` | Coefficients, permutation importance, SHAP; calibration and recalibrated agreement |
| 14 | `R/21_table1_cohorts.R`, `R/22_table2_results.R` | Tables 1 and 2 (and the XGBoost table) |
| 15 | `python/23_bvulgatus_targeted_lookup.py` | Targeted lookup of *Bacteroides vulgatus* (lumbar spine BMD used only here) |
| 16 | `R/24`-`R/26`, `python/27` | External PRS at alternative thresholds (set `SENS_PRSQ`, `SENS_OUT`; script 25 adds bootstrap CIs) |
| 17 | `R/28_ten_pcs_holdout.R`, `R/29_ten_pcs_lc.R` (+ `R/25` with `SENS_OUT`) | Ten instead of two genetic PCs |
| 18 | `R/30_reduced_model_holdout.R`, `R/31_species_elastic_net_and_mrs16.R`, `R/32_variant_restriction_cost.R`, `R/33_lookup_li2025.R`, `R/35_balance_loadings.R` | Remaining sensitivity analyses, the Li et al. (2025) lookup and the composition of balances 1 and 136 |
| 19 | `python/34_power_and_calibration_slopes.py` | Minimum detectable increments; calibration slopes across bootstrap resamples |
| 20 | `R/figures/*.R`, `python/fig*.py` | Figures |

## Where each result comes from

| Paper item | Script(s) |
|---|---|
| Figure 1 | `R/figures/fig1_study_design.R` (inputs from steps 1-5) |
| Figure 2 | `R/05a`, `R/05b`, `R/06`; `R/figures/fig2_prs_thresholds.R` |
| Figure 3 | `R/09`-`R/13`; `R/figures/fig3_prediction_performance.R` |
| Figure 4 | `R/14`, `R/15`, `R/16`; `R/figures/fig4_microbiome.R` |
| Figure 5 | `R/20`; `R/figures/fig5_calibration.R` |
| Table 1 | `R/21_table1_cohorts.R` |
| Table 2 | `R/09`-`R/13`, `R/20`; `R/22_table2_results.R` |
| Figure S1 | `python/figS1_genomic_qc.py` (summary files exported from the PLINK QC and GWAS runs) |
| Figure S2 | `python/figS2_cohort_distributions.py` |
| Figure S3 | `python/figS3_relatedness.py` (PI_HAT and KING summaries from the QC runs) |
| Figure S4 | `R/figures/figS4_prs_deciles.R` |
| Figure S5 | `R/09`-`R/13`; `R/figures/figS5_rmse.R` |
| Figure S6 | `R/figures/figS6_microbiome_by_cohort.R` |
| Figure S7 | `R/figures/figS7_host_variables.R` |
| Figure S8 | `R/35_balance_loadings.R`; `R/figures/figS8_balance_taxa.R` |
| Figure S9 | `R/19`; `R/figures/figS9_model_interpretation.R` |
| Tables S1-S2 | Compiled manually from the GWAS-sample description and the output of `plink/02_gwas.sh` |
| Table S3 | `R/05a`, `R/05b` |
| Table S4 | `R/09`-`R/13`, `R/22` (XGBoost rows) |
| Table S5 | `R/08_mrs_threshold.R` |
| Table S6 | `R/09`, `R/10`, `R/11` |
| Table S7 | A, G: `python/34`; B-C: `R/28`, `R/29`, `R/25`; D: `R/30`; E-F: `R/31`, `R/07`; H: `R/24`-`R/27` |
| Table S8 | `R/15_balance_covariate_decomposition.R` |
| Table S9 | `R/14`, `R/18` |
| Table S10 | `R/14`, `R/16`, `R/17` |
| Table S11 | `R/20_calibration_metrics.R` |
| Table S12 | `plink/06_prs_restricted.sh`, `R/32_variant_restriction_cost.R` |
| Table S13 | `R/19_model_interpretation.R` |
| Table S14 | `R/33_lookup_li2025.R` |
| Table S15 | `python/23_bvulgatus_targeted_lookup.py` |

## Expected data layout (not included)

Individual-level data are not distributed with this repository. The scripts expect:

```
data/genotypes/    los_qc.{bed,bim,fam}, los_qc_pca.eigenvec, guangzhou_chr{1..22}.{bed,bim,fam}
data/phenotypes/   los_dxa_bmd.xlsx, los_clinical_imputed.csv, los_lc_clinical_imputed.csv, los_questionnaire.xlsx,
                   guangzhou_clinical.csv, guangzhou_dxa.xlsx (sheets Spine, Hip), guangzhou_questionnaire.csv
data/microbiome/   los_bgi_species_abundance.csv, los_lc_species_abundance.csv, guangzhou_species_abundance.csv,
                   metagenome_sample_ids.xlsx, metagenome_sample_filter.xlsx, metagenome_id_map.xlsx,
                   metagenome_id_corrections.csv (optional; columns old, new)
data/external/     ukb_total_femur_bmd_sumstats.tsv.gz
derived/qc/        summary files exported from the PLINK QC and GWAS runs, used only by Figures S1 and S3
```

`metagenome_id_map.xlsx` links participants to metagenome samples and must contain two sheets: `bgi`
(columns `BGI_1794` = participant ID and `Metagenome_ID`) and `lc` (columns `LC_495` = participant ID and
`Metagenome_ID`).

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
pls, MASS, ggplot2 and patchwork; Python 3.12 with pandas, numpy, scipy, statsmodels, matplotlib and openpyxl;
PLINK v1.9 and PLINK 2. Random seeds are fixed (2026).
