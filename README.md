# Trading sleep for moderate-to-vigorous physical activity and gallstone disease (UK Biobank)

Analysis code for:

> [Authors]. Trading sleep for moderate-to-vigorous physical activity and risk of gallstone disease: a prospective
> cohort study with disease-specificity tests and abdominal MRI. *BMC Medicine* [year; volume:article]. doi:[DOI]

The study uses wrist accelerometry from UK Biobank (2013–2015). Each participant's day is split into four parts: sleep,
inactivity, light activity and moderate-to-vigorous physical activity (MVPA). The analysis estimates how 9.8-year
cumulative risks change when 15–60 min/day is moved between sleep and each type of activity. Compositions are
expressed as isometric log-ratio (ILR) coordinates. Risks come from a parametric g-formula with pooled logistic models
for the outcome and for death, which is treated as a competing event.

- **Primary outcome:** incident gallstone disease.
- **Secondary outcomes:** peptic ulcer, gastro-oesophageal reflux disease (GERD), diverticular disease and irritable
  bowel syndrome (IBS).
- **Additional checks:** a cholecystectomy endpoint, two negative control outcomes (follicular cyst, L72; other
  disorders of the external ear, H61) and abdominal MRI analyses of visceral fat, subcutaneous fat and liver fat.

The g-formula code is adapted from Yiallourou et al., *BMC Medicine* 2026;24:11
([github.com/BeaudanBrown/coda-dementia](https://github.com/BeaudanBrown/coda-dementia)).

## Data access

UK Biobank data cannot be shared by the authors. Bona fide researchers can apply for access at
<https://www.ukbiobank.ac.uk>. This research was conducted using the UK Biobank Resource under application number
[application number].

**This repository contains no participant-level data or results derived from them.** The `.gitignore` excludes
`data/`, `results/`, the figure and table outputs and `config.sh`.

The pipeline reads four participant-level extracts and one input table for FreeUKB. Set their paths in `config.sh`.
In the extracts, column names follow the UK Biobank Research Analysis Platform (RAP) format, for example `eid`, `p31`,
`p53_i0` and `p20003_i0_a0`.

| Variable in `config.sh` | Fields used | Read by |
|---|---|---|
| `UKB_ACCELEROMETER` | 40034–40041 (Oxford machine-learning behaviours, weekday and weekend); 90002, 90010–90012, 90015–90017, 90051, 90113, 90127 (wear, quality control, acceleration fractions) | 01 |
| `UKB_POPULATION` | 31, 34, 52, 191, 21022, 22189 | 01, 02 |
| `UKB_ASSESSMENT` | recruitment visit (instance 0): 53, 21000, 6138, 21001, 20116, 20117, diet 1289–1468, 738, 6142, 826, 2734, 2784, 2814, 6153, 6177, 20003, 20004, 2178, 2188, 137, 2345, 2674, 2365, 3393, 2247; imaging visit (instance 2): 53, 54, 50, 21001, 22407, 22408, 23357, 23358, 40061 | 02, 04, 07 |
| `UKB_HEALTH_OUTCOMES` | first-occurrence dates of the outcomes and of the chronic diseases (130706–131794, listed in the scripts); 41272, 41282 (OPCS-4); 40005, 40006 (cancer registry) | 02, 05 |
| `FREEUKB_INPUT` | fields 53, 54, 190, 191, 20002, 40000, 40001, 40002, 41270 and 41280, exported with the RAP Table Exporter and named as listed by `FreeUKB::ukb_required_fields("outcome_followup")` | 05 |

The definitions and codings of all variables are given in the paper's Additional file 1 (Supplementary Methods,
Table M1). The medication code lists are in `code_lists/`.

## Software

- **R 4.2.2:** data.table, survival, mice, ks and parallel. Install with `Rscript install_packages.R`.
- **FreeUKB 0.2.2:** used only for outcome ascertainment (`scripts/05_define_outcomes.R`). It is not on CRAN; see
  `install_packages.R`.
- **Python 3:** pandas, numpy, matplotlib and Pillow. Install with `pip install -r requirements.txt`.
- **Operating system:** the analyses were run on macOS. The bootstrap uses `parallel::mclapply`, which forks R
  processes. On Windows, set `N_CORES=1`.

## Running the pipeline

```bash
cp config.example.sh config.sh    # then edit the paths
bash run_all.sh                   # every step, in order
bash run_all.sh data main         # or selected steps
```

`run_all.sh` lists the exact commands and settings used for the paper. It uses 500 bootstrap replicates for the
primary and MRI analyses and 200 for the others. Each bootstrap replicate is saved as it finishes, and replicates that
already exist are skipped, so any step can be interrupted and restarted.

On 8 cores:

- a point estimate takes about 4 minutes;
- a 200-replicate run takes about 2.5 hours;
- the composition search takes about 4.5 hours per outcome;
- the whole pipeline takes several days.

### Scripts

| Step | Script | What it does | Main output |
|---|---|---|---|
| data | `scripts/01_accelerometry.py` | device quality control, 24-hour composition, day 0 | `data/accel_base.csv` |
| | `scripts/02_extract_fields.R` | baseline and health-record fields for accelerometer participants | `data/raw_*_fields.csv` |
| | `scripts/03_baseline_covariates.py` | covariates (demographics, lifestyle, diet score, BMI, chronic diseases) | `data/baseline_covariates.csv` |
| | `scripts/04_extract_healthuser.R` | healthy-user variables; hearing-aid use | `data/healthuser_covariates.csv` |
| | `scripts/05_define_outcomes.R` | outcome ascertainment with FreeUKB (eight definitions) | `data/outcomes/<name>/out.csv` |
| | `scripts/06_build_analysis_datasets.py` | exclusions, follow-up and covariates for each outcome | `data/analysis_<outcome>.csv` |
| | `scripts/07_extract_mri.R` | abdominal MRI phenotypes (imaging visit) | `data/raw_mri_fields.csv` |
| main | `scripts/10_gformula.R` | g-formula, all reallocations, sleep-duration groups and sex; models `primary`, `lag3`, `s2`, `s3`; exposure `thr` or `ml` | `results/.../<model>_range_thr/` |
| | `scripts/11_compare_groups.R` | RRs by group and ratios of RRs between groups (paired bootstrap P) | `table_ratio_between_groups.csv` |
| | `scripts/12_absolute_risk_evalues.R` | risk differences per 1,000 and E-values | `absolute_risk_and_evalues.csv` |
| sensitivity, groups | `scripts/13_sensitivity_models.R` | variants `healthuser`, `knots4`, `reduced` (main-effects model), `mrisub`, `noacid`, `nsaid`, `strat` (separate models by sleep-duration group) | `results/.../robust_<variant>/` |
| compositions | `scripts/14_composition_search.R` | lowest-risk, typical and highest-risk 24-hour compositions (train and test halves) | `results/.../ideal_composition/` |
| natural_course | `scripts/15_natural_course.R` | predicted vs Aalen–Johansen cumulative incidence | `results/natural_course/<outcome>.csv` |
| mri | `scripts/16_mri_gcomputation.R` | reallocation and VAT, ASAT and liver PDFF (linear g-computation); `postindex` excludes participants scanned before the end of accelerometer wear | `results/mri_<model>[_postindex]/` |
| | `scripts/17_fat_depots_cox.R` | fat depots and disease after imaging (Cox) | `fat_depots_*_cox.csv` |
| heterogeneity | `scripts/18_heterogeneity.R` | Cochran's Q, pairwise ratios of RRs (Holm), adiposity gradient, meta-regression | `results/heterogeneity/` |
| | `scripts/19_vat_metaregression.R` | meta-regression on the visceral-fat gradient | `results/heterogeneity/` |
| additional | `scripts/20_additional_estimates.R` | ratios to the negative controls, MVPA vs light activity, attenuation, 99% CIs, feasibility | `results/additional_estimates/` |
| figures | `figures/data/*.R`, `figures/*.py` | source data and figures | `figures/output/` |
| tables | `tables/table1.py`, `tables/supp_tables_data.R` | Table 1 and the data for Tables S3–S5 | `tables/output/` |

`scripts/helpers.R`, `scripts/gformula_core.R` and `scripts/mri_core.R` hold shared functions.

### Conventions

- The outcome is selected with the `OUTCOME` environment variable: `gallstone` (default), `gallstone_surgery`
  (cholecystectomy endpoint), `peptic_ulcer`, `gerd`, `diverticular`, `ibs`, `nco_cyst` or `nco_earwax`.
- For every outcome, the event and the follow-up time are stored in the columns `gs` and `time_to_gs`. These names
  come from the original gallstone analysis.
- Gallstone results are written to `results/<analysis>/` and other outcomes to `results/outcomes/<outcome>/<analysis>/`.
- `primary_range_thr` is the main analysis:
  - `range` means a reallocation is applied only when every behaviour stays within its observed range (the minimum to
    maximum after trimming);
  - `thr` is the exposure used in the paper: machine-learning sleep, with waking time divided at 30 and 100 mg;
  - `primary_range` is the sensitivity analysis with all four parts taken from the machine-learning classification.
- An estimate is reported only when at least 75% of a group can take the reallocation (`reported` column).

## Paper items and the code that produces them

| Item | Code |
|---|---|
| Fig. 1 | `figures/fig1_design.py` (participant flow from `data/*flow*.csv`) |
| Fig. 2 | `scripts/10`; `figures/data/fig2_data.R`, `fig2_main_data.R`; `figures/fig2_main.py` |
| Fig. 3 | sensitivity analyses (`scripts/10`, `scripts/13`); `figures/data/fig3_data.R`; `figures/fig3.py` |
| Fig. 4 | `scripts/16`, `scripts/17`; `figures/data/fig4_figS5_data.R`; `figures/fig4_mri.py` |
| Table 1 | `tables/table1.py` |
| Table 2 | `absolute_risk_and_evalues.csv` (scripts 10 and 12), the sensitivity analyses (`figures/data/fig3_robustness.csv`), `results/additional_estimates/negative_control_ratios.csv` (script 20) |
| Fig. S1 | `figures/figS1_dag.py` |
| Fig. S2 | `scripts/15_natural_course.R`; `figures/figS2_natural_course.py` |
| Fig. S3 | `scripts/14_composition_search.R`; `figures/figS3_compositions.py` |
| Fig. S4 | `figures/data/fig2_data.R`; `figures/figS4_five_diseases.py` |
| Fig. S5 | `scripts/18_heterogeneity.R`; `figures/data/fig4_figS5_data.R`; `figures/figS5_bmi_gradient.py` |
| Fig. S6 | `scripts/10`, `scripts/13` (`strat`, `reduced`); `figures/data/figS6_data.R`; `figures/figS6_sleep_groups.py` |
| Tables S1–S2 | design and outcome definitions (`scripts/05`, `scripts/06`) |
| Tables S3–S5 | `tables/supp_tables_data.R` and `figures/data/figS6_*.csv` |
| Table S6 | `results/mri_*/mean_differences.csv` and `figures/data/fig4_fatcox.csv` |
| Table M3 | `code_lists/` |

## Reproducibility notes

- **Random seeds:**
  - point estimates use `set.seed(2026)`;
  - bootstrap replicate *b* uses `set.seed(10000 + b)` in script 10, `20000 + b` in 16, `30000 + b` in 14 and
    `50000 + b` in 13;
  - the composition search splits the cohort with `set.seed(2026)` and evaluates the test half with
    `set.seed(2027)`.

  Results do not depend on the number of cores.
- **Missing covariates:** they are imputed once by predictive mean matching (mice, 10 iterations). Imputation is
  repeated within every bootstrap replicate.
- **Outcome ascertainment:** `scripts/05_define_outcomes.R` calls `FreeUKB::ukb_outcome_custom` with first-occurrence,
  hospital inpatient and death records. It replaces FreeUKB's internal censoring-date function for that R session, so
  that follow-up ends at the earliest end date of the three sources in each region:
  - England: 31 March 2023;
  - Scotland: 31 August 2022;
  - Wales: 31 May 2022.
- **Data releases:** UK Biobank records are updated over time, so extracts made at other dates can give slightly
  different counts.
- **Checks before release:**
  - **Data preparation:** the scripts were re-run from the raw extracts. They reproduced the analysis datasets used for
    the paper: the same participants in the same order, with identical outcome dates, follow-up and covariates.
  - **Point estimates:** from these datasets, the point estimates of the primary g-formula, the healthy-user model,
    the MRI primary model and the natural-course check were identical to the published ones.
  - **Figures and tables:** from the stored bootstrap results, the figure scripts reproduced every figure
    pixel for pixel, and the table scripts reproduced Table 1 and the supplementary table data.

## Licence

MIT (see `LICENSE`). FreeUKB is distributed by its developers under its own MIT licence.

## Contact

[Corresponding author, email]
