#!/usr/bin/env bash
# The complete pipeline, with the settings used for the paper. Each bootstrap replicate is saved to
# results/.../boots/<b>.rds and existing replicates are skipped, so any step can be stopped and restarted.
# On 8 cores a 200-replicate g-formula run takes about 2.5 hours and the composition search about 4.5 hours per
# outcome; the whole pipeline takes several days.
#
# Usage: bash run_all.sh [step ...]
#   steps: data main sensitivity groups compositions natural_course mri heterogeneity additional figures tables
#   (no argument runs every step in this order)
set -euo pipefail
cd "$(dirname "$0")"
source config.sh
C=${N_CORES:-8}
DISEASES="gallstone peptic_ulcer gerd diverticular ibs"
mkdir -p data results logs

data() {
  python3 scripts/01_accelerometry.py
  Rscript scripts/02_extract_fields.R
  python3 scripts/03_baseline_covariates.py
  Rscript scripts/04_extract_healthuser.R
  Rscript scripts/05_define_outcomes.R
  python3 scripts/06_build_analysis_datasets.py
  Rscript scripts/07_extract_mri.R
}

# primary model for one outcome, then group comparisons, risk differences and E-values
primary() {  # $1 outcome, $2 bootstrap replicates
  OUTCOME=$1 Rscript scripts/10_gformula.R "$2" "$C" primary range thr > "logs/$1_primary.log" 2>&1
  OUTCOME=$1 Rscript scripts/11_compare_groups.R
  OUTCOME=$1 Rscript scripts/12_absolute_risk_evalues.R
}

main() {
  for o in $DISEASES; do primary "$o" 500; done
  for o in gallstone_surgery nco_cyst nco_earwax; do primary "$o" 200; done  # cholecystectomy endpoint, negative controls
}

sensitivity() {
  for o in $DISEASES; do
    for m in lag3 s2; do  # excluding cases in the first 3 years; + BMI and chronic disease
      OUTCOME=$o Rscript scripts/10_gformula.R 200 "$C" "$m" range thr > "logs/${o}_${m}.log" 2>&1
    done
    OUTCOME=$o Rscript scripts/13_sensitivity_models.R healthuser 200 "$C" > "logs/${o}_healthuser.log" 2>&1
    OUTCOME=$o Rscript scripts/11_compare_groups.R robust_healthuser
  done
  # gallstone disease: gallstone-specific factors, machine-learning exposure, four-knot splines, main-effects model,
  # MRI subsample
  Rscript scripts/10_gformula.R 200 "$C" s3 range thr > logs/gallstone_s3.log 2>&1
  Rscript scripts/10_gformula.R 200 "$C" primary range ml > logs/gallstone_ml.log 2>&1
  for v in knots4 reduced mrisub; do
    Rscript scripts/13_sensitivity_models.R "$v" 200 "$C" > "logs/gallstone_${v}.log" 2>&1
  done
  # peptic ulcer and GERD: excluding users of acid-suppressing drugs or antacids; peptic ulcer: + aspirin and NSAIDs
  for o in peptic_ulcer gerd; do
    OUTCOME=$o Rscript scripts/13_sensitivity_models.R noacid 200 "$C" > "logs/${o}_noacid.log" 2>&1
    OUTCOME=$o Rscript scripts/11_compare_groups.R robust_noacid
  done
  OUTCOME=peptic_ulcer Rscript scripts/13_sensitivity_models.R nsaid 200 "$C" > logs/peptic_ulcer_nsaid.log 2>&1
}

groups() {  # separate models within each sleep-duration group
  for o in $DISEASES; do
    OUTCOME=$o Rscript scripts/13_sensitivity_models.R strat 200 "$C" > "logs/${o}_strat.log" 2>&1
  done
}

compositions() {  # lowest-risk and highest-risk 24-hour compositions
  for o in $DISEASES; do
    OUTCOME=$o Rscript scripts/14_composition_search.R 200 "$C" > "logs/${o}_compositions.log" 2>&1
    OUTCOME=$o Rscript scripts/12_absolute_risk_evalues.R  # adds E-values for the composition contrasts
  done
}

natural_course() {
  for o in $DISEASES gallstone_surgery nco_cyst nco_earwax; do
    OUTCOME=$o Rscript scripts/15_natural_course.R > "logs/${o}_natural_course.log" 2>&1
  done
}

mri() {
  for m in primary bmi; do
    for s in all postindex; do
      Rscript scripts/16_mri_gcomputation.R 500 "$C" "$m" "$s" > "logs/mri_${m}_${s}.log" 2>&1
    done
  done
  for o in $DISEASES; do OUTCOME=$o Rscript scripts/17_fat_depots_cox.R > "logs/${o}_fat_cox.log" 2>&1; done
}

heterogeneity() {  # needs the primary analyses and the fat-depot Cox models
  Rscript scripts/18_heterogeneity.R > logs/heterogeneity.log 2>&1
  Rscript scripts/19_vat_metaregression.R > logs/vat_metaregression.log 2>&1
}

additional() {  # ratios to the negative controls, MVPA vs light activity, attenuation, 99% CIs, feasibility
  Rscript scripts/20_additional_estimates.R > logs/additional_estimates.log 2>&1
}

figures() {
  for f in fig2_data fig2_main_data fig3_data fig4_figS5_data figS6_data; do Rscript "figures/data/$f.R"; done
  for f in fig1_design fig2_main fig3 fig4_mri figS1_dag figS2_natural_course figS3_compositions \
           figS4_five_diseases figS5_bmi_gradient figS6_sleep_groups; do
    python3 "figures/$f.py"
  done
}

tables() {
  python3 tables/table1.py
  Rscript tables/supp_tables_data.R
}

steps=${*:-"data main sensitivity groups compositions natural_course mri heterogeneity additional figures tables"}
for s in $steps; do
  echo "=== $s $(date)"
  "$s"
done
echo "=== done $(date)"
