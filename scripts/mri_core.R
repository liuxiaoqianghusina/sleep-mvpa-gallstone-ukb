# Abdominal MRI sample and linear-model formula, used by 16_mri_gcomputation.R, 17_fat_depots_cox.R and the mrisub
# variant of 13_sensitivity_models.R.

# Returns list(df, comp_limits) for participants with VAT or liver PDFF at the imaging visit.
load_mri_sample <- function(proj) {
  cohort <- load_cohort(proj, exposure = "thr", limits = "range")
  mri <- fread(file.path(proj, "data", "raw_mri_fields.csv"),
               colClasses = list(character = c("p23357_i2", "p23358_i2")))
  mri[, `:=`(
    mri_date = as.IDate(p53_i2),
    mri_centre = factor(p54_i2),
    height = p50_i2,
    bmi_mri = p21001_i2,
    vat = fifelse(is.na(p23357_i2) | p23357_i2 == "", p22407_i2, NA_real_),
    asat = fifelse(is.na(p23358_i2) | p23358_i2 == "", p22408_i2, NA_real_),
    log_pdff = log(p40061_i2)
  )]
  df <- merge(cohort$df, mri[, .(eid, mri_date, mri_centre, height, bmi_mri, vat, asat, log_pdff)], by = "eid")
  df <- df[!is.na(mri_date) & (!is.na(vat) | !is.na(log_pdff))]
  df[, age_mri := age_accel + as.numeric(mri_date - as.IDate(index_date)) / 365.25]
  list(df = df, comp_limits = cohort$comp_limits)
}

mri_outcomes <- c("vat", "asat", "log_pdff")

mri_imp_vars <- c("age_mri", "sex", "ethnic", "tdi", "edu", "smoke", "alcohol", "diet_score", "income",
                  "retired", "shift", "psych_meds", "height", "bmi_mri", "mri_centre", "R1", "R2", "R3")

make_mri_formula <- function(dt, adjust_bmi = FALSE) {
  env <- new.env(parent = globalenv())
  env$rcs3 <- rcs3
  q <- function(v) unname(quantile(dt[[v]], c(0.1, 0.5, 0.9), na.rm = TRUE))
  env$kA <- q("age_mri"); env$kD <- q("tdi"); env$kH <- q("height")
  env$kR1 <- q("R1"); env$kR2 <- q("R2"); env$kR3 <- q("R3")
  rhs <- "
    rcs3(R1, kR1) * income + rcs3(R2, kR2) * income + rcs3(R3, kR3) * income +
    rcs3(R1, kR1) * sex + rcs3(R2, kR2) * sex + rcs3(R3, kR3) * sex +
    rcs3(R1, kR1) * retired + rcs3(R2, kR2) * retired + rcs3(R3, kR3) * retired +
    rcs3(R1, kR1) * smoke + rcs3(R2, kR2) * smoke + rcs3(R3, kR3) * smoke +
    rcs3(R1, kR1) * rcs3(age_mri, kA) + rcs3(R2, kR2) * rcs3(age_mri, kA) + rcs3(R3, kR3) * rcs3(age_mri, kA) +
    rcs3(age_mri, kA) * sex +
    edu + ethnic + alcohol + diet_score + rcs3(tdi, kD) + shift + psych_meds +
    rcs3(height, kH) + mri_centre"
  if (adjust_bmi) {
    env$kB <- q("bmi_mri")
    rhs <- paste(rhs, "+ rcs3(bmi_mri, kB)")
  }
  as.formula(paste("~", rhs), env = env)
}
