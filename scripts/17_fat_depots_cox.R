# Abdominal MRI fat depots and disease occurring after the scan: cause-specific Cox models, follow-up from the
# imaging visit. MRI fields were extracted for the gallstone cohort only (07_extract_mri.R), so analyses of other
# outcomes are restricted to participants who are also in the gallstone cohort.
#
# Usage: OUTCOME=<outcome> Rscript scripts/17_fat_depots_cox.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))
source(file.path(proj, "scripts", "mri_core.R"))

df <- load_mri_sample(proj)$df
df[, end_date := as.IDate(index_date) + time_to_gs]
df <- df[end_date > mri_date]  # still at risk and gallstone-free at the scan
df[, `:=`(time = as.numeric(end_date - mri_date), event = gs)]

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
df[, `:=`(vat_z = z(vat), asat_z = z(asat), pdff_z = z(log_pdff), bmi_z = z(bmi_mri))]

base <- "age_mri + sex + ethnic + tdi + edu + smoke + alcohol + diet_score + factor(income) +
         retired + psych_meds + height + strata(mri_centre)"
fits <- list(
  "VAT (per SD)"                    = "vat_z",
  "ASAT (per SD)"                   = "asat_z",
  "Liver PDFF, log (per SD)"        = "pdff_z",
  "VAT, mutually adjusted for ASAT" = "vat_z + asat_z",
  "VAT, adjusted for BMI"           = "vat_z + bmi_z",
  "Liver PDFF, adjusted for BMI"    = "pdff_z + bmi_z"
)
res <- rbindlist(lapply(names(fits), function(lbl) {
  f <- as.formula(paste("Surv(time, event) ~", fits[[lbl]], "+", base))
  m <- coxph(f, data = df)
  term <- strsplit(fits[[lbl]], " \\+ ")[[1]][1]
  ci <- exp(confint(m)[term, ])
  data.table(model = lbl, n = m$n, events = m$nevent,
             hr = exp(coef(m)[term]), lower = ci[1], upper = ci[2],
             p = summary(m)$coefficients[term, "Pr(>|z|)"])
}))
dir.create(results_root(proj), recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(results_root(proj), if (outcome_name() == "gallstone") "fat_depots_gallstone_cox.csv"
                      else "fat_depots_outcome_cox.csv"))
cat(sprintf("At risk at imaging: %d, gallstone after scan: %d, median follow-up %.2f y\n",
            nrow(df), sum(df$event), median(df$time) / 365.25))
print(res, digits = 3)
