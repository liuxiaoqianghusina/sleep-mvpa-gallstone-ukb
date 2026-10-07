# Isotemporal substitution of 24-h behaviours and abdominal MRI phenotypes
# (visceral fat, subcutaneous fat, liver fat), mirroring the brain-MRI analysis of
# Yiallourou et al. 2025: linear models + g-computation of mean differences.
#
# Usage: Rscript scripts/16_mri_gcomputation.R [n_boot=0] [cores=4] [model=primary|bmi] [subset=all|postindex]
#   bmi = additionally adjusts for BMI at imaging (sensitivity; BMI is on the causal path)
#   postindex = only participants imaged on or after the end of accelerometer wear (sensitivity analysis)

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(parallel)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
n_boot <- if (length(args) >= 1) as.integer(args[1]) else 0L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
model_name <- if (length(args) >= 3) args[3] else "primary"
subset_name <- if (length(args) >= 4) args[4] else "all"
out_dir <- file.path(proj, "results", paste0("mri_", model_name, if (subset_name != "all") paste0("_", subset_name)))
dir.create(file.path(out_dir, "boots"), recursive = TRUE, showWarnings = FALSE)

source(file.path(proj, "scripts", "mri_core.R"))
mri_sample <- load_mri_sample(proj)
df <- mri_sample$df
if (subset_name == "postindex") {
  n_all <- nrow(df)
  df <- df[mri_date >= as.IDate(index_date)]
  cat(sprintf("post-index subset: %d of %d (excluded %d imaged before the end of accelerometer wear)\n",
              nrow(df), n_all, n_all - nrow(df)))
}
comp_limits <- mri_sample$comp_limits  # limits from the full cohort, as in the original
outcomes <- mri_outcomes
imp_vars <- mri_imp_vars
make_formula <- function(dt) make_mri_formula(dt, adjust_bmi = model_name == "bmi")

run_once <- function(dt) {
  dt <- add_ilr(copy(dt))
  dt <- impute_pmm(dt, imp_vars)
  f <- make_formula(dt)
  subs <- subs_grid()
  rbindlist(lapply(outcomes, function(y) {
    fit_dt <- dt[!is.na(get(y))]
    m <- lm(update(f, as.formula(paste(y, "~ ."))), data = fit_dt, model = FALSE)
    ref <- predict(m, newdata = fit_dt)
    in_cohort <- cohort_index(fit_dt)
    rbindlist(lapply(seq_len(nrow(subs)), function(i) {
      s <- subs[i]
      sub <- substitute(fit_dt, s$from, s$to, s$duration / 60, comp_limits)
      est <- predict(m, newdata = sub$data)
      rbindlist(lapply(names(in_cohort), function(c) {
        idx <- in_cohort[[c]]
        data.table(outcome = y, cohort = c, from = s$from, to = s$to, duration = s$duration,
                   n = sum(idx), ref_mean = mean(ref[idx]), md = mean(est[idx] - ref[idx]),
                   prop_sub = mean(sub$ok[idx]))
      }))
    }))
  }))
}

point_file <- file.path(out_dir, "point_estimate.rds")
if (!file.exists(point_file)) {
  set.seed(2026)
  saveRDS(run_once(df), point_file)
}
point <- readRDS(point_file)

if (n_boot > 0) {
  todo <- setdiff(seq_len(n_boot), as.integer(sub("\\.rds$", "", list.files(file.path(out_dir, "boots")))))
  invisible(mclapply(todo, function(b) {
    set.seed(20000 + b)
    res <- tryCatch(run_once(df[sample(.N, .N, replace = TRUE)]), error = function(e) NULL)
    if (!is.null(res)) saveRDS(res[, boot := b], file.path(out_dir, "boots", paste0(b, ".rds")))
    NULL
  }, mc.cores = n_cores, mc.preschedule = FALSE))
}

summary_dt <- copy(point)
boot_files <- list.files(file.path(out_dir, "boots"), full.names = TRUE)
if (length(boot_files) > 0) {
  boots <- rbindlist(lapply(boot_files, readRDS))
  ci <- boots[, .(lower = quantile(md, 0.025), upper = quantile(md, 0.975), n_boot = .N),
              by = .(outcome, cohort, from, to, duration)]
  summary_dt <- merge(summary_dt, ci, by = c("outcome", "cohort", "from", "to", "duration"))
}
summary_dt[, reported := prop_sub >= intervention_threshold]
setorder(summary_dt, outcome, cohort, from, duration)
fwrite(summary_dt, file.path(out_dir, "mean_differences.csv"))

cat(sprintf("\nMRI sample: %d (VAT %d, ASAT %d, liver PDFF %d); median accel->MRI gap %.2f y\n",
            nrow(df), sum(!is.na(df$vat)), sum(!is.na(df$asat)), sum(!is.na(df$log_pdff)),
            median(as.numeric(df$mri_date - as.IDate(df$index_date)) / 365.25)))
show <- summary_dt[abs(duration) == 30 & cohort %in% c("short", "normal", "long", "full")]
print(show[, c("outcome", "cohort", "from", "duration", "ref_mean", "md",
               intersect(c("lower", "upper"), names(show)), "prop_sub"), with = FALSE], digits = 3)
