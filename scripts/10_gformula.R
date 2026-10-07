# Isotemporal reallocation of 24-h behaviours and incident disease: parametric g-formula with pooled logistic hazard
# models for the outcome and for death (competing event), adapted from Yiallourou et al., BMC Medicine 2026;24:11
# (github.com/BeaudanBrown/coda-dementia). The outcome is chosen with the OUTCOME environment variable (default
# gallstone); the event and follow-up are stored in the columns gs / time_to_gs for every outcome.
#
# Usage: OUTCOME=<outcome> Rscript scripts/10_gformula.R [n_boot=0] [cores=6] [model] [limits=range] [exposure=thr]
#   model:    primary | s2 (+ BMI, chronic disease) | s3 (+ lipid-lowering drugs, HRT, OCP, parity)
#             | lag3 (exclude cases diagnosed in the first 3 years of follow-up)
#   limits:   range (observed sample range) | p1_99 (1st-99th percentiles)
#   exposure: thr (ML sleep + 30/100 mg acceleration thresholds) | ml (Oxford ML 4-class)
#   n_boot = 0 runs the point estimate on the full (singly imputed) data only.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(mice)
  library(parallel)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
n_boot <- if (length(args) >= 1) as.integer(args[1]) else 0L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 6L
model_name <- if (length(args) >= 3) args[3] else "primary"
limits <- if (length(args) >= 4) args[4] else "range"
exposure <- if (length(args) >= 5) args[5] else "thr"

out_dir <- file.path(results_root(proj), paste(c(model_name, limits, if (exposure != "ml") exposure), collapse = "_"))
dir.create(file.path(out_dir, "boots"), recursive = TRUE, showWarnings = FALSE)

cohort <- load_cohort(proj, exposure = exposure, limits = limits)
df <- cohort$df
comp_limits <- cohort$comp_limits
if (model_name == "lag3") {
  df <- df[!(gs == 1 & time_to_gs < 3 * 365.25)]
}

max_fu <- max(df$time_to_gs)
n_groups <- ceiling(max_fu / 365)
cuts <- seq(0, max_fu, length.out = n_groups + 1)

imp_vars <- c("age_accel", "sex", "ethnic", "tdi", "edu", "bmi", "smoke", "alcohol", "diet_score",
              "CVD", "diabetes", "hypertension", "Dyslipidemia", "Cancer", "anxiety", "depression",
              "sleep_disorder", "income", "retired", "sick_disabled", "shift", "psych_meds",
              "chol_med", "hrt", "ocp", "parity", "R1", "R2", "R3", "gs", "death", "time_to_gs")

make_formula <- function(long, model) {
  env <- new.env(parent = globalenv())
  env$rcs3 <- rcs3
  q <- function(v, p = c(0.1, 0.5, 0.9)) unname(quantile(long[[v]], p))
  env$kT <- q("timegroup", c(0.01, 0.5, 0.99))
  env$kA <- q("age_accel")
  env$kD <- q("tdi")
  env$kR1 <- q("R1"); env$kR2 <- q("R2"); env$kR3 <- q("R3")
  rhs <- "
    rcs3(R1, kR1) * rcs3(timegroup, kT) + rcs3(R2, kR2) * rcs3(timegroup, kT) + rcs3(R3, kR3) * rcs3(timegroup, kT) +
    rcs3(R1, kR1) * sex + rcs3(R2, kR2) * sex + rcs3(R3, kR3) * sex +
    rcs3(R1, kR1) * smoke + rcs3(R2, kR2) * smoke + rcs3(R3, kR3) * smoke +
    rcs3(R1, kR1) * rcs3(age_accel, kA) + rcs3(R2, kR2) * rcs3(age_accel, kA) + rcs3(R3, kR3) * rcs3(age_accel, kA) +
    rcs3(R1, kR1) * income + rcs3(R2, kR2) * income + rcs3(R3, kR3) * income +
    rcs3(R1, kR1) * retired + rcs3(R2, kR2) * retired + rcs3(R3, kR3) * retired +
    rcs3(timegroup, kT) * rcs3(age_accel, kA) + rcs3(timegroup, kT) * sex +
    rcs3(age_accel, kA) * sex +
    edu + ethnic + alcohol + diet_score + rcs3(tdi, kD) + shift + psych_meds"
  if (model == "s2") {
    env$kB <- q("bmi")
    rhs <- paste(rhs, "+ rcs3(bmi, kB) + CVD + diabetes + hypertension + Dyslipidemia + Cancer +
                       anxiety + depression + sleep_disorder + sick_disabled")
  }
  if (model == "s3") {
    # female-specific factors are coded 0 for men; sex is already in the model
    rhs <- paste(rhs, "+ chol_med + hrt + ocp + parity")
  }
  as.formula(paste("~", rhs), env = env)
}

fit_models <- function(dt, model) {
  dt <- copy(dt)[, fu_time := time_to_gs]
  long <- survSplit(Surv(time_to_gs, gs) ~ ., data = as.data.frame(dt),
                    cut = cuts[-c(1, length(cuts))], episode = "timegroup",
                    start = "start", end = "end", event = "gs")
  setDT(long)
  long[, last := end == fu_time]
  long[, death := as.integer(last & death == 1)]
  f <- make_formula(long, model)
  m_gs <- glm(update(f, gs ~ .), data = long[death == 0], family = binomial(),
              model = FALSE, y = FALSE)
  m_death <- glm(update(f, death ~ .), data = long[gs == 0], family = binomial(),
                 model = FALSE, y = FALSE)
  strip <- function(m) {
    m[c("residuals", "fitted.values", "effects", "linear.predictors", "weights",
        "prior.weights", "data")] <- NULL
    m$qr$qr <- NULL
    m
  }
  list(gs = strip(m_gs), death = strip(m_death))
}

# Cumulative incidence of gallstone by the end of follow-up, accounting for competing death
final_risk <- function(dt, models) {
  n <- nrow(dt)
  long <- dt[rep(seq_len(n), each = n_groups)]
  long[, timegroup := rep(seq_len(n_groups), n)]
  long[, h_gs := predict(models$gs, newdata = .SD, type = "response")]
  long[, h_d := predict(models$death, newdata = .SD, type = "response")]
  long[, risk := cumsum(h_gs * cumprod((1 - shift(h_gs, fill = 0)) * (1 - h_d))), by = row_id]
  long[timegroup == n_groups, risk]
}

subs <- subs_grid()

run_once <- function(dt, model) {
  dt <- add_ilr(copy(dt))
  dt <- impute_pmm(dt, imp_vars)
  dt[, row_id := .I]
  models <- fit_models(dt, model)
  ref <- final_risk(dt, models)
  in_cohort <- cohort_index(dt)
  res <- rbindlist(lapply(seq_len(nrow(subs)), function(i) {
    s <- subs[i]
    sub <- substitute(dt, s$from, s$to, s$duration / 60, comp_limits)
    risk <- final_risk(sub$data, models)
    rbindlist(lapply(names(in_cohort), function(c) {
      idx <- in_cohort[[c]]
      data.table(cohort = c, from = s$from, to = s$to, duration = s$duration,
                 ref_risk = mean(ref[idx]), sub_risk = mean(risk[idx]),
                 prop_sub = mean(sub$ok[idx]))
    }))
  }))
  res[, rr := sub_risk / ref_risk]
  res
}

# ---- point estimate ----
point_file <- file.path(out_dir, "point_estimate.rds")
if (!file.exists(point_file)) {
  set.seed(2026)
  t0 <- Sys.time()
  point <- run_once(df, model_name)
  saveRDS(point, point_file)
  message(sprintf("point estimate done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
point <- readRDS(point_file)

# ---- bootstrap ----
if (n_boot > 0) {
  todo <- setdiff(seq_len(n_boot), as.integer(sub("\\.rds$", "", list.files(file.path(out_dir, "boots")))))
  invisible(mclapply(todo, function(b) {
    set.seed(10000 + b)
    bs <- df[sample(.N, .N, replace = TRUE)]
    res <- tryCatch(run_once(bs, model_name), error = function(e) NULL)
    if (!is.null(res)) saveRDS(res[, boot := b], file.path(out_dir, "boots", paste0(b, ".rds")))
    NULL
  }, mc.cores = n_cores, mc.preschedule = FALSE))
}

# ---- summarise ----
boot_files <- list.files(file.path(out_dir, "boots"), full.names = TRUE)
summary_dt <- copy(point)
if (length(boot_files) > 0) {
  boots <- rbindlist(lapply(boot_files, readRDS))
  ci <- boots[, .(lower = quantile(rr, 0.025), upper = quantile(rr, 0.975), n_boot = .N),
              by = .(cohort, from, to, duration)]
  summary_dt <- merge(summary_dt, ci, by = c("cohort", "from", "to", "duration"))
}
summary_dt[, reported := prop_sub >= intervention_threshold]
setorder(summary_dt, cohort, from, duration)
fwrite(summary_dt, file.path(out_dir, "risk_ratios.csv"))

cat(sprintf("\nN = %d, gallstone = %d, deaths = %d, time groups = %d\n",
            nrow(df), sum(df$gs), sum(df$death), n_groups))
cat("Cohort sizes:\n"); print(df[, .(n = .N, gallstone = sum(gs)), by = sleep_group])
show <- summary_dt[abs(duration) == 30 & cohort %in% c("short", "normal", "long", "full")]
print(show[, c("cohort", "from", "duration", "ref_risk", "rr",
               intersect(c("lower", "upper"), names(show)), "prop_sub"), with = FALSE], digits = 3)
