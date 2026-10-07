# Natural-course check of the parametric g-formula.
# Model-predicted cumulative incidence under the observed compositions (the natural course) is compared with the
# nonparametric Aalen-Johansen estimate (death as a competing event) at the end of each of the 10 follow-up
# intervals. Point estimates only. Seed, imputation, cut points and the primary outcome and death models are the
# same as in 10_gformula.R, so the predicted risk at 9.81 years reproduces ref_risk of the primary analysis.
#
# Usage: OUTCOME=<outcome> Rscript scripts/15_natural_course.R
# Output: results/natural_course/<outcome>.csv

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(mice)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))
out_dir <- file.path(proj, "results", "natural_course")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
o <- outcome_name()

cohort <- load_cohort(proj, exposure = "thr", limits = "range")
df <- cohort$df
max_fu <- max(df$time_to_gs)
n_groups <- ceiling(max_fu / 365)
cuts <- seq(0, max_fu, length.out = n_groups + 1)

imp_vars <- c("age_accel", "sex", "ethnic", "tdi", "edu", "bmi", "smoke", "alcohol", "diet_score",
              "CVD", "diabetes", "hypertension", "Dyslipidemia", "Cancer", "anxiety", "depression",
              "sleep_disorder", "income", "retired", "sick_disabled", "shift", "psych_meds",
              "chol_med", "hrt", "ocp", "parity", "R1", "R2", "R3", "gs", "death", "time_to_gs")

# primary model of 10_gformula.R
make_formula <- function(long) {
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
  as.formula(paste("~", rhs), env = env)
}

fit_models <- function(dt) {
  dt <- copy(dt)[, fu_time := time_to_gs]
  long <- survSplit(Surv(time_to_gs, gs) ~ ., data = as.data.frame(dt),
                    cut = cuts[-c(1, length(cuts))], episode = "timegroup",
                    start = "start", end = "end", event = "gs")
  setDT(long)
  long[, last := end == fu_time]
  long[, death := as.integer(last & death == 1)]
  f <- make_formula(long)
  list(gs = glm(update(f, gs ~ .), data = long[death == 0], family = binomial(), model = FALSE, y = FALSE),
       death = glm(update(f, death ~ .), data = long[gs == 0], family = binomial(), model = FALSE, y = FALSE))
}

# same steps as run_once() in 10_gformula.R, keeping the risk at every interval
set.seed(2026)
dt <- add_ilr(copy(df))
dt <- impute_pmm(dt, imp_vars)
dt[, row_id := .I]
models <- fit_models(dt)
n <- nrow(dt)
long <- dt[rep(seq_len(n), each = n_groups)]
long[, timegroup := rep(seq_len(n_groups), n)]
long[, h_gs := predict(models$gs, newdata = .SD, type = "response")]
long[, h_d := predict(models$death, newdata = .SD, type = "response")]
long[, risk := cumsum(h_gs * cumprod((1 - shift(h_gs, fill = 0)) * (1 - h_d))), by = row_id]
pred <- long[, .(pred_risk = mean(risk)), by = timegroup]

# observed: Aalen-Johansen cumulative incidence with death as a competing event
dt[, state := factor(fifelse(gs == 1, 1L, fifelse(death == 1, 2L, 0L)), levels = 0:2,
                     labels = c("censored", "outcome", "death"))]
aj <- survfit(Surv(time_to_gs, state) ~ 1, data = dt)
s <- summary(aj, times = cuts[-1], extend = TRUE)
k <- match("outcome", aj$states)  # pstate columns follow the order of aj$states
stopifnot(!is.na(k))
obs <- data.table(timegroup = seq_len(n_groups), obs_risk = s$pstate[, k], obs_lower = s$lower[, k],
                  obs_upper = s$upper[, k])

res <- merge(pred, obs, by = "timegroup")
res[, `:=`(outcome = o, years = cuts[-1][timegroup] / 365.25, n = n, events = sum(dt$gs == 1))]
setcolorder(res, c("outcome", "timegroup", "years"))
fwrite(res, file.path(out_dir, paste0(o, ".csv")))

# check against the primary analysis (ref_risk = natural-course risk at the last interval)
pt <- readRDS(file.path(results_root(proj), "primary_range_thr", "point_estimate.rds"))
ref <- unique(pt[cohort == "full", ref_risk])[1]
cat(sprintf("%s: predicted %.2f vs primary ref %.2f per 1,000; observed %.2f (%.2f-%.2f); max |pred-obs| %.2f per 1,000\n",
            o, 1000 * res[.N, pred_risk], 1000 * ref, 1000 * res[.N, obs_risk], 1000 * res[.N, obs_lower],
            1000 * res[.N, obs_upper], 1000 * max(abs(res$pred_risk - res$obs_risk))))
