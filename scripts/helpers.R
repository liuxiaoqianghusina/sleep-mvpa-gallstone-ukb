# Shared helpers for the 24-h composition analyses.

behaviours <- c("sleep", "inactivity", "light", "mvpa")
short_sleep <- 7
long_sleep <- 9
intervention_threshold <- 0.75

# Outcome selected with the OUTCOME environment variable (default gallstone, the original analysis).
outcome_name <- function() Sys.getenv("OUTCOME", "gallstone")

# gallstone results stay in results/; other outcomes go to results/outcomes/<outcome>/
results_root <- function(proj) {
  o <- outcome_name()
  if (o == "gallstone") file.path(proj, "results") else file.path(proj, "results", "outcomes", o)
}

# Load the analysis dataset and prepare the composition (threshold-based by default).
load_cohort <- function(proj, exposure = "thr", limits = "range") {
  df <- fread(file.path(proj, "data", paste0("analysis_", outcome_name(), ".csv")))
  setnames(df, c("Sex", "Ethnic", "Townsend.deprivation.index.at.recruitment", "Education",
                 "Body.mass.index..BMI....Instance.0"),
           c("sex", "ethnic", "tdi", "edu", "bmi"))
  factor_vars <- c("sex", "ethnic", "edu", "smoke", "alcohol",
                   "CVD", "diabetes", "hypertension", "Dyslipidemia", "Cancer",
                   "anxiety", "depression", "sleep_disorder",
                   "retired", "sick_disabled", "shift", "psych_meds", "chol_med", "hrt", "ocp")
  df[, (factor_vars) := lapply(.SD, factor), .SDcols = factor_vars]
  if (exposure == "thr") {
    df[, `:=`(inactivity = inactivity_thr, light = light_thr, mvpa = mvpa_thr)]
  }
  min_mvpa <- min(df$mvpa[df$mvpa > 0])
  df[mvpa == 0, mvpa := 0.65 * min_mvpa]
  df[, (behaviours) := lapply(.SD, function(x) x / rowSums(.SD) * 24), .SDcols = behaviours]
  df[, sleep_group := fcase(sleep < short_sleep, "short", sleep <= long_sleep, "normal", default = "long")]
  # sensitivity groupings: sample tertiles, and NSF age-specific range (7-8 h at age >= 65)
  tert <- quantile(df$sleep, c(1 / 3, 2 / 3))
  df[, sleep_tertile := fcase(sleep < tert[1], "t1", sleep < tert[2], "t2", default = "t3")]
  upper <- ifelse(df$age_accel >= 65, 8, long_sleep)
  df[, sleep_group_age := fcase(sleep < short_sleep, "short", sleep <= upper, "normal", default = "long")]
  df[, eid := as.integer(eid)]
  probs <- if (limits == "range") c(0, 1) else c(0.01, 0.99)
  comp_limits <- lapply(behaviours, function(b) quantile(df[[b]], probs))
  names(comp_limits) <- behaviours
  list(df = df, comp_limits = comp_limits)
}

# ILR with the original sequential binary partition:
# (sleep, inactivity | light, mvpa), (sleep | inactivity), (light | mvpa)
add_ilr <- function(dt) {
  dt[, R1 := log(sqrt(sleep * inactivity) / sqrt(light * mvpa))]
  dt[, R2 := sqrt(1 / 2) * log(sleep / inactivity)]
  dt[, R3 := sqrt(1 / 2) * log(light / mvpa)]
  dt
}

# Harrell restricted cubic spline with 3 knots (linear + 1 non-linear term)
rcs3 <- function(x, k) {
  p3 <- function(u) pmax(u, 0)^3
  nl <- (p3(x - k[1]) - p3(x - k[2]) * (k[3] - k[1]) / (k[3] - k[2]) +
           p3(x - k[3]) * (k[2] - k[1]) / (k[3] - k[2])) / (k[3] - k[1])^2
  m <- cbind(x, nl)
  colnames(m) <- c("", "nl")
  m
}

# Single imputation by predictive mean matching (one dataset per bootstrap sample, as in the original)
impute_pmm <- function(dt, imp_vars) {
  x <- as.data.frame(dt[, imp_vars, with = FALSE])
  miss <- names(which(colSums(is.na(x)) > 0))
  if (length(miss) == 0) return(copy(dt))
  meth <- make.method(x)
  meth[meth != ""] <- "pmm"
  imp <- mice(x, m = 1, maxit = 10, method = meth, printFlag = FALSE)
  out <- copy(dt)
  done <- as.data.table(complete(imp, 1))
  out[, (miss) := done[, miss, with = FALSE]]
  out
}

# Move `hours` from `from` to `to`; only applied when both stay within comp_limits.
substitute <- function(dt, from, to, hours, comp_limits) {
  new_from <- dt[[from]] - hours
  new_to <- dt[[to]] + hours
  ok <- new_from >= comp_limits[[from]][1] & new_from <= comp_limits[[from]][2] &
    new_to >= comp_limits[[to]][1] & new_to <= comp_limits[[to]][2]
  sub <- copy(dt)
  set(sub, j = from, value = ifelse(ok, new_from, dt[[from]]))
  set(sub, j = to, value = ifelse(ok, new_to, dt[[to]]))
  add_ilr(sub)
  list(data = sub, ok = ok)
}

subs_grid <- function() {
  s <- CJ(from = c("mvpa", "light", "inactivity"), duration = c(-60, -45, -30, -15, 15, 30, 45, 60))
  s[, to := "sleep"]
  s
}

cohort_names <- c(outer(c("short", "normal", "long", "full"), c("", "_female", "_male"), paste0))
extra_cohorts <- list(
  tertile1 = c("sleep_tertile", "t1"), tertile2 = c("sleep_tertile", "t2"), tertile3 = c("sleep_tertile", "t3"),
  agespec_short = c("sleep_group_age", "short"), agespec_normal = c("sleep_group_age", "normal"),
  agespec_long = c("sleep_group_age", "long")
)

# Named list of logical vectors: sleep group x sex, plus the sensitivity groupings.
cohort_index <- function(dt) {
  idx <- lapply(cohort_names, function(c) {
    parts <- strsplit(c, "_")[[1]]
    i <- if (parts[1] == "full") rep(TRUE, nrow(dt)) else dt$sleep_group == parts[1]
    if (length(parts) == 2) i <- i & dt$sex == ifelse(parts[2] == "female", 0, 1)
    i
  })
  names(idx) <- cohort_names
  c(idx, lapply(extra_cohorts, function(e) dt[[e[1]]] == e[2]))
}

project_dir <- function() {
  normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
}
