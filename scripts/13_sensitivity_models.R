# Robustness of the sleep-group effect modification and of the MRI triangulation.
#
# Usage: OUTCOME=<outcome> Rscript scripts/13_sensitivity_models.R <variant> [n_boot=200] [cores=8]
#   strat   separate g-formula models within each sleep group (main-effects formula,
#           as the short-sleeper group has few events), so group-specific RRs no longer
#           borrow strength from the other groups
#   knots4  primary model with 4 spline knots (5/35/65/95th pct) for the ILR coordinates
#   mrisub  primary model restricted to participants with abdominal MRI
#   reduced pooled model with the same main-effects formula as `strat`, to separate the
#           effect of stratifying from the effect of simplifying the formula
#   nsaid   primary model + baseline aspirin and NSAID use (peptic ulcer)
#   noacid  primary model excluding baseline users of acid suppressants or antacids
#           (GERD, peptic ulcer: people already treated for reflux / dyspepsia)
#   healthuser  primary model + self-rated health, long-standing illness, number of
#           medications and cancer screening uptake (healthy-user / health-seeking bias)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(mice)
  library(parallel)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))
source(file.path(proj, "scripts", "gformula_core.R"))

args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]
n_boot <- if (length(args) >= 2) as.integer(args[2]) else 200L
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 8L
stopifnot(variant %in% c("strat", "knots4", "mrisub", "reduced", "nsaid", "noacid", "healthuser"))
out_dir <- file.path(results_root(proj), paste0("robust_", variant))
dir.create(file.path(out_dir, "boots"), recursive = TRUE, showWarnings = FALSE)

cohort <- load_cohort(proj, exposure = "thr", limits = "range")
df <- cohort$df
comp_limits <- cohort$comp_limits
if (variant == "mrisub") {
  source(file.path(proj, "scripts", "mri_core.R"))
  df <- df[eid %in% load_mri_sample(proj)$df$eid]
}
if (variant == "noacid") df <- df[acid_rx == 0 & antacid == 0]
imp_vars <- imp_vars_gs
if (variant == "healthuser") {
  hu <- fread(file.path(proj, "data", "healthuser_covariates.csv"))
  hu_vars <- setdiff(names(hu), "eid")
  df[hu, on = "eid", (hu_vars) := mget(paste0("i.", hu_vars))]  # update join keeps row order
  imp_vars <- c(imp_vars_gs, hu_vars)
}
tc <- time_cuts(cohort$df)  # same time grid as the primary analysis
r_probs <- if (variant == "knots4") c(0.05, 0.35, 0.65, 0.95) else c(0.1, 0.5, 0.9)
subs <- subs_grid()

# risks at end of follow-up for the reference and every substitution, averaged by cohort
sub_results <- function(dt, models, cohorts) {
  dt[, row_id := .I]
  final <- function(d) risk_by_time(d, models, tc$n_groups)[timegroup == tc$n_groups, risk]
  ref <- final(dt)
  rbindlist(lapply(seq_len(nrow(subs)), function(i) {
    s <- subs[i]
    sub <- substitute(dt, s$from, s$to, s$duration / 60, comp_limits)
    risk <- final(sub$data)
    rbindlist(lapply(names(cohorts), function(c) {
      idx <- cohorts[[c]]
      data.table(cohort = c, from = s$from, to = s$to, duration = s$duration,
                 ref_risk = mean(ref[idx]), sub_risk = mean(risk[idx]), prop_sub = mean(sub$ok[idx]))
    }))
  }))
}

run_once <- function(dt) {
  dt <- impute_pmm(add_ilr(copy(dt)), imp_vars)
  if (variant == "strat") {
    res <- rbindlist(lapply(c("short", "normal", "long"), function(g) {
      d <- dt[sleep_group == g]
      models <- fit_models_gs(d, tc$cuts, "reduced", r_probs)
      idx <- cohort_index(d)
      sub_results(d, models, idx[c(g, paste0(g, c("_female", "_male")))])
    }))
  } else {
    models <- fit_models_gs(dt, tc$cuts, if (variant %in% c("reduced", "nsaid", "healthuser")) variant else "primary", r_probs)
    idx <- cohort_index(dt)
    res <- sub_results(dt, models, idx[cohort_names])
  }
  res[, rr := sub_risk / ref_risk]
  res
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
    set.seed(50000 + b)
    res <- tryCatch(run_once(df[sample(.N, .N, replace = TRUE)]), error = function(e) NULL)
    if (!is.null(res)) saveRDS(res[, boot := b], file.path(out_dir, "boots", paste0(b, ".rds")))
    NULL
  }, mc.cores = n_cores, mc.preschedule = FALSE))
}

summary_dt <- copy(point)
boot_files <- list.files(file.path(out_dir, "boots"), full.names = TRUE)
if (length(boot_files) > 0) {
  boots <- rbindlist(lapply(boot_files, readRDS))
  ci <- boots[, .(lower = quantile(rr, 0.025), upper = quantile(rr, 0.975), n_boot = .N),
              by = .(cohort, from, to, duration)]
  summary_dt <- merge(summary_dt, ci, by = c("cohort", "from", "to", "duration"))
}
summary_dt[, reported := prop_sub >= intervention_threshold]
setorder(summary_dt, cohort, from, duration)
fwrite(summary_dt, file.path(out_dir, "risk_ratios.csv"))

cat(sprintf("%s: N = %d, gallstone = %d\n", variant, nrow(df), sum(df$gs)))
print(summary_dt[abs(duration) == 30 & cohort %in% c("short", "normal", "long", "full"),
                 c("cohort", "from", "duration", "rr", intersect(c("lower", "upper"), names(summary_dt)), "prop_sub"),
                 with = FALSE], digits = 3)
