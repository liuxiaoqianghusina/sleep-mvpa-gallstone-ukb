# "Lowest risk", "typical" and "highest risk" 24-h compositions for gallstone disease,
# following Yiallourou et al. 2025 (cum_targets.R / ideal_comp.R):
#   1. grid of synthetic compositions in 15-min steps within the 2.5th-97.5th percentile
#      of each behaviour, summing to 24 h; drop implausible ones (KDE density below the
#      5th percentile of observed compositions);
#   2. training half: fit the g-formula models and estimate risk under every synthetic
#      composition, standardised to the covariate distribution; pick lowest / highest risk
#      and the most common (highest density) composition;
#   3. testing half: bootstrap cumulative incidence under each of the three compositions.
#
# Usage: OUTCOME=<outcome> Rscript scripts/14_composition_search.R [n_boot=200] [cores=3]

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(mice)
  library(parallel)
  library(ks)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(proj, "scripts", "helpers.R"))
source(file.path(proj, "scripts", "gformula_core.R"))

args <- commandArgs(trailingOnly = TRUE)
n_boot <- if (length(args) >= 1) as.integer(args[1]) else 200L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
out_dir <- file.path(results_root(proj), "ideal_composition")
dir.create(file.path(out_dir, "boots"), recursive = TRUE, showWarnings = FALSE)

step_min <- 15
uni_threshold <- 0.025
synth_threshold <- 0.05

df <- add_ilr(load_cohort(proj, exposure = "thr", limits = "range")$df)
tc <- time_cuts(df)

# ---- 1. synthetic compositions ----
mins <- df[, lapply(.SD, function(x) x * 60), .SDcols = behaviours]
rng <- lapply(behaviours, function(b) {
  q <- quantile(mins[[b]], c(uni_threshold, 1 - uni_threshold))
  seq(ceiling(q[1] / step_min) * step_min, floor(q[2] / step_min) * step_min, by = step_min)
})
names(rng) <- behaviours
synth <- do.call(CJ, rng)[sleep + inactivity + light + mvpa == 1440]
synth[, (behaviours) := lapply(.SD, function(x) x / 60), .SDcols = behaviours]
add_ilr(synth)

obs <- as.matrix(df[, .(R1, R2, R3)])
centre <- colMeans(obs)
scl <- apply(obs, 2, sd)
fit_kde <- kde(scale(obs, centre, scl), gridsize = rep(100, 3))
dens_threshold <- quantile(log(predict(fit_kde, x = scale(obs, centre, scl))), synth_threshold)
synth[, dens := log(predict(fit_kde, x = scale(as.matrix(synth[, .(R1, R2, R3)]), centre, scl)))]
n_grid <- nrow(synth)
synth <- synth[is.finite(dens) & dens > dens_threshold]
message(sprintf("synthetic compositions: %d on the grid, %d plausible", n_grid, nrow(synth)))

# ---- 2. training half ----
set.seed(2026)
train_idx <- sample(nrow(df), floor(nrow(df) / 2))
train <- impute_pmm(df[train_idx], imp_vars_gs)
test <- df[-train_idx]
train_models <- fit_models_gs(train, tc$cuts)

mode_of <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
train[, age_round := round(age_accel)]
covar <- train[, .(prop = .N / nrow(train)), by = .(income, sex, retired, smoke, age_accel = age_round)]
covar[, `:=`(edu = mode_of(train$edu), ethnic = mode_of(train$ethnic), alcohol = mode_of(train$alcohol),
             shift = mode_of(train$shift), psych_meds = mode_of(train$psych_meds),
             diet_score = mean(train$diet_score), tdi = mean(train$tdi))]
covar[, row_id := .I]

synth_risk_file <- file.path(out_dir, "synthetic_composition_risk.csv")
if (!file.exists(synth_risk_file)) {
  chunks <- split(seq_len(nrow(synth)), ceiling(seq_len(nrow(synth)) / 40))
  synth[, risk := unlist(mclapply(chunks, function(rows) {
    vapply(rows, function(i) {
      pd <- copy(covar)[, `:=`(R1 = synth$R1[i], R2 = synth$R2[i], R3 = synth$R3[i])]
      r <- risk_by_time(pd, train_models, tc$n_groups)[timegroup == tc$n_groups]
      weighted.mean(r$risk, pd$prop[r$row_id])
    }, numeric(1))
  }, mc.cores = n_cores))]
  fwrite(synth, synth_risk_file)
}
synth <- fread(synth_risk_file)

ref <- rbind(
  synth[order(risk)][1][, composition := "Lowest risk"],
  synth[order(-dens)][1][, composition := "Typical"],
  synth[order(-risk)][1][, composition := "Highest risk"]
)
fwrite(ref, file.path(out_dir, "reference_compositions.csv"))

# ---- 3. testing half ----
cum_risk <- function(dt) {
  dt <- impute_pmm(dt, imp_vars_gs)
  models <- fit_models_gs(dt, tc$cuts)
  dt[, row_id := .I]
  rbindlist(lapply(seq_len(nrow(ref)), function(k) {
    d <- copy(dt)[, `:=`(R1 = ref$R1[k], R2 = ref$R2[k], R3 = ref$R3[k])]
    r <- risk_by_time(d, models, tc$n_groups)
    r[, .(risk = mean(risk)), by = timegroup][, composition := ref$composition[k]]
  }))
}

point_file <- file.path(out_dir, "test_point.rds")
if (!file.exists(point_file)) {
  set.seed(2027)
  saveRDS(cum_risk(test), point_file)
}
point <- readRDS(point_file)

if (n_boot > 0) {
  todo <- setdiff(seq_len(n_boot), as.integer(sub("\\.rds$", "", list.files(file.path(out_dir, "boots")))))
  invisible(mclapply(todo, function(b) {
    set.seed(30000 + b)
    res <- tryCatch(cum_risk(test[sample(.N, .N, replace = TRUE)]), error = function(e) NULL)
    if (!is.null(res)) saveRDS(res[, boot := b], file.path(out_dir, "boots", paste0(b, ".rds")))
    NULL
  }, mc.cores = n_cores, mc.preschedule = FALSE))
}

out <- copy(point)
boot_files <- list.files(file.path(out_dir, "boots"), full.names = TRUE)
if (length(boot_files) > 0) {
  boots <- rbindlist(lapply(boot_files, readRDS))
  out <- merge(out, boots[, .(lower = quantile(risk, 0.025), upper = quantile(risk, 0.975)),
                          by = .(composition, timegroup)], by = c("composition", "timegroup"))
  final <- dcast(boots[timegroup == tc$n_groups], boot ~ composition, value.var = "risk")
  contrasts <- data.table(
    contrast = c("Lowest vs typical", "Highest vs typical"),
    rr = c(point[timegroup == tc$n_groups & composition == "Lowest risk", risk] /
             point[timegroup == tc$n_groups & composition == "Typical", risk],
           point[timegroup == tc$n_groups & composition == "Highest risk", risk] /
             point[timegroup == tc$n_groups & composition == "Typical", risk]),
    lower = c(quantile(final$`Lowest risk` / final$Typical, 0.025),
              quantile(final$`Highest risk` / final$Typical, 0.025)),
    upper = c(quantile(final$`Lowest risk` / final$Typical, 0.975),
              quantile(final$`Highest risk` / final$Typical, 0.975)),
    n_boot = nrow(final)
  )
  fwrite(contrasts, file.path(out_dir, "contrasts_vs_typical.csv"))
  print(contrasts, digits = 3)
}
out[, years := timegroup * tc$years_per_group]
fwrite(out, file.path(out_dir, "cumulative_incidence.csv"))

cat("\nReference compositions (hours/day):\n")
print(ref[, .(composition, sleep, inactivity, light, mvpa, risk_train = round(risk, 4), dens = round(dens, 2))])
cat("\nCumulative incidence at end of follow-up (testing half):\n")
print(out[timegroup == tc$n_groups], digits = 3)
