# Pooled logistic g-formula components used by 13_sensitivity_models.R and 14_composition_search.R.
# 10_gformula.R and 15_natural_course.R define the same primary model inline (kept as run for the published results).

imp_vars_gs <- c("age_accel", "sex", "ethnic", "tdi", "edu", "bmi", "smoke", "alcohol", "diet_score",
                 "CVD", "diabetes", "hypertension", "Dyslipidemia", "Cancer", "anxiety", "depression",
                 "sleep_disorder", "income", "retired", "sick_disabled", "shift", "psych_meds",
                 "chol_med", "hrt", "ocp", "parity", "R1", "R2", "R3", "gs", "death", "time_to_gs")

time_cuts <- function(df) {
  max_fu <- max(df$time_to_gs)
  n_groups <- ceiling(max_fu / 365)
  list(cuts = seq(0, max_fu, length.out = n_groups + 1), n_groups = n_groups,
       years_per_group = max_fu / n_groups / 365.25)
}

# Harrell restricted cubic spline with any number of knots (k knots -> k - 1 columns);
# identical to rcs3() when three knots are given.
rcsk <- function(x, k) {
  n <- length(k)
  p3 <- function(u) pmax(u, 0)^3
  nl <- sapply(seq_len(n - 2), function(j) {
    (p3(x - k[j]) - p3(x - k[n - 1]) * (k[n] - k[j]) / (k[n] - k[n - 1]) +
       p3(x - k[n]) * (k[n - 1] - k[j]) / (k[n] - k[n - 1])) / (k[n] - k[1])^2
  })
  m <- cbind(x, nl)
  colnames(m) <- c("", paste0("nl", seq_len(n - 2)))
  m
}

# model: primary | s2 | s3 | nsaid (+ aspirin/NSAID use) | healthuser | reduced (main effects only; for small within-group models)
# r_probs: quantiles used as spline knots for the ILR coordinates
make_formula_gs <- function(long, model = "primary", r_probs = c(0.1, 0.5, 0.9)) {
  env <- new.env(parent = globalenv())
  env$rcs3 <- rcs3
  env$rcsk <- rcsk
  q <- function(v, p = c(0.1, 0.5, 0.9)) unname(quantile(long[[v]], p))
  env$kT <- q("timegroup", c(0.01, 0.5, 0.99))
  env$kA <- q("age_accel")
  env$kD <- q("tdi")
  env$kR1 <- q("R1", r_probs); env$kR2 <- q("R2", r_probs); env$kR3 <- q("R3", r_probs)
  if (model == "reduced") {
    rhs <- "rcsk(R1, kR1) + rcsk(R2, kR2) + rcsk(R3, kR3) + rcs3(timegroup, kT) + rcs3(age_accel, kA) + sex +
            smoke + income + retired + edu + ethnic + alcohol + diet_score + rcs3(tdi, kD) + shift + psych_meds"
    return(as.formula(paste("~", rhs), env = env))
  }
  rhs <- "
    rcsk(R1, kR1) * rcs3(timegroup, kT) + rcsk(R2, kR2) * rcs3(timegroup, kT) + rcsk(R3, kR3) * rcs3(timegroup, kT) +
    rcsk(R1, kR1) * sex + rcsk(R2, kR2) * sex + rcsk(R3, kR3) * sex +
    rcsk(R1, kR1) * smoke + rcsk(R2, kR2) * smoke + rcsk(R3, kR3) * smoke +
    rcsk(R1, kR1) * rcs3(age_accel, kA) + rcsk(R2, kR2) * rcs3(age_accel, kA) + rcsk(R3, kR3) * rcs3(age_accel, kA) +
    rcsk(R1, kR1) * income + rcsk(R2, kR2) * income + rcsk(R3, kR3) * income +
    rcsk(R1, kR1) * retired + rcsk(R2, kR2) * retired + rcsk(R3, kR3) * retired +
    rcs3(timegroup, kT) * rcs3(age_accel, kA) + rcs3(timegroup, kT) * sex +
    rcs3(age_accel, kA) * sex +
    edu + ethnic + alcohol + diet_score + rcs3(tdi, kD) + shift + psych_meds"
  if (model == "s2") {
    env$kB <- q("bmi")
    rhs <- paste(rhs, "+ rcs3(bmi, kB) + CVD + diabetes + hypertension + Dyslipidemia + Cancer +
                       anxiety + depression + sleep_disorder + sick_disabled")
  }
  if (model == "s3") rhs <- paste(rhs, "+ chol_med + hrt + ocp + parity")
  if (model == "nsaid") rhs <- paste(rhs, "+ aspirin + nsaid")
  if (model == "healthuser") rhs <- paste(rhs, "+ factor(srh) + longill + n_meds + bowel_screen + sex_screen")
  as.formula(paste("~", rhs), env = env)
}

fit_models_gs <- function(dt, cuts, model = "primary", r_probs = c(0.1, 0.5, 0.9)) {
  dt <- copy(dt)[, fu_time := time_to_gs]
  long <- survSplit(Surv(time_to_gs, gs) ~ ., data = as.data.frame(dt),
                    cut = cuts[-c(1, length(cuts))], episode = "timegroup",
                    start = "start", end = "end", event = "gs")
  setDT(long)
  long[, last := end == fu_time]
  long[, death := as.integer(last & death == 1)]
  f <- make_formula_gs(long, model, r_probs)
  strip <- function(m) {
    m[c("residuals", "fitted.values", "effects", "linear.predictors", "weights",
        "prior.weights", "data")] <- NULL
    m$qr$qr <- NULL
    m
  }
  list(gs = strip(glm(update(f, gs ~ .), data = long[death == 0], family = binomial(), model = FALSE, y = FALSE)),
       death = strip(glm(update(f, death ~ .), data = long[gs == 0], family = binomial(), model = FALSE, y = FALSE)))
}

# Cumulative incidence of gallstone (competing death) at every time group.
# Returns a long table with row_id, timegroup, risk; `dt` must carry `row_id`.
risk_by_time <- function(dt, models, n_groups) {
  n <- nrow(dt)
  long <- dt[rep(seq_len(n), each = n_groups)]
  long[, timegroup := rep(seq_len(n_groups), n)]
  long[, h_gs := predict(models$gs, newdata = .SD, type = "response")]
  long[, h_d := predict(models$death, newdata = .SD, type = "response")]
  long[, risk := cumsum(h_gs * cumprod((1 - shift(h_gs, fill = 0)) * (1 - h_d))), by = row_id]
  long[, .(row_id, timegroup, risk)]
}
