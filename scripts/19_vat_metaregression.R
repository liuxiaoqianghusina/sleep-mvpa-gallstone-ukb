# Meta-regression of the 30-min substitution log RRs (full cohort) on each disease's log HR per SD of
# MRI visceral fat (fat-depot Cox models), with a dispersion-adjusted p on k - 2 df.
# Needs the fat-depot Cox results of 17_fat_depots_cox.R for every disease.
# Output: results/heterogeneity/metareg_on_vat_gradient.csv
# Usage: Rscript scripts/19_vat_metaregression.R

suppressMessages(library(data.table))
proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
ds <- c("gallstone","peptic_ulcer","gerd","diverticular","ibs")
rd <- function(o) file.path(if (o=="gallstone") file.path(proj,"results") else file.path(proj,"results","outcomes",o))
x <- rbindlist(lapply(ds, function(o) {
  c <- fread(file.path(rd(o), if (o=="gallstone") "fat_depots_gallstone_cox.csv" else "fat_depots_outcome_cox.csv"))
  b <- rbindlist(lapply(list.files(file.path(rd(o),"primary_range_thr","boots"), full.names=TRUE), readRDS))[cohort=="full"]
  p <- readRDS(file.path(rd(o),"primary_range_thr","point_estimate.rds"))[cohort=="full"]
  rbindlist(lapply(c("VAT (per SD)","VAT, adjusted for BMI"), function(m)
    merge(p[, .(from, duration, rr)], b[, .(se = sd(log(rr))), by=.(from, duration)], by=c("from","duration"))[
      , `:=`(disease=o, vat_model=m, log_hr=log(c[model==m, hr]))]))
}))
out <- x[abs(duration)==30, {
  w <- 1/se^2; f <- lm(log(rr) ~ log_hr, weights=w); X <- cbind(1, log_hr)
  s <- sqrt(solve(t(X) %*% diag(w) %*% X)[2,2]); q <- sum(w*resid(f)^2); phi <- max(1, q/3)
  .(slope=coef(f)[2], p_fixed=2*pnorm(-abs(coef(f)[2]/s)), p_disp=2*pt(-abs(coef(f)[2]/(s*sqrt(phi))),3))
}, by=.(vat_model, from, duration)]
print(out[, lapply(.SD, function(v) if (is.numeric(v)) signif(v,2) else v)])
fwrite(out, file.path(proj,"results","heterogeneity","metareg_on_vat_gradient.csv"))
