# Health-status and health-seeking indicators for the healthy-user sensitivity analysis:
# self-rated health, long-standing illness, number of medications, and cancer screening uptake
# (baseline assessment, instance 0), restricted to the accelerometer participants.
#
# Input (config.sh): UKB_ASSESSMENT
# Output: data/healthuser_covariates.csv (eid, srh, longill, n_meds, bowel_screen, sex_screen)
#         data/hearing_baseline.csv (eid, hearing_aid, hearing_difficulty), used to exclude
#         baseline hearing-aid users from the external ear disorders (H61) negative control outcome
# Usage:  Rscript scripts/04_extract_healthuser.R

suppressPackageStartupMessages(library(data.table))

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
ids <- fread(file.path(proj, "data", "accel_base.csv"), select = "eid")$eid

cols <- c("eid",
          "p2178_i0",   # overall health rating: 1 excellent, 2 good, 3 fair, 4 poor
          "p2188_i0",   # long-standing illness, disability or infirmity
          "p137_i0",    # number of treatments / medications taken
          "p2345_i0",   # ever had bowel cancer screening
          "p2674_i0",   # ever had breast cancer screening / mammogram (women)
          "p2365_i0",   # ever had PSA test (men)
          "p3393_i0",   # hearing aid user
          "p2247_i0")   # hearing difficulty / problems (99 = completely deaf)
dt <- fread(Sys.getenv("UKB_ASSESSMENT"), select = cols, nThread = 2)
dt <- dt[eid %in% ids]
fwrite(dt, file.path(proj, "data", "raw_healthuser_fields_all.csv"))

# negative codes (-1 do not know, -3 prefer not to answer) -> missing
neg_na <- function(x) replace(x, x < 0, NA)
out <- dt[, .(eid,
              srh = neg_na(p2178_i0),
              longill = neg_na(p2188_i0),
              n_meds = pmin(neg_na(p137_i0), 10),
              bowel_screen = neg_na(p2345_i0),
              # the mammogram question was asked of women and the PSA question of men
              sex_screen = fcoalesce(neg_na(p2674_i0), neg_na(p2365_i0)))]
fwrite(out, file.path(proj, "data", "healthuser_covariates.csv"))
fwrite(dt[, .(eid, hearing_aid = neg_na(p3393_i0),
              hearing_difficulty = neg_na(fifelse(p2247_i0 == 99, 1L, p2247_i0)))],
       file.path(proj, "data", "hearing_baseline.csv"))
message("done: ", nrow(out), " participants")
print(summary(out[, -"eid"]))
