# Extract the baseline and health-record fields used for covariates, exclusions and sensitivity analyses, restricted
# to participants with valid accelerometer data (data/accel_base.csv from 01_accelerometry.py).
#
# Inputs (config.sh): UKB_POPULATION, UKB_ASSESSMENT, UKB_HEALTH_OUTCOMES
# Outputs (all values kept as text; recoding happens in 03 and 06):
#   data/raw_population_fields.csv      age at recruitment, sex, Townsend deprivation index
#   data/raw_assessment_fields.csv      recruitment visit (instance 0): date, ethnicity, education, BMI, smoking,
#                                       alcohol, diet, income, employment, shift work, medications (fields 20003, 6153,
#                                       6177), reproductive factors and self-reported operations (field 20004)
#   data/raw_health_outcome_fields.csv  first-occurrence dates of the chronic diseases and of K50, K51 and K90,
#                                       OPCS-4 operations with dates (fields 41272, 41282) and the cancer registry
#                                       (fields 40005, 40006)
# Usage: Rscript scripts/02_extract_fields.R

suppressPackageStartupMessages(library(data.table))

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
ids <- fread(file.path(proj, "data", "accel_base.csv"), select = "eid")$eid

pick <- function(file, patterns) {
  hdr <- names(fread(file, nrows = 0))
  cols <- c("eid", grep(paste(patterns, collapse = "|"), hdr, value = TRUE))
  message(basename(file), ": reading ", length(cols), " columns")
  dt <- fread(file, select = cols, nThread = 4, colClasses = "character")
  dt[eid %in% as.character(ids)]
}

po <- pick(Sys.getenv("UKB_POPULATION"), c("^p21022$", "^p31$", "^p22189$"))
fwrite(po, file.path(proj, "data", "raw_population_fields.csv"))

ac <- pick(Sys.getenv("UKB_ASSESSMENT"), c(
  "^p53_i0$", "^p21000_i0$", "^p6138_i0", "^p21001_i0$",    # recruitment date, ethnicity, education, BMI
  "^p20116_i0$", "^p20117_i0$",                             # smoking, alcohol
  "^p1289_i0$", "^p1299_i0$", "^p1309_i0$", "^p1319_i0$",   # vegetables (cooked, raw), fruit (fresh, dried)
  "^p1329_i0$", "^p1339_i0$", "^p1349_i0$",                 # oily fish, non-oily fish, processed meat
  "^p1369_i0$", "^p1379_i0$", "^p1389_i0$",                 # beef, lamb, pork
  "^p1438_i0$", "^p1448_i0$", "^p1458_i0$", "^p1468_i0$",   # bread intake and type, cereal intake and type
  "^p738_i0$", "^p6142_i0$", "^p826_i0$",                   # income, employment, shift work
  "^p2734_i0$", "^p2784_i0$", "^p2814_i0$",                 # live births, oral contraceptive, HRT
  "^p6153_i0$", "^p6177_i0$",                               # medications (women / men)
  "^p20003_i0_", "^p20004_i0_"                              # medication codes, self-reported operations
))
fwrite(ac, file.path(proj, "data", "raw_assessment_fields.csv"))

chronic <- c(131296, 131298, 131300, 131302, 131304, 131306,          # I20-I25
             131360, 131362, 131364, 131366, 131368, 131378,          # I60-I64, I69
             130706, 130708, 130710, 130712, 130714,                  # E10-E14
             131286, 130814, 130904, 130906, 130894, 130896, 131060)  # I10, E78, F40-F41, F32-F33, G47
ho <- pick(Sys.getenv("UKB_HEALTH_OUTCOMES"), c(
  paste0("^p", chronic, "$"),
  "^p131626$", "^p131628$", "^p131688$",                    # K50, K51, K90
  "^p41272$", "^p41282_",                                   # OPCS-4 operations and dates
  "^p40005_i", "^p40006_i"                                  # cancer registry: date, ICD-10 type
))
fwrite(ho, file.path(proj, "data", "raw_health_outcome_fields.csv"))
message("done: ", nrow(po), " / ", nrow(ac), " / ", nrow(ho), " participants")
