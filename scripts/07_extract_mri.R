# Abdominal MRI phenotypes at the imaging visit (instance 2) for the gallstone cohort (data/analysis_gallstone.csv).
# Input (config.sh): UKB_ASSESSMENT
# Output: data/raw_mri_fields.csv
# Usage:  Rscript scripts/07_extract_mri.R

suppressPackageStartupMessages(library(data.table))

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
ids <- fread(file.path(proj, "data", "analysis_gallstone.csv"), select = "eid")$eid

fields <- c(
  "p53_i2",     # date of imaging visit
  "p54_i2",     # imaging centre
  "p50_i2",     # standing height
  "p21001_i2",  # BMI at imaging
  "p22407_i2",  # visceral adipose tissue volume (L)
  "p22408_i2",  # abdominal subcutaneous adipose tissue volume (L)
  "p22415_i2",  # total adipose tissue volume
  "p22434_i2",  # abdominal fat ratio
  "p22435_i2",  # muscle fat infiltration
  "p22436_i2",  # 10P liver PDFF (%)
  "p40061_i2",  # liver PDFF (%), alternative pipeline
  "p21090_i2",  # pancreas PDFF (%)
  "p23357_i2",  # VAT error indicator
  "p23358_i2"   # ASAT error indicator
)
file <- Sys.getenv("UKB_ASSESSMENT")
hdr <- names(fread(file, nrows = 0))
cols <- c("eid", intersect(fields, hdr))
message("missing fields: ", paste(setdiff(fields, hdr), collapse = ", "))
mri <- fread(file, select = cols, nThread = 4)
mri <- mri[eid %in% ids]
fwrite(mri, file.path(proj, "data", "raw_mri_fields.csv"))
message("done: ", nrow(mri))
