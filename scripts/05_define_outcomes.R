# Outcome ascertainment with the FreeUKB R package (version 0.2.2, function ukb_outcome_custom), applied in the same
# way to every outcome:
#   - events from first-occurrence fields, hospital inpatient records (ICD-10) and death registrations (underlying and
#     contributory causes); the date of onset is the earliest record;
#   - participants with the condition at recruitment (a dated record before the recruitment visit, or a self-reported
#     illness code in field 20002) are excluded;
#   - follow-up ends at region-specific administrative censoring dates, taken as the earliest end date of the three
#     sources in each region (see fix_censor_date below).
#
#   name            ICD-10  self-report code (coding 6)
#   gallstone       K80     1162   primary outcome
#   gerd            K21     1138   secondary outcome
#   diverticular    K57     1458   secondary outcome
#   ibs             K58     1154   secondary outcome
#   gastric_ulcer   K25     1142   peptic ulcer = first of gastric or duodenal ulcer (06_build_analysis_datasets.py)
#   duodenal_ulcer  K26     1457
#   cyst            L72     -      negative control outcome (primary)
#   earwax          H61     -      negative control outcome (baseline hearing-aid users are excluded in 06)
#
# L72 and H61 have no self-report code in coding 6, but FreeUKB requires one; "-999" matches no participant, so only
# dated records before recruitment define baseline cases.
#
# Inputs (paths set in config.sh):
#   FREEUKB_INPUT        table exported with the UK Biobank RAP Table Exporter, with the columns listed by
#                        FreeUKB::ukb_required_fields("outcome_followup"): fields 53, 54, 190, 191, 20002, 40000,
#                        40001, 40002, 41270 and 41280
#   UKB_HEALTH_OUTCOMES  extract with the first-occurrence date fields listed below
# Output: data/outcomes/<name>/out.csv
# Usage:  Rscript scripts/05_define_outcomes.R [name ...]

suppressPackageStartupMessages({
  library(FreeUKB)
  library(data.table)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
outcomes <- list(
  gallstone      = list(code = "K80", field = "p131674", self_report = "1162", desc = "cholelithiasis"),
  gerd           = list(code = "K21", field = "p131584", self_report = "1138", desc = "gastro-oesophageal reflux disease"),
  diverticular   = list(code = "K57", field = "p131636", self_report = "1458", desc = "diverticular disease of intestine"),
  ibs            = list(code = "K58", field = "p131638", self_report = "1154", desc = "irritable bowel syndrome"),
  gastric_ulcer  = list(code = "K25", field = "p131590", self_report = "1142", desc = "gastric ulcer"),
  duodenal_ulcer = list(code = "K26", field = "p131592", self_report = "1457", desc = "duodenal ulcer"),
  cyst           = list(code = "L72", field = "p131794", self_report = "-999",
                        desc = "follicular cysts of skin and subcutaneous tissue"),
  earwax         = list(code = "H61", field = "p131224", self_report = "-999", desc = "other disorders of external ear")
)
args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else names(outcomes)
stopifnot(all(todo %in% names(outcomes)))

dt1 <- fread(Sys.getenv("FREEUKB_INPUT"))
fields <- unname(sapply(outcomes[todo], `[[`, "field"))
fo <- fread(Sys.getenv("UKB_HEALTH_OUTCOMES"), select = c("eid", fields), nThread = 2)
# first-occurrence columns named as FreeUKB expects them: Date.<ICD10>.first.reported..<description>.
setnames(fo, c("eid", fields), c("Participant.ID", unname(sapply(outcomes[todo], function(o)
  make.names(sprintf("Date %s first reported (%s)", o$code, o$desc))))))
codes <- unname(sapply(outcomes[todo], `[[`, "code"))
drop <- grep(paste0("^Date\\.(", paste(codes, collapse = "|"), ")\\.first\\.reported"), names(dt1), value = TRUE)
if (length(drop)) dt1[, (drop) := NULL]  # one first-occurrence column per outcome, taken from the extract
dat <- merge(dt1, fo, by = "Participant.ID", all.x = TRUE)

# Administrative censoring: the earliest end date of the requested sources in each region. Replaces FreeUKB's
# internal .freeukb_resolve_censor_date for this session only.
ns <- asNamespace("FreeUKB")
fix_censor_date <- function(region, source = NULL, censor_dates = NULL) {
  region <- toupper(trimws(as.character(region)))
  source <- tolower(source)
  censor_tbl <- if (!is.null(censor_dates)) censor_dates else .freeukb_censoring_dates()
  if (is.data.frame(censor_tbl)) {
    censor_tbl$region <- toupper(trimws(as.character(censor_tbl$region)))
    censor_tbl$source <- tolower(trimws(as.character(censor_tbl$source)))
    censor_tbl$end_date <- .freeukb_clean_dates(.freeukb_parse_date(censor_tbl$end_date))
    if (!is.null(source)) censor_tbl <- censor_tbl[censor_tbl$source %in% source, , drop = FALSE]
    by_region <- tapply(censor_tbl$end_date, censor_tbl$region, function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0L) return(as.Date(NA))
      min(x)
    })
    end_dates <- by_region[match(region, names(by_region))]
    return(as.Date(end_dates, origin = "1970-01-01"))
  }
  .freeukb_abort("censor_dates must be a data.frame or named vector.", class = "freeukb_error_input")
}
environment(fix_censor_date) <- ns
unlockBinding(".freeukb_resolve_censor_date", ns)
assign(".freeukb_resolve_censor_date", fix_censor_date, envir = ns)
lockBinding(".freeukb_resolve_censor_date", ns)

for (n in todo) {
  o <- outcomes[[n]]
  out <- ukb_outcome_custom(data = dat, outcome = setNames(list(o$code), n),
                            source = c("first_occurrence", "inpatient", "death"),
                            self_report_codes = setNames(list(o$self_report), n), keep = "all")
  dir.create(file.path(proj, "data", "outcomes", n), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file.path(proj, "data", "outcomes", n, "out.csv"), row.names = FALSE)
  ev <- grep("_event$", names(out), value = TRUE)
  message(n, ": ", nrow(out), " participants without the condition at recruitment, ",
          sum(out[[ev]] == 1, na.rm = TRUE), " events")
}
