# Copy this file to config.sh and set the paths. config.sh is listed in .gitignore.
# UK Biobank data are not part of this repository; access is by application to UK Biobank.
#
# Participant-level extracts (one row per participant, column names as written by the UK Biobank Research Analysis
# Platform, e.g. eid, p31, p53_i0, p20003_i0_a0). The fields read from each file are listed in README.md.
export UKB_POPULATION="/path/to/population_characteristics.csv"   # fields 31, 34, 52, 191, 21022, 22189
export UKB_ASSESSMENT="/path/to/assessment_centre.csv"             # recruitment (instance 0) and imaging visit (instance 2)
export UKB_HEALTH_OUTCOMES="/path/to/health_outcomes.csv"          # first occurrences, operations, cancer registry
export UKB_ACCELEROMETER="/path/to/accelerometer.csv"              # fields 40034-40041 and 90002-90127

# Input table for FreeUKB::ukb_outcome_custom, exported with the RAP Table Exporter
# (columns listed by FreeUKB::ukb_required_fields("outcome_followup"); see scripts/05_define_outcomes.R)
export FREEUKB_INPUT="/path/to/freeukb_outcome_followup.csv"

# Parallel workers for the bootstrap
export N_CORES=8
