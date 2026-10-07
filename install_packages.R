# R packages used by the pipeline (versions used for the paper in brackets; R 4.2.2).
#   data.table (1.15.4), survival (3.4.0), mice (3.16.0), ks (1.14.2); parallel is part of base R.
install.packages(c("data.table", "survival", "mice", "ks"))

# FreeUKB (0.2.2, MIT licence) is used only by scripts/05_define_outcomes.R. It is not on CRAN; install it from the
# source archive distributed by its developers, for example:
#   install.packages("FreeUKB_0.2.2.tar.gz", repos = NULL, type = "source")
