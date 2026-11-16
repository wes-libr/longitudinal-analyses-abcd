# =============================================================================
# tests/run_mock.R : run the whole pipeline on simulated data (no ABCD data)
#   setwd("<...>/pipeline"); source("tests/run_mock.R")
# Outputs go to output_mock/ and a temporary derived folder.
# =============================================================================
options(lp.mock = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
source("tests/mock_nbdc.R")
# Use the real NBDCtools event filters: from the installed package if present,
# otherwise from a source checkout (NBDCTOOLS_SRC = path to NBDCtools/R).
if (requireNamespace("NBDCtools", quietly = TRUE)) {
  filter_events_abcd <- NBDCtools::filter_events_abcd
  filter_empty_rows  <- function(data) NBDCtools::filter_empty_rows(data, id_cols = c("participant_id", "session_id"))
} else {
  source(file.path(Sys.getenv("NBDCTOOLS_SRC"), "filter.R"))
  .fer <- filter_empty_rows
  filter_empty_rows <- function(data) .fer(data, id_cols = c("participant_id", "session_id"))
}
.mock$data <- simulate_abcd(n_fam = as.integer(Sys.getenv("LP_MOCK_NFAM", "1300")))
.lp_args <- list(dir_derived = Sys.getenv("LP_TEST_DERIVED", file.path(tempdir(), "lp_derived")), dir_output = "output_mock",
                 fast = TRUE, jm_iter = 2500L, jm_burnin = 500L, mi_m = 5L)
if (nzchar(Sys.getenv("LP_MODULES"))) .lp_args$modules <- strsplit(Sys.getenv("LP_MODULES"), ",")[[1]]
source("run_all.R")
