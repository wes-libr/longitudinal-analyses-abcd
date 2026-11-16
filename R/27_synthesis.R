# =============================================================================
# 27_synthesis.R : one table of headline estimates across modules
# =============================================================================
lp_log("27_synthesis")
fs <- list.files(cfg$dir_output, pattern = "^zz_key_.*\\.csv$", full.names = TRUE)
fs <- fs[!grepl("suppressed", fs)]
key <- bind_rows(lapply(fs, function(f) utils::read.csv(f, stringsAsFactors = FALSE) |> mutate(across(everything(), as.character))))
# add the LMM baseline slopes for the primary marker
lmm <- file.path(cfg$dir_output, "11_lmm_fixed.csv")
if (file.exists(lmm)) key <- bind_rows(
  utils::read.csv(lmm, stringsAsFactors = FALSE) |> filter(model == "M2_random_slope", term == "time") |>
    transmute(module = "11", model, outcome, term = "time (per year)", estimate, conf.low, conf.high, p.value, n_id) |>
    mutate(across(everything(), as.character)),
  key)
safe_write_csv(key |> relocate(module, model), "99_key_estimates")
