# =============================================================================
# 02_describe_missing.R : retention, missingness patterns, predictors of dropout
# =============================================================================
lp_log("02_describe_missing")
img <- load_derived("panel_img"); annual <- load_derived("panel_annual"); baseline <- load_derived("baseline")
waves <- c("ses-00A", "ses-02A", "ses-04A", "ses-06A")

# Imaging pattern per participant (1 = QC-passed scan at that wave)
pat <- img |> distinct(id, session_id) |> mutate(v = 1L) |>
  pivot_wider(names_from = session_id, values_from = v, values_fill = 0L)
for (w in setdiff(waves, names(pat))) pat[[w]] <- 0L
pat <- pat |> mutate(pattern = paste0(.data[["ses-00A"]], .data[["ses-02A"]], .data[["ses-04A"]], .data[["ses-06A"]]),
                     n_scans = .data[["ses-00A"]] + .data[["ses-02A"]] + .data[["ses-04A"]] + .data[["ses-06A"]],
                     last_obs = apply(cbind(.data[["ses-00A"]], .data[["ses-02A"]], .data[["ses-04A"]], .data[["ses-06A"]]), 1, function(r) max(which(r == 1)) * 2 - 2),
                     monotone = vapply(strsplit(pattern, ""), function(s) { s <- as.integer(s); all(diff(s) <= 0) }, logical(1)))
save_derived(pat, "missing_patterns")
safe_write_csv(pat |> count(pattern, name = "n") |> arrange(desc(n)), "02_imaging_patterns", count_cols = "n")
safe_write_csv(pat |> count(monotone, name = "n"), "02_monotone_vs_intermittent", count_cols = "n")

# Predictors of missing the next imaging wave (discrete-time attrition model),
# family-clustered GEE. Uses the participant's current marker and externalizing.
ann_ext <- annual |> select(id, session_id, ext)
nxt <- img |> select(id, session_id, wave_yr, all_of(cfg$primary), icv_z) |>
  left_join(ann_ext, by = c("id", "session_id")) |>
  left_join(baseline |> select(id, fam, female, age_bl, edu_bl), by = "id") |>
  filter(wave_yr < 6) |>
  mutate(next_obs = as.integer(paste(id, sprintf("ses-%02dA", wave_yr + 2)) %in% paste(img$id, img$session_id)),
         marker_z = as.numeric(scale(.data[[cfg$primary]])), ext_z = as.numeric(scale(ext))) |>
  filter(!is.na(marker_z), !is.na(ext_z), !is.na(edu_bl)) |>
  arrange(fam, id, wave_yr)
if (need("geepack")) {
  fit <- try_fit("attrition GEE", geepack::geeglm(next_obs ~ factor(wave_yr) + marker_z + ext_z + female + age_bl + edu_bl,
                                                 id = factor(fam, levels = unique(fam)), data = nxt,
                                                 family = binomial, corstr = "exchangeable"))
  if (!is.null(fit)) safe_write_csv(tidy_fixed(fit, "02", "attrition_gee", "observed_next_wave") |>
                                      mutate(n_obs = nrow(nxt)), "02_attrition_predictors", count_cols = "n_obs")
}

# Binge assessment completeness by annual wave
safe_write_csv(annual |> group_by(session_id) |>
                 summarise(n_rows = n(), n_binge_assessed = sum(!is.na(binge_any)), .groups = "drop"),
               "02_binge_assessment_by_wave", count_cols = c("n_rows", "n_binge_assessed"))
