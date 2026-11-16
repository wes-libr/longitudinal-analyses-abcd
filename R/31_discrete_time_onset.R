# =============================================================================
# 31_discrete_time_onset.R : discrete-time hazard model for binge onset
# =============================================================================
# Onset is known only to fall between consecutive assessments. One row per
# at-risk interval (t_{k-1}, t_k] between annual/mid-year assessments; the
# event is a positive assessment at t_k. A complementary log-log model with
# offset log(interval length) is a grouped-time proportional-hazards model
# (piecewise-constant hazard), so exp(coef) is a hazard ratio on the
# continuous-time scale and each interval's exposure can be aligned with what
# actually happened during it:
#   * marker value: last scan before the interval starts (z), the predictive
#     counterpart of the joint model without measurement-error correction
#   * pandemic: the FRACTION of the interval in the restricted era, from
#     calendar dates (fixes the misdating that defeated 28_temporal_trends T4)
#   * age at interval start (1-year bands = baseline hazard), sex, education,
#     birth cohort
# Family-clustered GEE (independence) for inference.
lp_log("31_discrete_time_onset")
if (!need("geepack")) stop("geepack required")
dec_year <- function(x) as.numeric(format(x, "%Y")) + (as.numeric(format(x, "%j")) - 1) / 365.25
c0 <- dec_year(as.Date(cfg$covid_start)); c1 <- dec_year(as.Date(cfg$covid_end))
core <- load_derived("panel_core") |> filter(!is.na(binge_any), !is.na(time)) |> arrange(id, time)
img <- load_derived("panel_img") |> select(id, time, mz = all_of(paste0(cfg$primary, "_z"))) |> filter(!is.na(mz))
pp <- core |> group_by(id) |>
  mutate(t_start = lag(time), cal_start = lag(cal), cal_end = cal, prev_pos = cumsum(lag(binge_any, default = 0))) |>
  ungroup() |> filter(!is.na(t_start), prev_pos == 0, time > t_start) |>
  mutate(len = time - t_start, event = binge_any,
         frac_restricted = ifelse(is.na(cal_start) | is.na(cal_end), NA,
                                  pmax(0, pmin(cal_end, c1) - pmax(cal_start, c0)) / pmax(cal_end - cal_start, 1e-6)),
         age_start = age - len, age_band = factor(pmin(17, pmax(12, floor(age_start)))))  # bands below 12 have almost no onsets: pooled into "12" (<= 12) so the baseline hazard is estimable
# last scan at or before the interval start
pp <- pp |> left_join(img, by = "id", relationship = "many-to-many") |>
  filter(is.na(time.y) | time.y <= t_start + 1e-8) |>
  group_by(id, t_start) |> slice_max(time.y, n = 1, with_ties = FALSE) |> ungroup() |>
  rename(time = time.x) |>
  mutate(fam = paste(site, fam, sep = "_"), cohort_c = birth_year - mean(birth_year, na.rm = TRUE)) |>
  filter(!is.na(mz), !is.na(edu_bl), !is.na(frac_restricted)) |> arrange(fam, id, t_start)
save_derived(pp, "person_period")
fit <- try_fit("cloglog discrete-time", geepack::geeglm(event ~ age_band + mz + frac_restricted + cohort_c + female + edu_bl + offset(log(len)),
                                                        id = factor(fam, levels = unique(fam)), data = pp,
                                                        family = binomial(link = "cloglog"), corstr = "independence"))
if (!is.null(fit)) {
  tb <- tidy_fixed(fit, "31", "discrete_time_cloglog", "binge_onset") |>
    mutate(HR = exp(estimate), HR_low = exp(conf.low), HR_high = exp(conf.high),
           n_intervals = nrow(pp), n_id = n_distinct(pp$id), n_events = sum(pp$event))
  safe_write_csv(tb, "31_discrete_time_onset", count_cols = c("n_id", "n_events"))
  key_write(tb |> filter(term %in% c("mz", "frac_restricted", "cohort_c")) |>
              transmute(model, term = ifelse(term == "frac_restricted", "restricted era (log HR, full interval)", term),
                        estimate, conf.low, conf.high, p.value, n_id, n_events), "31")
}
safe_write_csv(pp |> group_by(age_band) |> summarise(n_intervals = n(), n_events = sum(event),
                                                     mean_length_yr = mean(len), .groups = "drop"),
               "31_intervals_by_age", count_cols = c("n_intervals", "n_events"))
