# =============================================================================
# 12_glmm_gee.R : binge drinking over annual waves (tutorial GLMM/GEE methods)
# =============================================================================
# Subject-specific (GLMM, glmmTMB) vs population-averaged (GEE, family-
# clustered) models for any binge episode; negative-binomial hurdle for binge
# days; and the GEE re-fitted with the ABCD raked propensity weight so the
# prevalence curve refers to the target population rather than the sample.
lp_log("12_glmm_gee")
a <- load_derived("panel_annual") |>
  filter(!is.na(binge_any), !is.na(time)) |>
  mutate(fam = paste(site, fam, sep = "_"), age_c = age - 12, wave = factor(wave_yr)) |>
  arrange(fam, id, wave_yr)
fx <- list()
if (need("glmmTMB")) {
  g1 <- try_fit("GLMM binge", glmmTMB::glmmTMB(binge_any ~ age_c + I(age_c^2) + female + edu_cat +
                                                 (1 | site) + (1 | site:fam) + (1 | site:fam:id),
                                               family = binomial, data = a))
  if (!is.null(g1)) fx$g1 <- tidy_fixed(g1, "12", "glmm_binary", "binge_any")
  g2 <- try_fit("hurdle NB binge days", glmmTMB::glmmTMB(binge_days ~ age_c + female + (1 | site) + (1 | site:fam:id),
                                                         ziformula = ~ age_c + female, family = glmmTMB::truncated_nbinom2,
                                                         data = a |> filter(!is.na(binge_days))))
  if (!is.null(g2)) fx$g2 <- tidy_fixed(g2, "12", "hurdle_nb_days", "binge_days")
}
if (need("geepack")) {
  # A wave factor separates completely when a wave has (almost) no binge cases,
  # which blows the GEE estimates up; fit the GEEs from the first wave with at
  # least min_cell cases onwards (earlier waves are reported in the prevalence table)
  ok_w <- a |> group_by(wave_yr) |> summarise(nb = sum(binge_any), .groups = "drop") |> filter(nb >= cfg$min_cell)
  w0 <- if (nrow(ok_w)) min(ok_w$wave_yr) else Inf
  lp_log("  GEE fitted from wave ", w0, " (earlier waves have < ", cfg$min_cell, " binge cases)")
  ag <- a |> filter(wave_yr >= w0) |> mutate(wave = droplevels(wave))
  idf <- factor(ag$fam, levels = unique(ag$fam))
  e1 <- try_fit("GEE binge", geepack::geeglm(binge_any ~ wave + female, id = idf, data = ag,
                                             family = binomial, corstr = "exchangeable"))
  if (!is.null(e1)) fx$e1 <- tidy_fixed(e1, "12", "gee_family_clustered", "binge_any")
  aw <- ag |> filter(!is.na(wt), wt > 0) |> mutate(wave = droplevels(wave))
  e2 <- try_fit("weighted GEE binge", geepack::geeglm(binge_any ~ wave + female, id = factor(aw$fam, levels = unique(aw$fam)),
                                                      data = aw, weights = wt / mean(wt), family = binomial, corstr = "independence"))
  if (!is.null(e2)) fx$e2 <- tidy_fixed(e2, "12", "gee_weighted_population", "binge_any")
  # Prevalence by wave, sample vs weighted
  prev <- a |> group_by(wave) |>
    summarise(n_assessed = n(), n_binge = sum(binge_any),
              prev_sample = mean(binge_any),
              prev_weighted = if (any(!is.na(wt))) weighted.mean(binge_any[!is.na(wt)], wt[!is.na(wt)]) else NA_real_,  # drop rows without a weight
              .groups = "drop")
  safe_write_csv(prev, "12_binge_prevalence_by_wave", count_cols = c("n_assessed", "n_binge"))
}
safe_write_csv(bind_rows(fx) |> mutate(n_id = n_distinct(a$id)), "12_binge_models", count_cols = "n_id")
