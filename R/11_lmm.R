# =============================================================================
# 11_lmm.R : linear mixed models for brain trajectories (tutorial methods,
# corrected: nested site/family/participant + crossed scanner, time in years)
# =============================================================================
# M1 random intercept; M2 random slope; M3 time-invariant covariates (sex,
# caregiver education); M4 time-varying puberty split into within/between
# (person-mean centering); M5 cohort vs ageing (baseline age x time) for the
# accelerated longitudinal design: if cross-sectional and longitudinal age
# effects converge, the age_bl_c and time coefficients agree.
lp_log("11_lmm")
if (!need("lme4")) stop("lme4 required")
d <- img_data() |>
  group_by(id) |> mutate(pds_pm = mean(pds, na.rm = TRUE), pds_cwc = pds - pds_pm) |> ungroup()
fx <- list(); vc <- list()
for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]))
  if (n_distinct(dm$id) < 50) next
  cov <- paste0(" + age_bl_c + female", icv_term(m))
  specs <- list(
    M1_random_intercept = paste0(m, " ~ time", cov, " + ", re_terms()),
    M2_random_slope     = paste0(m, " ~ time", cov, " + ", re_terms("time")),
    M3_time_invariant   = paste0(m, " ~ time * (female + edu_cat) + age_bl_c", icv_term(m), " + ", re_terms("time")),
    M4_time_varying_pds = paste0(m, " ~ time + pds_cwc + pds_pm", cov, " + ", re_terms("time")),
    M5_cohort_vs_ageing = paste0(m, " ~ time * age_bl_c + female", icv_term(m), " + ", re_terms("time")))
  for (s in names(specs)) {
    dat <- if (s == "M4_time_varying_pds") dm |> filter(!is.na(pds_cwc)) else if (s == "M3_time_invariant") dm |> filter(!is.na(edu_cat)) else dm
    fit <- try_fit(paste(m, s), lme4::lmer(as.formula(specs[[s]]), data = dat, REML = TRUE,
                                           control = lme4::lmerControl(calc.derivs = FALSE)))
    if (is.null(fit)) next
    if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
    fx[[paste(m, s)]] <- tidy_fixed(fit, "11", s, m) |> mutate(n_id = n_distinct(dat$id), n_obs = nrow(dat))
    vc[[paste(m, s)]] <- tidy_vc(fit, "11", s, m)
    if (s == "M2_random_slope") save_derived(fit, paste0("lmm_", m))
  }
}
safe_write_csv(bind_rows(fx), "11_lmm_fixed", count_cols = c("n_id", "n_obs"))
safe_write_csv(bind_rows(vc), "11_lmm_variance_components")
