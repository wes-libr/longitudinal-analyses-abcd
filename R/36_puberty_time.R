# =============================================================================
# 36_puberty_time.R : pubertal development as the time axis
# =============================================================================
# Pubertal maturation may drive both brain change and the onset of drinking.
# Uses the continuous youth PDS mean score (cfg$pds_mean, sex-specific items).
#  P1 brain: which clock fits better? gamm4 models with s(age), s(PDS), and
#     both, compared by AIC; smooth tests for PDS given age.
#  P2 pubertal TIMING: each youth's mean residual from a sex-specific smooth
#     of PDS on age (positive = earlier maturer), standardized. LMM:
#     marker ~ time * timing; does earlier maturation mean faster thinning?
#  P3 onset: Cox models with (a) pubertal timing and (b) PDS as a time-varying
#     covariate (counting process, value at each assessment); family-clustered
#     SEs. Fitted separately because status = age + timing makes them collinear.
# Truth in tests/mock_nbdc.R: log HR +0.3 per SD of earlier timing.
lp_log("36_puberty_time")
if (!need(c("lme4", "survival", "mgcv"))) stop("lme4, survival, mgcv required")
core <- load_derived("panel_core") |> filter(!is.na(pds_mean), !is.na(age))
if (!nrow(core)) stop("No continuous PDS scores (cfg$pds_mean)")
tim <- bind_rows(lapply(c(0, 1), function(s) {
  dd <- core |> filter(female == s)
  g <- mgcv::gam(pds_mean ~ s(age, k = 6), data = dd)
  dd |> mutate(res = pds_mean - predict(g, dd))
})) |> group_by(id) |> summarise(timing = mean(res), n_pds = n(), .groups = "drop") |>
  mutate(timing_z = as.numeric(scale(timing)))
save_derived(tim, "puberty_timing")
d <- img_data() |> left_join(tim |> select(id, timing_z), by = "id")
out <- list()
## P1
if (need("gamm4")) for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]), !is.na(pds_mean))
  base <- paste0(m, " ~ female", icv_term(m))
  re <- ~ (1 | site) + (1 | site:fam) + (1 + time | site:fam:id)
  fits <- list(age = paste(base, "+ s(age, k = 6)"), puberty = paste(base, "+ s(pds_mean, k = 6)"),
               both = paste(base, "+ s(age, k = 6) + s(pds_mean, k = 6)"))
  for (nm in names(fits)) {
    g <- try_fit(paste("P1", m, nm), gamm4::gamm4(as.formula(fits[[nm]]), random = re, data = dm, REML = FALSE))
    if (is.null(g)) next
    st <- summary(g$gam)$s.table
    out[[paste("P1", m, nm)]] <- tibble(module = "36", model = paste0("P1_clock_", nm), outcome = m, term = rownames(st),
                                        estimate = st[, "edf"], p.value = st[, "p-value"], AIC = AIC(g$mer), n_id = n_distinct(dm$id))
  }
}
## P2
for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]), !is.na(timing_z))
  fit <- try_fit(paste("P2", m), lme4::lmer(as.formula(paste0(m, " ~ time * timing_z + age_bl_c + female", icv_term(m), " + ", re_terms("time"))),
                                            data = dm, control = lme4::lmerControl(calc.derivs = FALSE)))
  if (is.null(fit)) next
  if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
  out[[paste("P2", m)]] <- tidy_fixed(fit, "36", "P2_timing_lmm", m) |> filter(grepl("timing", term)) |> mutate(n_id = n_distinct(dm$id))
}
## P3
ev <- load_derived("events") |> filter(!is.na(edu_bl)) |> left_join(tim |> select(id, timing_z), by = "id") |>
  filter(!is.na(timing_z)) |> mutate(age_bl_c = age_bl - 10, fam = paste(site, fam, sep = "_"))
tm <- survival::tmerge(ev |> select(id, fam, female, age_bl_c, edu_bl, timing_z, t_event, event), ev, id = id, onset = event(t_event, event))
tm <- survival::tmerge(tm, core |> filter(id %in% ev$id) |> select(id, time, pds_mean), id = id, pds_tv = tdc(time, pds_mean))
tm <- tm |> filter(!is.na(pds_tv))
# P3a timing only (baseline covariate); P3b current pubertal status only. Kept
# separate: status at time t is age + timing, so both together are collinear.
p3 <- list(P3a_cox_puberty_timing = survival::Surv(tstart, tstop, onset) ~ timing_z + female + age_bl_c + edu_bl,
           P3b_cox_timevarying_pds = survival::Surv(tstart, tstop, onset) ~ pds_tv + female + age_bl_c + edu_bl)
for (nm in names(p3)) {
  cx <- try_fit(nm, survival::coxph(p3[[nm]], data = tm, cluster = fam))
  if (!is.null(cx)) out[[nm]] <- broom::tidy(cx, conf.int = TRUE) |>
    transmute(module = "36", model = nm, outcome = "binge_onset", term, estimate, std.error = robust.se,
              conf.low, conf.high, p.value, n_id = nrow(ev), n_events = sum(ev$event))
}
res <- bind_rows(out)
safe_write_csv(res, "36_puberty_time", count_cols = c("n_id", "n_events"))
if (nrow(res)) key_write(res |> filter(grepl("^P3", model) & term %in% c("timing_z", "pds_tv") | (model == "P2_timing_lmm" & outcome == cfg$primary)) |>
                           select(model, outcome, term, estimate, conf.low, conf.high, p.value, n_id), "36")
