# =============================================================================
# 20_gamm_age.R : nonlinear age trajectories (GAMM via gamm4)
# =============================================================================
# Age (not wave) is the time axis. Smooth terms are estimated as lme4 random
# effects (gamm4), so the same nested structure as 11_lmm is used:
# site / family / participant (intercept + slope in time) + scanner.
# (mgcv::bam with participant "re" smooths was ~30 s per 1.6k participants in
# testing and scales poorly to 11k; gamm4 was ~12 s.)
#  G1  s(age) by sex                         normative curves + growth rates
#  G2  G1 + s(age_bl)                        cohort/selection check for the
#      accelerated longitudinal design: a non-flat s(age_bl) means the
#      between-person age gradient differs from within-person change.
#  G3  s(age) by eventual binge onset        do eventual bingers already differ
#      (and diverge)? descriptive companion to 24/25.
lp_log("20_gamm_age")
if (!need(c("gamm4", "mgcv"))) stop("gamm4 and mgcv required")
d <- img_data() |>
  mutate(sex_f = factor(ifelse(female == 1, "female", "male")),
         ev_f  = factor(ifelse(eventual == 1, "eventual_binge", "no_binge")))
re <- if (cfg$scanner_re) ~ (1 | site) + (1 | site:fam) + (1 + time | site:fam:id) + (1 | scanner) else
  ~ (1 | site) + (1 | site:fam) + (1 + time | site:fam:id)
curves <- list(); smooth_tests <- list()
for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]))
  if (n_distinct(dm$id) < 100) next
  base <- paste0(m, " ~ ", if (icv_term(m) != "") "icv_z + " else "")
  specs <- list(
    G1_sex      = paste0(base, "sex_f + s(age, by = sex_f, k = 8)"),
    G2_cohort   = paste0(base, "sex_f + s(age, by = sex_f, k = 8) + s(age_bl, k = 5)"),
    G3_eventual = paste0(base, "sex_f + ev_f + s(age, by = ev_f, k = 8)"))
  for (s in names(specs)) {
    fit <- try_fit(paste(m, s), gamm4::gamm4(as.formula(specs[[s]]), random = re, data = dm, REML = TRUE))
    if (is.null(fit)) next
    g <- fit$gam; st <- summary(g)$s.table
    smooth_tests[[paste(m, s)]] <- tibble(module = "20", model = s, marker = m, term = rownames(st),
                                          edf = st[, "edf"], F = st[, "F"], p.value = st[, "p-value"], n_id = n_distinct(dm$id))
    grp <- if (s == "G3_eventual") "ev_f" else "sex_f"
    rng <- dm |> group_by(g = .data[[grp]]) |>
      summarise(lo = quantile(age, .025), hi = quantile(age, .975), n_obs = n(), .groups = "drop")
    grid <- bind_rows(lapply(seq_len(nrow(rng)), function(i)
      tibble(g = rng$g[i], age = seq(ceiling(rng$lo[i] * 4) / 4, floor(rng$hi[i] * 4) / 4, by = .25), n_obs = rng$n_obs[i])))
    nd <- tibble(age = grid$age, icv_z = 0, age_bl = 10,
                 sex_f = if (grp == "sex_f") grid$g else factor("female", levels = levels(d$sex_f)),
                 ev_f  = if (grp == "ev_f") grid$g else factor("no_binge", levels = levels(d$ev_f)))
    p0 <- predict(g, nd, se.fit = TRUE); p1 <- predict(g, nd |> mutate(age = age + 0.01))
    curves[[paste(m, s)]] <- tibble(module = "20", model = s, marker = m, group = as.character(grid$g), age = nd$age,
                                    fit = as.numeric(p0$fit), se = as.numeric(p0$se.fit),
                                    rate_per_year = as.numeric(p1 - p0$fit) / 0.01, n_obs = grid$n_obs)
  }
}
safe_write_csv(bind_rows(smooth_tests), "20_gamm_smooth_tests", count_cols = "n_id")
safe_write_csv(bind_rows(curves), "20_gamm_curves", count_cols = "n_obs")
