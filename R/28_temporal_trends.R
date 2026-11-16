# =============================================================================
# 28_temporal_trends.R : age, period (calendar time, pandemic) and cohort
# =============================================================================
# ABCD visits span ~2016-2025, so "trajectory with age" is entangled with
# calendar time: the pandemic (remote visits, closures, reported accelerated
# brain maturation), secular declines in adolescent drinking, and scanner /
# software drift. Birth years span ~3 years, so cohort is also estimable.
# Age = period - cohort, so linear APC terms are not jointly identified; we use
# (i) nonlinear period smooths with age smooths, (ii) pandemic exposure terms
# that are non-linear in calendar time, and (iii) explicit two-way models
# (age+cohort, age+period) compared by AIC.
#
#  T1 brain: gamm4  y ~ s(age) + s(cal)            period smooth (drift/pandemic)
#  T2 brain: LMM    y ~ time + age_bl + pand_yrs    years of pandemic exposure at
#                   scan (accelerated-maturation test), + cohort (birth year)
#  T3 binge: GLMM   binge ~ age + covid_era + remote + birth cohort
#  T4 onset: Cox, counting process split at the pandemic boundaries
#            (time-varying era), + birth cohort
#  T5 joint: JMbayes2 joint model (as J1) with the time-varying pandemic era
#            in the survival submodel (counting-process Surv)
lp_log("28_temporal_trends")
if (!need(c("lme4", "glmmTMB", "survival"))) stop("lme4, glmmTMB, survival required")
d <- img_data() |> mutate(cohort_c = birth_year - mean(birth_year, na.rm = TRUE), cal_c = cal - 2019)
if (all(is.na(d$cal))) stop("No visit dates (", cfg$visit_dt, "): temporal-trend module skipped")
dec_year <- function(x) as.numeric(format(x, "%Y")) + (as.numeric(format(x, "%j")) - 1) / 365.25
c0 <- dec_year(as.Date(cfg$covid_start)); c1 <- dec_year(as.Date(cfg$covid_end))
out <- list(); key <- list()

## T1 period smooth for each marker
if (need("gamm4")) for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]), !is.na(cal))
  f <- as.formula(paste0(m, " ~ female", icv_term(m), " + s(age, k = 8) + s(cal, k = 6)"))
  g <- try_fit(paste("T1 period smooth", m), gamm4::gamm4(f, random = ~ (1 | site) + (1 | site:fam) + (1 + time | site:fam:id) + (1 | scanner), data = dm))
  if (is.null(g)) next
  st <- summary(g$gam)$s.table
  out[[paste("T1s", m)]] <- tibble(module = "28", model = "T1_period_smooth", outcome = m, term = rownames(st),
                                   estimate = st[, "edf"], statistic = st[, "F"], p.value = st[, "p-value"], n_id = n_distinct(dm$id))
  cg <- seq(ceiling(quantile(dm$cal, .02) * 4) / 4, floor(quantile(dm$cal, .98) * 4) / 4, by = .25)
  nd <- tibble(cal = cg, age = median(dm$age), female = 0, icv_z = 0)
  p <- predict(g$gam, nd, se.fit = TRUE)
  out[[paste("T1c", m)]] <- tibble(module = "28", model = "T1_period_curve", outcome = m, term = "calendar_year",
                                   cal = cg, estimate = as.numeric(p$fit) - mean(p$fit), std.error = as.numeric(p$se.fit),
                                   n_id = n_distinct(dm$id))
}

## T2 pandemic exposure and cohort in brain LMMs (age-cohort / age-period / APC-step)
for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]), !is.na(cal))
  cov <- paste0(" + female", icv_term(m), " + ", re_terms("time"))
  specs <- list(T2a_age_cohort = paste0(m, " ~ time + age_bl_c + cohort_c", cov),
                T2b_age_pandemic = paste0(m, " ~ time + age_bl_c + pand_yrs", cov),
                T2c_age_cohort_pandemic = paste0(m, " ~ time + age_bl_c + cohort_c + pand_yrs", cov))
  for (s in names(specs)) {
    fit <- try_fit(paste(m, s), lme4::lmer(as.formula(specs[[s]]), data = dm, REML = FALSE, control = lme4::lmerControl(calc.derivs = FALSE)))
    if (is.null(fit)) next
    if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
    out[[paste(m, s)]] <- tidy_fixed(fit, "28", s, m) |> filter(term %in% c("time", "age_bl_c", "cohort_c", "pand_yrs")) |>
      mutate(AIC = AIC(fit), n_id = n_distinct(dm$id))
  }
}

## T3 binge over calendar time
a <- load_derived("panel_annual") |>
  filter(!is.na(binge_any), !is.na(cal)) |>
  mutate(fam = paste(site, fam, sep = "_"), age_c = age - 12, cohort_c = birth_year - mean(birth_year, na.rm = TRUE))
if (nrow(a)) {
  t3 <- try_fit("T3 binge APC-step", glmmTMB::glmmTMB(binge_any ~ age_c + I(age_c^2) + covid_era + remote + cohort_c + female +
                                                        (1 | site) + (1 | site:fam:id), family = binomial, data = a))
  if (!is.null(t3)) out$T3 <- tidy_fixed(t3, "28", "T3_binge_glmm_period", "binge_any") |> mutate(n_id = n_distinct(a$id))
  bins <- a |> mutate(half_year = floor(cal * 2) / 2) |> group_by(half_year) |>
    summarise(n_assessed = n(), n_binge = sum(binge_any), prevalence = mean(binge_any), mean_age = mean(age), .groups = "drop")
  safe_write_csv(bins, "28_binge_by_calendar_halfyear", count_cols = c("n_assessed", "n_binge"))
}

## T4 onset hazard with time-varying pandemic era (counting process)
ev <- load_derived("events") |> filter(!is.na(cal_bl), !is.na(edu_bl)) |>
  mutate(age_bl_c = age_bl - 10, cohort_c = birth_year - mean(birth_year, na.rm = TRUE),
         s0 = pmax(0, c0 - cal_bl), s1 = pmax(0, c1 - cal_bl))
split_era <- function(ev, tcol = "t_event") {
  base <- ev |> mutate(tt = .data[[tcol]])
  tm <- survival::tmerge(base |> select(-any_of(c("t_left", "t_right"))), base, id = id, onset = event(tt, event))
  tm <- survival::tmerge(tm, base |> filter(s0 > 0) |> transmute(id, s0, v = 1L), id = id, era_r = tdc(s0, v))
  tm <- survival::tmerge(tm, base |> filter(s1 > 0) |> transmute(id, s1, v = 1L), id = id, era_p = tdc(s1, v))
  tm |> mutate(era_r = ifelse(is.na(era_r), 0L, era_r), era_p = ifelse(is.na(era_p), 0L, era_p),
               era = factor(ifelse(era_p == 1, "post", ifelse(era_r == 1, "restricted", "pre")), levels = c("pre", "restricted", "post")))
}
tm <- try_fit("T4 split", split_era(ev))
if (!is.null(tm)) {
  t4 <- try_fit("T4 Cox era", survival::coxph(survival::Surv(tstart, tstop, onset) ~ era + cohort_c + female + age_bl_c + edu_bl, data = tm, cluster = fam))
  if (!is.null(t4)) out$T4 <- broom::tidy(t4, conf.int = TRUE) |>
      transmute(module = "28", model = "T4_cox_timevarying_era", outcome = "binge_onset", term, estimate, std.error = robust.se,
                conf.low, conf.high, p.value, n_id = nrow(ev), n_events = sum(ev$event))
}

## T5 joint model with time-varying era in the survival submodel
if (!is.null(tm) && need(c("JMbayes2", "nlme"))) {
  mz <- paste0(cfg$primary, "_z")
  dl <- d |> filter(id %in% tm$id) |>
    filter(time <= t_event, !is.na(.data[[mz]])) |> select(id, time, female, age_bl_c, all_of(mz)) |> arrange(id, time) |> as.data.frame()
  tmj <- tm |> filter(id %in% dl$id) |> arrange(id, tstart) |> as.data.frame()
  lme5 <- try_fit("T5 lme", nlme::lme(as.formula(paste(mz, "~ time + female + age_bl_c")), random = ~ time | id, data = dl,
                                      control = nlme::lmeControl(opt = "optim")))
  cox5 <- try_fit("T5 cox", survival::coxph(survival::Surv(tstart, tstop, onset) ~ era + cohort_c + female + age_bl_c + edu_bl, data = tmj, x = TRUE))
  if (!is.null(lme5) && !is.null(cox5)) {
    j5 <- try_fit("T5 joint model", JMbayes2::jm(cox5, lme5, time_var = "time", data_Surv = tmj, id_var = "id",
                                                  control = list(n_iter = cfg$jm_iter, n_burnin = cfg$jm_burnin, n_chains = 3L, cores = 1L)))
    if (!is.null(j5)) {
      sv <- as.data.frame(summary(j5)$Survival); sv$term <- rownames(sv)
      out$T5 <- as_tibble(sv) |> transmute(module = "28", model = "T5_joint_timevarying_era", outcome = "binge_onset", term,
                                           estimate = Mean, std.error = StDev, conf.low = `2.5%`, conf.high = `97.5%`, p.value = P,
                                           Rhat, n_id = n_distinct(tmj$id), n_events = sum(tmj$onset))
    }
  }
}
res <- bind_rows(out)
safe_write_csv(res, "28_temporal_trends", count_cols = c("n_id", "n_events"))
if (nrow(res)) key_write(res |> filter(term %in% c("pand_yrs", "cohort_c") & model == "T2c_age_cohort_pandemic" & outcome == cfg$primary |
                                       grepl("^era|covid_era|remote", term) | grepl("value\\(", term)) |>
                           select(model, outcome, term, estimate, conf.low, conf.high, p.value, n_id), "28")
