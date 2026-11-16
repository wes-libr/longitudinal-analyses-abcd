# =============================================================================
# 24_joint_model.R : FLAGSHIP joint model of brain trajectories and binge onset
# =============================================================================
# Question: does a youth's (error-free) brain trajectory, its current level and
# its current rate of change, predict the hazard of binge-drinking onset?
#
# Longitudinal submodel (nlme::lme): marker_z ~ time + female + age_bl_c,
#   random intercept + slope by participant. Marker standardized by baseline SD,
#   so association parameters are per baseline SD.
# Survival submodel (coxph): time = years since baseline to the first binge-
#   positive assessment; adjusted for sex, baseline age, caregiver education.
# Association (JMbayes2): value(marker) and slope(marker) at time t.
# Only scans at or before onset enter the longitudinal submodel (scans after
# onset belong to the "consequences" question in 25_onset_anchored).
#
# Variants: J1 value; J2 value + slope; J3 multivariate (cfg$joint_multi,
# value of each); J4 = J1 with midpoint onset times (interval-censoring
# sensitivity). Comparison: naive time-varying Cox with the last observed
# (error-prone) scan carried forward, which the joint model should improve on
# through measurement-error correction.
#
# Limitations noted in the output: family/site clustering is not in the
# longitudinal submodel (lme takes one grouping level); onset is interval-
# censored between assessments (J4 brackets the effect of the convention).
lp_log("24_joint_model")
if (!need(c("JMbayes2", "nlme", "survival"))) stop("JMbayes2, nlme, survival required")
d <- img_data(); ev0 <- load_derived("events") |> filter(!is.na(edu_bl)) |>
  mutate(age_bl_c = age_bl - 10, t_mid = ifelse(event == 1, (t_left + t_right) / 2, t_right))
ctrl <- if (isTRUE(cfg$fast) || is.null(cfg$jm24_iter)) list(n_iter = cfg$jm_iter, n_burnin = cfg$jm_burnin, n_chains = 3L, cores = 1L) else
  list(n_iter = cfg$jm24_iter, n_burnin = cfg$jm24_burnin, n_thin = cfg$jm24_thin, n_chains = 3L,
       cores = max(1L, min(3L, as.integer(cfg$jm_cores), parallel::detectCores() - 1L)))
lp_log("  JMbayes2 control: ", paste(names(ctrl), unlist(ctrl), sep = "=", collapse = ", "))

prep <- function(markers, tcol = "t_event") {
  ev <- ev0 |> mutate(t_jm = .data[[tcol]])
  dl <- d |> filter(id %in% ev$id) |> left_join(ev |> select(id, t_jm), by = "id") |>
    filter(time <= t_jm) |> select(id, fam, time, female, age_bl_c, all_of(paste0(markers, "_z"))) |>
    filter(if_all(all_of(paste0(markers, "_z")), ~ !is.na(.x))) |> arrange(id, time) |> as.data.frame()
  ev <- ev |> filter(id %in% dl$id) |> arrange(id) |> as.data.frame()
  list(dl = dl, ev = ev)
}
fit_lme <- function(dl, mz) nlme::lme(as.formula(paste(mz, "~ time + female + age_bl_c")), random = ~ time | id,
                                      data = dl, control = nlme::lmeControl(opt = "optim"), na.action = na.omit)
fit_cox <- function(ev) survival::coxph(survival::Surv(t_jm, event) ~ female + age_bl_c + edu_bl, data = ev, x = TRUE)
jm_tab <- function(j, label, dat) {
  s <- summary(j)
  sv <- as.data.frame(s$Survival); sv$term <- rownames(sv)
  lg <- bind_rows(lapply(seq_along(s$Outcomes), function(k) { o <- as.data.frame(s$Outcomes[[k]]); o$term <- paste0("long", k, ": ", rownames(o)); o }))
  bind_rows(sv, lg) |> as_tibble() |>
    transmute(module = "24", model = label, term, estimate = Mean, sd = StDev, conf.low = `2.5%`, conf.high = `97.5%`,
              p.value = P, Rhat = Rhat, n_id = nrow(dat$ev), n_events = sum(dat$ev$event), n_scans = nrow(dat$dl)) |>
    (\(tb) { bad <- tb |> filter(!grepl("^long", term), Rhat > 1.05)
      if (nrow(bad)) lp_log("  WARNING ", label, ": Rhat > 1.05 for ", paste0(bad$term, " (", round(bad$Rhat, 3), ")", collapse = ", "),
                            "; increase cfg$jm24_iter") else lp_log("  ", label, ": all survival-submodel Rhat <= 1.05")
      tb })()
}
res <- list()
m1 <- cfg$primary; mz <- paste0(m1, "_z")
P <- prep(m1)
lme1 <- try_fit("lme", fit_lme(P$dl, mz)); cox1 <- try_fit("cox", fit_cox(P$ev))
if (!is.null(lme1) && !is.null(cox1)) {
  j1 <- try_fit("J1 value", JMbayes2::jm(cox1, lme1, time_var = "time", data_Surv = P$ev, id_var = "id", control = ctrl))
  if (!is.null(j1)) { res$J1 <- jm_tab(j1, "J1_value", P); save_derived(j1, "jm_J1") }
  ff <- as.formula(paste0("~ value(", mz, ") + slope(", mz, ")"))
  j2 <- try_fit("J2 value+slope", JMbayes2::jm(cox1, lme1, time_var = "time", functional_forms = ff, data_Surv = P$ev, id_var = "id", control = ctrl))
  if (!is.null(j2)) { res$J2 <- jm_tab(j2, "J2_value_slope", P); save_derived(j2, "jm_J2") }
}
# J3 multivariate
mm <- intersect(cfg$joint_multi, names(cfg$brain))
if (length(mm) >= 2) {
  P3 <- prep(mm)
  lmes <- lapply(paste0(mm, "_z"), function(v) try_fit(paste("lme", v), fit_lme(P3$dl, v)))
  cox3 <- try_fit("cox (J3)", fit_cox(P3$ev))
  if (all(!vapply(lmes, is.null, logical(1))) && !is.null(cox3)) {
    j3 <- try_fit("J3 multivariate", JMbayes2::jm(cox3, lmes, time_var = "time", data_Surv = P3$ev, id_var = "id", control = ctrl))
    if (!is.null(j3)) { res$J3 <- jm_tab(j3, "J3_multivariate_value", P3); save_derived(j3, "jm_J3") }
  }
}
# J4 midpoint onset
P4 <- prep(m1, "t_mid")
lme4_ <- try_fit("lme (J4)", fit_lme(P4$dl, mz)); cox4 <- try_fit("cox (J4)", fit_cox(P4$ev))
if (!is.null(lme4_) && !is.null(cox4)) {
  j4 <- try_fit("J4 midpoint", JMbayes2::jm(cox4, lme4_, time_var = "time", data_Surv = P4$ev, id_var = "id", control = ctrl))
  if (!is.null(j4)) res$J4 <- jm_tab(j4, "J4_value_midpoint_onset", P4)
}
# Naive comparison: time-varying Cox with last observed scan carried forward
base <- P$ev |> select(id, fam, female, age_bl_c, edu_bl, t_jm, event)
tm <- survival::tmerge(base, base, id = id, onset = event(t_jm, event))
tm <- survival::tmerge(tm, P$dl |> select(id, time, mk = all_of(mz)), id = id, mk = tdc(time, mk))
nv <- try_fit("naive TV Cox", survival::coxph(survival::Surv(tstart, tstop, onset) ~ mk + female + age_bl_c + edu_bl, data = tm, cluster = fam))
if (!is.null(nv)) res$naive <- broom::tidy(nv, conf.int = TRUE) |>
  transmute(module = "24", model = "naive_tv_cox_locf", term = ifelse(term == "mk", paste0("value(", mz, ")"), term),
            estimate, sd = std.error, conf.low, conf.high, p.value, n_id = nrow(base), n_events = sum(base$event))
out <- bind_rows(res) |> mutate(scale = "log hazard ratio per baseline SD of the marker (association rows)")
safe_write_csv(out, "24_joint_models", count_cols = c("n_id", "n_events"))
if (nrow(out)) key_write(out |> filter(grepl("value\\(|slope\\(", term)) |> select(model, term, estimate, conf.low, conf.high, p.value, any_of("Rhat"), n_id, n_events), "24")
