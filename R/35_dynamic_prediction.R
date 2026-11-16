# =============================================================================
# 35_dynamic_prediction.R : individualized binge-onset risk from the joint model
# =============================================================================
# The joint model of 24 (value association) gives each youth a predicted risk
# of onset in (t_L, t_L + dt] given the scans observed up to the landmark t_L
# (cfg$dp_tstart years since baseline; horizon cfg$dp_dt). Discrimination
# (time-dependent AUC) and calibration/accuracy (time-dependent Brier score)
# are estimated out of sample by K-fold cross-validation with folds formed by
# family, so siblings never sit in both training and test sets. The same
# metrics for a covariates-only Cox model show what the brain marker adds.
lp_log("35_dynamic_prediction")
if (!need(c("JMbayes2", "nlme", "survival"))) stop("JMbayes2, nlme, survival required")
mz <- paste0(cfg$primary, "_z")
ev <- load_derived("events") |> filter(!is.na(edu_bl)) |> mutate(age_bl_c = age_bl - 10)
d <- img_data() |> filter(id %in% ev$id, !is.na(.data[[mz]])) |>
  select(id, time, female, age_bl_c, all_of(mz)) |> inner_join(ev |> select(id, fam, t_event, event, edu_bl), by = "id") |>
  filter(time <= t_event) |> arrange(id, time) |> as.data.frame()
ev <- ev |> filter(id %in% d$id) |> arrange(id) |> as.data.frame()
K <- if (isTRUE(cfg$fast)) 2L else cfg$cv_folds
fams <- unique(ev$fam); fold_of <- setNames(sample(rep(seq_len(K), length.out = length(fams))), fams)
ev$fold <- fold_of[ev$fam]; d$fold <- fold_of[d$fam]
ctrl <- list(n_iter = cfg$jm_iter, n_burnin = cfg$jm_burnin, n_chains = 3L, cores = 1L)
res <- list()
for (k in seq_len(K)) {
  tr_l <- d[d$fold != k, ]; te_l <- d[d$fold == k, ]; tr_s <- ev[ev$fold != k, ]; te_s <- ev[ev$fold == k, ]
  lm_k <- try_fit(paste("fold", k, "lme"), nlme::lme(as.formula(paste(mz, "~ time + female + age_bl_c")), random = ~ time | id,
                                                     data = tr_l, control = nlme::lmeControl(opt = "optim")))
  cx_k <- try_fit(paste("fold", k, "cox"), survival::coxph(survival::Surv(t_event, event) ~ female + age_bl_c + edu_bl, data = tr_s, x = TRUE))
  if (is.null(lm_k) || is.null(cx_k)) next
  j_k <- try_fit(paste("fold", k, "jm"), JMbayes2::jm(cx_k, lm_k, time_var = "time", data_Surv = tr_s, id_var = "id", control = ctrl))
  if (is.null(j_k)) next
  auc <- tryCatch(JMbayes2::tvAUC(j_k, te_l, Tstart = cfg$dp_tstart, Dt = cfg$dp_dt)$auc, error = function(e) { lp_log("  tvAUC: ", conditionMessage(e)); NA })
  bri <- tryCatch(JMbayes2::tvBrier(j_k, te_l, Tstart = cfg$dp_tstart, Dt = cfg$dp_dt)$Brier, error = function(e) { lp_log("  tvBrier: ", conditionMessage(e)); NA })
  auc0 <- tryCatch(JMbayes2::tvAUC(cx_k, te_s, Tstart = cfg$dp_tstart, Dt = cfg$dp_dt)$auc, error = function(e) NA)
  bri0 <- tryCatch(JMbayes2::tvBrier(cx_k, te_s, Tstart = cfg$dp_tstart, Dt = cfg$dp_dt)$Brier, error = function(e) NA)
  at_risk <- sum(te_s$t_event > cfg$dp_tstart)
  ev_win <- sum(te_s$event == 1 & te_s$t_event > cfg$dp_tstart & te_s$t_event <= cfg$dp_tstart + cfg$dp_dt)
  res[[k]] <- tibble(module = "35", fold = k, model = c("joint_model_marker", "cox_covariates_only"),
                     AUC = c(auc, auc0), Brier = c(bri, bri0), n_at_risk = at_risk, n_events_window = ev_win)
}
out <- bind_rows(res)
safe_write_csv(out, "35_dynamic_prediction_folds", count_cols = c("n_at_risk", "n_events_window"))
if (nrow(out)) {
  sm <- out |> group_by(model) |> summarise(AUC = mean(AUC, na.rm = TRUE), Brier = mean(Brier, na.rm = TRUE),
                                           n_folds = n(), n_events_window = sum(n_events_window), .groups = "drop") |>
    mutate(module = "35", landmark_yr = cfg$dp_tstart, horizon_yr = cfg$dp_dt)
  safe_write_csv(sm, "35_dynamic_prediction", count_cols = "n_events_window")
  key_write(sm |> transmute(model, term = sprintf("AUC (landmark %g y, horizon %g y)", cfg$dp_tstart, cfg$dp_dt), estimate = AUC, n_events = n_events_window), "35")
}
