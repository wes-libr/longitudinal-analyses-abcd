# =============================================================================
# 30_target_trial.R : target-trial emulation of binge initiation -> brain change
# =============================================================================
# Target trial. Eligible: youth with no binge episode to date, a QC-passed scan
# at or before the trial's start, and an annual assessment at start k.
# Strategies: initiate binge drinking by the next annual assessment (A = 1) vs
# not (A = 0). Outcome: annualized change in the marker from the last scan
# before the trial to the first scan after the next assessment (within 3 years).
# A sequence of nested trials (one per annual visit k = 0..5) is stacked; a
# youth can enter several trials until initiation. Family-clustered inference.
#
# Confounders at trial start L_k: last marker value (z), CBCL externalizing,
# pubertal stage, sipping (lifetime alcohol indicator), age, sex, caregiver
# education, trial number.
# Estimators:
#   IPTW : stabilized treatment weights x inverse-probability-of-censoring
#          weights for a missing outcome scan; weighted GEE (marginal structural
#          model for the point treatment in each trial).
#   g-formula : outcome regression E[Y | A, L] standardized over L for A = 1
#          and A = 0; family-cluster bootstrap (cfg$tt_boot replicates).
# Truth in tests/mock_nbdc.R: extra thinning of 0.010 mm/yr after onset.
lp_log("30_target_trial")
if (!need(c("geepack"))) stop("geepack required")
m <- cfg$primary
ann <- load_derived("panel_annual") |> arrange(id, wave_yr)
img <- load_derived("panel_img") |> select(id, time, y = all_of(m)) |> filter(!is.na(y)) |> arrange(id, time)
base <- load_derived("baseline")
trials <- list()
for (k in 0:5) {
  a_k <- ann |> filter(wave_yr == k) |> select(id, t0 = time, ext, pds, alc_lt, age)
  a_k1 <- ann |> filter(wave_yr == k + 1) |> select(id, t1 = time, A = binge_any)
  hist_ok <- ann |> filter(wave_yr <= k) |> group_by(id) |> summarise(clean = all(binge_any %in% 0), .groups = "drop") |> filter(clean)
  tk <- a_k |> inner_join(hist_ok, by = "id") |> inner_join(a_k1, by = "id") |> filter(!is.na(A), !is.na(t0), !is.na(t1))
  pre <- img |> inner_join(tk |> select(id, t0), by = "id") |> filter(time <= t0 + 0.25) |>
    group_by(id) |> slice_max(time, n = 1, with_ties = FALSE) |> ungroup() |> transmute(id, t_pre = time, y_pre = y)
  postscan <- img |> inner_join(tk |> select(id, t1), by = "id") |> filter(time > t1, time <= t1 + 3) |>
    group_by(id) |> slice_min(time, n = 1, with_ties = FALSE) |> ungroup() |> transmute(id, t_post = time, y_post = y)
  trials[[k + 1]] <- tk |> inner_join(pre, by = "id") |> left_join(postscan, by = "id") |>
    mutate(trial = k, Y = (y_post - y_pre) / (t_post - t_pre))
}
tt <- bind_rows(trials) |> left_join(base |> select(id, fam, site, female, edu_bl), by = "id") |>
  mutate(fam = paste(site, fam, sep = "_"), y_pre_z = as.numeric(scale(y_pre)), observed = as.integer(!is.na(Y)),
         ext = ifelse(is.na(ext), mean(ext, na.rm = TRUE), ext), pds = ifelse(is.na(pds), mean(pds, na.rm = TRUE), pds),
         alc_lt = ifelse(is.na(alc_lt), 0, alc_lt)) |>
  filter(!is.na(edu_bl), !is.na(age)) |> arrange(fam, id, trial)
save_derived(tt, "target_trial_data")
safe_write_csv(tt |> group_by(trial) |> summarise(n_eligible = n(), n_initiate = sum(A), n_outcome = sum(observed), .groups = "drop"),
               "30_trial_flow", count_cols = c("n_eligible", "n_initiate", "n_outcome"))
Lf <- "y_pre_z + ext + pds + alc_lt + age + female + edu_bl + trial"
res <- list()
# IPTW + IPCW
ps <- glm(as.formula(paste("A ~", Lf)), family = binomial, data = tt)
pA <- mean(tt$A); p <- predict(ps, type = "response")
sw <- ifelse(tt$A == 1, pA / p, (1 - pA) / (1 - p))
cm <- glm(as.formula(paste("observed ~ A +", Lf)), family = binomial, data = tt)
cw <- mean(tt$observed) / predict(cm, type = "response")
tt$w <- sw * cw; q <- quantile(tt$w, c(.01, .99)); tt$w <- pmin(pmax(tt$w, q[1]), q[2])
to <- tt |> filter(observed == 1)
g_un <- try_fit("unadjusted", geepack::geeglm(Y ~ A, id = factor(fam, levels = unique(fam)), data = to, corstr = "independence"))
g_ip <- try_fit("IPTW MSM", geepack::geeglm(Y ~ A, id = factor(fam, levels = unique(fam)), weights = w, data = to, corstr = "independence"))
for (g in list(list("unadjusted", g_un), list("IPTW_IPCW_msm", g_ip))) if (!is.null(g[[2]])) {
  s <- summary(g[[2]])$coefficients["A", ]
  res[[g[[1]]]] <- tibble(module = "30", model = g[[1]], term = "initiate_binge (mm/yr)", estimate = s$Estimate, std.error = s$Std.err,
                          conf.low = s$Estimate - 1.96 * s$Std.err, conf.high = s$Estimate + 1.96 * s$Std.err)
}
# parametric g-formula with family bootstrap
gform <- function(dat) {
  om <- lm(as.formula(paste("Y ~ A * (", Lf, ")")), data = dat |> filter(observed == 1))
  mean(predict(om, dat |> mutate(A = 1))) - mean(predict(om, dat |> mutate(A = 0)))
}
est <- gform(tt)
fams <- unique(tt$fam); B <- if (isTRUE(cfg$fast)) 30L else cfg$tt_boot
sp <- split(seq_len(nrow(tt)), tt$fam)
boots <- vapply(seq_len(B), function(b) { f <- sample(fams, replace = TRUE); tryCatch(gform(tt[unlist(sp[f]), ]), error = function(e) NA_real_) }, 0)
res$gf <- tibble(module = "30", model = "g_formula", term = "initiate_binge (mm/yr)", estimate = est, std.error = sd(boots, na.rm = TRUE),
                 conf.low = quantile(boots, .025, na.rm = TRUE), conf.high = quantile(boots, .975, na.rm = TRUE))
out <- bind_rows(res) |> mutate(p.value = 2 * pnorm(-abs(estimate / std.error)), n_trials = nrow(tt), n_id = n_distinct(tt$id),
                                n_initiators = sum(tt$A), max_weight = max(tt$w))
safe_write_csv(out, "30_target_trial", count_cols = c("n_id", "n_initiators"))
key_write(out |> select(model, term, estimate, conf.low, conf.high, p.value, n_id, n_initiators), "30")
