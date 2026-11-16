# =============================================================================
# 26_missing_data.R : missing imaging data under MAR and MNAR
# =============================================================================
# Estimands (primary marker, baseline-SD units, participants with a baseline
# scan):  E1 = mean annual change;  E2 = difference in annual change between
# eventual bingers and never-bingers (eventual:time).
#
#  A  MAR, likelihood       LMM on observed data (reference)
#  B  MAR, IPW              inverse-probability-of-observation weights per
#                           follow-up wave (baseline marker, externalizing, sex,
#                           age, education, eventual onset), stabilized and
#                           truncated at the 99th percentile; weighted GEE,
#                           family-clustered
#  C  MAR, MI               multiple imputation (mice, PMM) of the wide panel
#                           with externalizing at each wave as auxiliary data
#  D  MNAR, delta           delta-adjusted MI: imputed scans shifted by delta
#                           baseline SDs (all dropouts, or eventual bingers
#                           only); the tipping point is the delta at which E2's
#                           interval first includes 0
#  E  MNAR, pattern-mixture patterns = last observed imaging wave; pattern-
#                           specific slopes averaged over pattern shares; slopes
#                           for baseline-only participants identified by
#                           complete-case (CCMV) or neighbouring-case (NCMV)
#                           restrictions
# The joint model in 24 is the shared-parameter (selection-type) counterpart.
lp_log("26_missing_data")
if (!need(c("lme4", "geepack", "mice"))) stop("lme4, geepack, mice required")
d <- img_data(); m <- cfg$primary; mz <- paste0(m, "_z"); waves <- c(0, 2, 4, 6)
ann <- load_derived("panel_annual") |> filter(wave_yr %in% waves) |> select(id, wave_yr, ext)
bl_ids <- unique(d$id[d$wave_yr == 0 & !is.na(d[[mz]])])
long <- d |> filter(id %in% bl_ids, wave_yr %in% waves, !is.na(.data[[mz]])) |>
  transmute(id, fam, wave_yr, time = wave_yr, y = .data[[mz]], female, age_bl_c, edu_bl, eventual)
covs <- long |> group_by(id) |> summarise(across(c(fam, female, age_bl_c, edu_bl, eventual), first), .groups = "drop")
wide <- covs |>
  left_join(long |> select(id, wave_yr, y) |> pivot_wider(names_from = wave_yr, values_from = y, names_prefix = "y_"), by = "id") |>
  left_join(ann |> filter(id %in% bl_ids) |> pivot_wider(names_from = wave_yr, values_from = ext, names_prefix = "ext_"), by = "id")
for (v in c(paste0("y_", waves), paste0("ext_", waves))) if (!v %in% names(wide)) wide[[v]] <- NA_real_
wide <- wide |> filter(!is.na(edu_bl))
res <- list(); n_id <- nrow(wide)
grab <- function(est, se, method, estimand, extra = list())
  tibble(module = "26", model = method, term = estimand, estimate = est, std.error = se,
         conf.low = est - 1.96 * se, conf.high = est + 1.96 * se, n_id = n_id, !!!extra)
lmm_fit <- function(dl) {
  f <- lme4::lmer(y ~ time * eventual + female + age_bl_c + (1 + time | id), data = dl, control = lme4::lmerControl(calc.derivs = FALSE))
  b <- lme4::fixef(f); V <- as.matrix(vcov(f))
  list(E1 = c(b["time"], sqrt(V["time", "time"])), E2 = c(b["time:eventual"], sqrt(V["time:eventual", "time:eventual"])))
}
to_long <- function(w) w |> select(id, fam, female, age_bl_c, eventual, starts_with("y_")) |>
  pivot_longer(starts_with("y_"), names_to = "w", values_to = "y") |> mutate(time = as.numeric(sub("y_", "", w))) |> filter(!is.na(y))

## A. MAR likelihood
A <- try_fit("A LMM", lmm_fit(to_long(wide)))
if (!is.null(A)) res$A <- bind_rows(grab(A$E1[1], A$E1[2], "A_MAR_lmm", "E1_mean_slope"), grab(A$E2[1], A$E2[2], "A_MAR_lmm", "E2_slope_diff_eventual"))

## B. IPW
ipw <- wide
for (w in waves[-1]) {
  r <- as.integer(!is.na(ipw[[paste0("y_", w)]]))
  dat <- ipw |> mutate(r = r, ext0 = ifelse(is.na(ext_0), mean(ext_0, na.rm = TRUE), ext_0))
  g <- glm(r ~ y_0 + ext0 + female + age_bl_c + edu_bl + eventual, family = binomial, data = dat)
  p <- predict(g, type = "response")
  sw <- mean(r) / p; sw <- pmin(sw, quantile(sw, .99))
  ipw[[paste0("w_", w)]] <- sw
}
ipw$w_0 <- 1
lw <- to_long(ipw) |> left_join(ipw |> select(id, starts_with("w_")) |> pivot_longer(-id, names_to = "w", values_to = "wt") |>
                                  mutate(time = as.numeric(sub("w_", "", w))) |> select(-w), by = c("id", "time")) |>
  arrange(fam, id, time)
for (lab in c("unweighted", "ipw")) {
  g <- try_fit(paste("B GEE", lab), geepack::geeglm(y ~ time * eventual + female + age_bl_c, id = factor(fam, levels = unique(fam)),
                                                    weights = if (lab == "ipw") wt else rep(1, nrow(lw)), data = lw, corstr = "independence"))
  if (is.null(g)) next
  s <- summary(g)$coefficients
  res[[paste0("B_", lab)]] <- bind_rows(grab(s["time", "Estimate"], s["time", "Std.err"], paste0("B_gee_", lab), "E1_mean_slope"),
                                        grab(s["time:eventual", "Estimate"], s["time:eventual", "Std.err"], paste0("B_gee_", lab), "E2_slope_diff_eventual"))
}
safe_write_csv(tibble(module = "26", wave = waves[-1], mean_weight = sapply(waves[-1], function(w) mean(ipw[[paste0("w_", w)]])),
                      max_weight = sapply(waves[-1], function(w) max(ipw[[paste0("w_", w)]])), n_id = n_id),
               "26_ipw_weights", count_cols = "n_id")

## C/D. MI and delta adjustment
imp_vars <- c(paste0("y_", waves), paste0("ext_", waves), "female", "age_bl_c", "edu_bl", "eventual")
mi_in <- wide |> select(all_of(imp_vars)) |> as.data.frame()
miss <- is.na(mi_in[, paste0("y_", waves)])
imp <- try_fit("C mice", mice::mice(mi_in, m = cfg$mi_m, method = "pmm", printFlag = FALSE, seed = cfg$seed))
rubin <- function(q, u) { m <- length(q); qb <- mean(q); ub <- mean(u); b <- var(q); t <- ub + (1 + 1 / m) * b
  df <- (m - 1) * (1 + ub / ((1 + 1 / m) * b))^2; c(est = qb, se = sqrt(t), df = df) }
if (!is.null(imp)) {
  sdy <- 1   # markers are already in baseline-SD units
  scen <- tidyr::expand_grid(delta = cfg$delta_grid, who = c("all_dropouts", "eventual_bingers"))
  scen <- scen[!(scen$delta == 0 & scen$who == "eventual_bingers"), ]
  for (k in seq_len(nrow(scen))) {
    est <- lapply(seq_len(cfg$mi_m), function(i) {
      ci <- mice::complete(imp, i)
      shift <- miss * (scen$delta[k] * sdy)
      if (scen$who[k] == "eventual_bingers") shift[wide$eventual != 1, ] <- 0
      ci[, paste0("y_", waves)] <- ci[, paste0("y_", waves)] + shift
      lmm_fit(to_long(bind_cols(wide |> select(id, fam), ci)))
    })
    for (E in c("E1", "E2")) {
      rb <- rubin(sapply(est, function(e) e[[E]][1]), sapply(est, function(e) e[[E]][2]^2))
      lab <- if (scen$delta[k] == 0) "C_MI_MAR" else paste0("D_MI_delta_", scen$who[k])
      res[[paste(k, E)]] <- grab(rb[["est"]], rb[["se"]], lab, c(E1 = "E1_mean_slope", E2 = "E2_slope_diff_eventual")[[E]],
                                 list(delta_sd = scen$delta[k], df = rb[["df"]]))
    }
  }
}

## E. Pattern-mixture
pat <- long |> filter(id %in% wide$id) |> group_by(id) |> summarise(last_obs = max(wave_yr), .groups = "drop")
lp_long <- to_long(wide) |> left_join(pat, by = "id")
pmm_marginal <- function(dl, pw) {
  dl2 <- dl |> filter(last_obs > 0) |> mutate(pf = factor(last_obs))
  f <- lme4::lmer(y ~ time * pf + female + age_bl_c + (1 + time | id), data = dl2, control = lme4::lmerControl(calc.derivs = FALSE))
  b <- lme4::fixef(f); V <- as.matrix(vcov(f)); lv <- levels(dl2$pf)
  # slope per pattern = time + time:pf<level>
  L <- t(sapply(lv, function(l) { v <- setNames(rep(0, length(b)), names(b)); v["time"] <- 1
    if (l != lv[1]) v[paste0("time:pf", l)] <- 1; v }))
  slopes <- as.numeric(L %*% b); names(slopes) <- lv
  out <- list()
  for (restr in c("CCMV", "NCMV")) {
    ref <- if (restr == "CCMV") as.character(max(as.integer(lv))) else as.character(min(as.integer(lv)))
    wts <- setNames(rep(0, length(lv)), lv)
    for (p in names(pw)) { tgt <- if (p == "0") ref else p; if (tgt %in% lv) wts[tgt] <- wts[tgt] + pw[[p]] }
    cvec <- colSums(L * wts)
    out[[restr]] <- c(est = sum(cvec * b), se = sqrt(as.numeric(t(cvec) %*% V %*% cvec)))
  }
  list(marg = out, slopes = slopes)
}
pw_all <- prop.table(table(pat$last_obs))
pm <- try_fit("E PMM", pmm_marginal(lp_long, pw_all))
if (!is.null(pm)) for (r in names(pm$marg))
  res[[paste0("E_", r)]] <- grab(pm$marg[[r]][["est"]], pm$marg[[r]][["se"]], paste0("E_PMM_", r), "E1_mean_slope")
pm_by <- lapply(c(0, 1), function(e) { dd <- lp_long |> filter(eventual == e); pp <- prop.table(table(pat$last_obs[pat$id %in% dd$id]))
  try_fit(paste("E PMM eventual", e), pmm_marginal(dd, pp)) })
if (!any(vapply(pm_by, is.null, logical(1)))) for (r in c("CCMV", "NCMV")) {
  dlt <- pm_by[[2]]$marg[[r]][["est"]] - pm_by[[1]]$marg[[r]][["est"]]
  sev <- sqrt(pm_by[[2]]$marg[[r]][["se"]]^2 + pm_by[[1]]$marg[[r]][["se"]]^2)
  res[[paste0("E2_", r)]] <- grab(dlt, sev, paste0("E_PMM_", r), "E2_slope_diff_eventual")
}
out <- bind_rows(res) |> mutate(p.value = 2 * pnorm(-abs(estimate / std.error)))
safe_write_csv(out, "26_missing_data_sensitivity", count_cols = "n_id")
if (nrow(out)) {
  e2 <- out |> filter(term == "E2_slope_diff_eventual", grepl("^D_", model) | model == "C_MI_MAR")
  # Tipping point: if E2's interval excludes 0 under MAR (delta = 0), the
  # first delta at which it includes 0; if it already includes 0, the first
  # delta at which it excludes 0 (a difference that MNAR dropout could create).
  mar <- e2 |> filter(model == "C_MI_MAR")
  sig0 <- nrow(mar) > 0 && (mar$conf.low[1] > 0 | mar$conf.high[1] < 0)
  tip <- e2 |> filter(model != "C_MI_MAR") |> bind_rows(mar |> mutate(model = "D_MI_delta_all_dropouts"),
                                                         mar |> mutate(model = "D_MI_delta_eventual_bingers")) |>
    group_by(model) |> arrange(desc(delta_sd), .by_group = TRUE) |>
    summarise(significant_under_MAR = sig0,
              tipping_delta_sd = { cover <- conf.low <= 0 & conf.high >= 0; k <- if (sig0) which(cover) else which(!cover)
                                   if (length(k)) delta_sd[k[1]] else NA_real_ },
              .groups = "drop")
  safe_write_csv(tip, "26_tipping_points")
  key_write(out |> select(model, term, estimate, conf.low, conf.high, p.value, n_id, any_of("delta_sd")), "26")
}
