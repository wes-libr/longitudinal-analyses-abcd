# =============================================================================
# 34_regional.R : regional (Desikan ROI) onset-anchored models with FDR, and a
# data-driven composite marker for the joint model
# =============================================================================
# R1 Mass-univariate: the piecewise model of 25 for every cortical-thickness
#    ROI in cfg$roi_table (68 in ABCD 7.0), Benjamini-Hochberg FDR across ROIs
#    separately for each onset term (eventual, eventual x time, post, tso).
# R2 Composite marker: PCA of baseline ROI thickness (each ROI standardized by
#    its baseline mean/SD); the PC1 loadings are applied to every scan, giving
#    a longitudinal composite on the baseline-SD scale. The composite then
#    enters the joint model of 24 (value association) as a low-dimensional
#    summary of regional development.
# Vertex-wise analysis (FEMA, BLMM) needs vertex data outside the tabulated
# release; the same design (eventual, post, tso, nested random effects) is
# what to specify there.
lp_log("34_regional")
if (!need("lme4")) stop("lme4 required")
rn <- load_derived("roi_names")$roi
d <- img_data()
rn <- rn[rn %in% names(d)]
if (!length(rn)) stop("No ROI columns (cfg$roi_table/roi_pattern matched nothing)")
fx <- list()
for (r in rn) {
  dm <- d |> filter(!is.na(.data[[r]]))
  f <- paste0("`", r, "` ~ time + eventual + eventual:time + post + tso + age_bl_c + female + ", re_terms("time"))
  fit <- try_fit(paste("ROI", r), lme4::lmer(as.formula(f), data = dm, control = lme4::lmerControl(calc.derivs = FALSE)))
  if (is.null(fit)) next
  if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
  fx[[r]] <- tidy_fixed(fit, "34", "R1_roi_piecewise", r) |> filter(term %in% c("eventual", "time:eventual", "eventual:time", "post", "tso")) |>
    mutate(n_id = n_distinct(dm$id))
}
roi <- bind_rows(fx) |> group_by(term) |> mutate(p.fdr = p.adjust(p.value, "BH"), fdr_sig = p.fdr < 0.05) |> ungroup()
safe_write_csv(roi, "34_roi_onset_anchored", count_cols = "n_id")
safe_write_csv(roi |> group_by(term) |> summarise(n_rois = n(), n_fdr_sig = sum(fdr_sig), .groups = "drop") |> mutate(module = "34"),
               "34_roi_fdr_summary")
key_write(roi |> filter(fdr_sig) |> transmute(model = "R1_roi_piecewise_fdr", outcome, term, estimate, conf.low, conf.high, p.value = p.fdr, n_id), "34")

## R2 composite marker
bl <- d |> filter(is_bl) |> select(id, all_of(rn)) |> filter(if_all(all_of(rn), ~ !is.na(.x)))
mu <- sapply(bl[rn], mean); sdv <- sapply(bl[rn], sd)
pc <- prcomp(scale(bl[rn], center = mu, scale = sdv))
load1 <- pc$rotation[, 1] * sign(sum(pc$rotation[, 1]))           # orient: higher = thicker overall
zall <- scale(as.matrix(d[rn]), center = mu, scale = sdv)
d$composite <- as.numeric(zall %*% load1)
d$composite <- (d$composite - mean(d$composite[d$is_bl], na.rm = TRUE)) / sd(d$composite[d$is_bl], na.rm = TRUE)
safe_write_csv(tibble(module = "34", roi = rn, pc1_loading = load1, var_explained_pc1 = pc$sdev[1]^2 / sum(pc$sdev^2)), "34_composite_loadings")
save_derived(d |> select(id, time, composite), "roi_composite")
if (need(c("JMbayes2", "nlme", "survival"))) {
  ev <- load_derived("events") |> filter(!is.na(edu_bl)) |> mutate(age_bl_c = age_bl - 10)
  dl <- d |> filter(id %in% ev$id, !is.na(composite)) |> left_join(ev |> select(id, t_ev = t_event), by = "id") |>
    filter(time <= t_ev) |> select(id, time, female, age_bl_c, composite) |> arrange(id, time) |> as.data.frame()
  ev <- ev |> filter(id %in% dl$id) |> arrange(id) |> as.data.frame()
  lm1 <- try_fit("lme composite", nlme::lme(composite ~ time + female + age_bl_c, random = ~ time | id, data = dl, control = nlme::lmeControl(opt = "optim")))
  cx1 <- try_fit("cox composite", survival::coxph(survival::Surv(t_event, event) ~ female + age_bl_c + edu_bl, data = ev, x = TRUE))
  if (!is.null(lm1) && !is.null(cx1)) {
    j <- try_fit("JM composite", JMbayes2::jm(cx1, lm1, time_var = "time", data_Surv = ev, id_var = "id",
                                              control = list(n_iter = cfg$jm_iter, n_burnin = cfg$jm_burnin, n_chains = 3L, cores = 1L)))
    if (!is.null(j)) {
      sv <- as.data.frame(summary(j)$Survival); sv$term <- rownames(sv)
      tb <- as_tibble(sv) |> transmute(module = "34", model = "R2_joint_model_roi_composite", term, estimate = Mean,
                                       conf.low = `2.5%`, conf.high = `97.5%`, p.value = P, Rhat, n_id = nrow(ev), n_events = sum(ev$event))
      safe_write_csv(tb, "34_composite_joint_model", count_cols = c("n_id", "n_events"))
      key_write(tb |> filter(grepl("value", term)) |> select(model, term, estimate, conf.low, conf.high, p.value, n_id, n_events), "34b")
    }
  }
}
