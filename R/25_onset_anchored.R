# =============================================================================
# 25_onset_anchored.R : does brain development change after binge onset?
# =============================================================================
# Piecewise ("interrupted") mixed model on all participants:
#   y ~ time + eventual + eventual:time      (pre-existing level/slope differences
#                                              of eventual bingers, BEFORE onset)
#       + post + tso                          (level shift and slope change AFTER
#                                              onset; tso = years since onset)
#       + age_bl_c + female [+ icv_z] + nested RE + scanner
# Never-bingers anchor the age trend; eventual:time absorbs selection on
# trajectory, so post/tso target change attributable to the onset period.
# Sensitivity: (S1) drop scans within 1 year after onset (reverse timing /
# measurement noise around onset), (S2) eventual bingers only.
lp_log("25_onset_anchored")
if (!need("lme4")) stop("lme4 required")
d <- img_data()
fx <- list()
for (m in names(cfg$brain)) {
  dm <- d |> filter(!is.na(.data[[m]]))
  f <- paste0(m, " ~ time + eventual + eventual:time + post + tso + age_bl_c + female", icv_term(m), " + ", re_terms("time"))
  sets <- list(main = dm, S1_drop_within_1yr = dm |> filter(!(post == 1 & tso < 1)), S2_bingers_only = dm |> filter(eventual == 1))
  for (s in names(sets)) {
    ff <- if (s == "S2_bingers_only") sub("eventual \\+ eventual:time \\+ ", "", f) else f
    fit <- try_fit(paste(m, s), lme4::lmer(as.formula(ff), data = sets[[s]], control = lme4::lmerControl(calc.derivs = FALSE)))
    if (is.null(fit)) next
    if (requireNamespace("lmerTest", quietly = TRUE)) fit <- if (isFALSE(getOption("lp.lmertest_ok"))) fit else tryCatch(lmerTest::as_lmerModLmerTest(fit), error = function(e) { lp_log("  lmerTest conversion failed (", conditionMessage(e), "); using lme4 fit, normal-approximation p-values"); fit })
    fx[[paste(m, s)]] <- tidy_fixed(fit, "25", s, m) |>
      mutate(n_id = n_distinct(sets[[s]]$id), n_bingers = n_distinct(sets[[s]]$id[sets[[s]]$eventual == 1]),
             n_post_scans = sum(sets[[s]]$post))
  }
}
out <- bind_rows(fx)
safe_write_csv(out, "25_onset_anchored", count_cols = c("n_id", "n_bingers", "n_post_scans"))
if (nrow(out)) key_write(out |> filter(term %in% c("eventual", "time:eventual", "eventual:time", "post", "tso")) |>
                           select(model, outcome, term, estimate, conf.low, conf.high, p.value, n_id, n_bingers), "25")
