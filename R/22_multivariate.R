# =============================================================================
# 22_multivariate.R : coupled development of several brain markers
# =============================================================================
#  MV-LMM   stacked multivariate mixed model (glmmTMB): marker-specific fixed
#           effects and residual variances; participant random intercepts and
#           slopes for all markers in one unstructured covariance. Output: the
#           correlation matrix of intercepts/slopes across markers (do regions
#           that start thinner/smaller also change faster? do regions co-develop?).
#  PP-LGCM  parallel-process latent growth (lavaan) for the first two markers,
#           with eventual binge onset predicting all growth factors.
lp_log("22_multivariate")
mk <- unique(c(cfg$joint_multi, "area"))
mk <- intersect(mk, names(cfg$brain))
d <- img_data()
long <- d |> select(id, fam, site, time, female, age_bl_c, icv_z, eventual, all_of(paste0(mk, "_z"))) |>
  pivot_longer(all_of(paste0(mk, "_z")), names_to = "marker", values_to = "y") |>
  filter(!is.na(y)) |> mutate(marker = factor(sub("_z$", "", marker), levels = mk))
out <- list()
if (need("glmmTMB")) {
  f <- y ~ 0 + marker + marker:time + marker:female + marker:age_bl_c + marker:icv_z + marker:eventual +
    us(0 + marker + marker:time | id)
  fit <- try_fit("MV-LMM", glmmTMB::glmmTMB(f, dispformula = ~ 0 + marker, data = long))
  if (!is.null(fit)) {
    save_derived(fit, "mvlmm")
    vc <- glmmTMB::VarCorr(fit)$cond$id
    R <- attr(vc, "correlation"); sds <- attr(vc, "stddev")
    rn <- gsub("marker", "", rownames(R)); rn <- ifelse(grepl(":time", rn), paste0(sub(":time", "", rn), "_slope"), paste0(rn, "_int"))
    dimnames(R) <- list(rn, rn)
    out$cor <- as_tibble(as.table(R), .name_repair = "minimal") |> setNames(c("re1", "re2", "correlation")) |>
      filter(as.integer(factor(re1, levels = rn)) < as.integer(factor(re2, levels = rn))) |>
      mutate(module = "22", model = "mvlmm_random_effect_correlations", n_id = n_distinct(long$id))
    out$sd <- tibble(module = "22", model = "mvlmm_random_effect_sd", re1 = rn, sd = sds, n_id = n_distinct(long$id))
    out$fx <- tidy_fixed(fit, "22", "mvlmm_fixed", "stacked_z") |> mutate(n_id = n_distinct(long$id))
    key_write(out$cor |> filter(grepl("_slope", re1) & grepl("_slope", re2)) |>
                transmute(model, term = paste(re1, re2, sep = " ~ "), estimate = correlation, n_id), "22")
  }
}
safe_write_csv(bind_rows(out$cor, out$sd), "22_mvlmm_random_effects", count_cols = "n_id")
if (!is.null(out$fx)) safe_write_csv(out$fx, "22_mvlmm_fixed", count_cols = "n_id")

if (need("lavaan") && length(mk) >= 2) {
  w <- img_wide(d |> mutate(y = .data[[paste0(mk[1], "_z")]], x = .data[[paste0(mk[2], "_z")]]), c("y", "x")) |>
    left_join(d |> distinct(id, eventual), by = "id")
  pp <- "iy =~ 1*y_0 + 1*y_2 + 1*y_4 + 1*y_6
         sy =~ 0*y_0 + 2*y_2 + 4*y_4 + 6*y_6
         ix =~ 1*x_0 + 1*x_2 + 1*x_4 + 1*x_6
         sx =~ 0*x_0 + 2*x_2 + 4*x_4 + 6*x_6
         iy + sy + ix + sx ~ eventual + female + age_bl"
  g <- try_fit("PP-LGCM", lavaan::growth(pp, data = w, missing = "fiml", estimator = "MLR", cluster = "fam"))
  if (!is.null(g)) {
    ss <- lavaan::standardizedSolution(g) |> as_tibble() |> filter(op %in% c("~~", "~"), lhs != rhs | op == "~") |>
      filter(lhs %in% c("iy", "sy", "ix", "sx")) |>
      transmute(module = "22", model = "parallel_process_lgcm", term = paste(lhs, op, rhs), std_estimate = est.std,
                conf.low = ci.lower, conf.high = ci.upper, p.value = pvalue, n_id = lavaan::lavInspect(g, "nobs"),
                note = paste0("y = ", mk[1], ", x = ", mk[2]))
    safe_write_csv(ss, "22_parallel_process_lgcm", count_cols = "n_id")
  }
}
