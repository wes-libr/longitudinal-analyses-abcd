# =============================================================================
# 10_change.R : difference scores and residualized change (tutorial methods)
# =============================================================================
# Difference score: paired change baseline -> Year 2 per marker.
# Residualized change: last scan ~ baseline scan + years between scans +
#   binge onset before the last scan, family-clustered robust SEs (GEE, identity).
lp_log("10_change")
d <- img_data()
out <- list(); res <- list()
for (m in names(cfg$brain)) {
  w <- d |> filter(wave_yr %in% c(0, 2)) |> select(id, wave_yr, y = all_of(m)) |>
    pivot_wider(names_from = wave_yr, values_from = y, names_prefix = "w") |> filter(!is.na(w0), !is.na(w2))
  if (nrow(w) < 20) next
  tt <- t.test(w$w2, w$w0, paired = TRUE)
  out[[m]] <- tibble(module = "10", model = "paired_difference_bl_to_y2", marker = m, n_id = nrow(w),
                     mean_change = unname(tt$estimate), conf.low = tt$conf.int[1], conf.high = tt$conf.int[2],
                     p.value = tt$p.value, cohen_dz = mean(w$w2 - w$w0) / sd(w$w2 - w$w0))

  last <- d |> filter(!is.na(.data[[m]])) |> group_by(id) |> filter(n() >= 2) |>
    summarise(fam = first(fam), y_bl = .data[[m]][which.min(time)], y_last = .data[[m]][which.max(time)],
              gap = max(time) - min(time), icv_z = icv_z[which.min(time)], female = first(female),
              age_bl_c = first(age_bl_c), t_last = max(time), event = first(event), t_event = first(t_event), .groups = "drop") |>
    mutate(binge_before_last = as.integer(!is.na(event) & event == 1 & t_event < t_last)) |>
    arrange(fam)
  if (need("geepack")) {
    f <- as.formula(paste0("y_last ~ y_bl + gap + binge_before_last + female + age_bl_c", icv_term(m)))
    g <- try_fit(paste("residualized change", m), geepack::geeglm(f, id = factor(fam, levels = unique(fam)), data = last, corstr = "independence"))
    if (!is.null(g)) res[[m]] <- tidy_fixed(g, "10", "residualized_change", m) |>
      mutate(n_id = nrow(last), n_binge = sum(last$binge_before_last))
  }
}
safe_write_csv(bind_rows(out), "10_difference_scores", count_cols = "n_id")
safe_write_csv(bind_rows(res), "10_residualized_change", count_cols = c("n_id", "n_binge"))
