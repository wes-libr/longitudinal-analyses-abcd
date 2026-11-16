# =============================================================================
# 37_evalues.R : E-values for unmeasured confounding
# =============================================================================
# For each onset-related effect on a continuous brain outcome (25 onset-
# anchored, 30 target trial, 32 within-family), the estimate and its interval
# limit closest to the null are converted to a standardized mean difference
# d = estimate / SD(outcome at baseline) and then to an approximate risk ratio
# RR = exp(0.91 * |d|) (VanderWeele 2017). The E-value, RR + sqrt(RR (RR - 1)),
# is the minimum strength of association (risk-ratio scale) an unmeasured
# confounder would need with both onset and the outcome to explain away the
# estimate (or its interval). E-value of the interval = 1 when it includes 0.
# Hazard ratios (24 joint model, 31 discrete time) use RR ~ HR (rare-event
# approximation; reported for completeness).
# Note: for slope terms (tso, mm/yr) d is per year of exposure.
lp_log("37_evalues")
img <- load_derived("panel_img")
sd_bl <- sapply(names(cfg$brain), function(m) sd(img[[m]][img$is_bl], na.rm = TRUE))
ev_fun <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
rd <- function(f) { p <- file.path(cfg$dir_output, f); if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
rows <- list()
a <- rd("25_onset_anchored.csv")
if (!is.null(a)) rows$a <- a |> filter(model == "main", term %in% c("eventual", "post", "tso")) |> mutate(source = "25_onset_anchored")
w <- rd("32_within_family.csv")
if (!is.null(w)) rows$w <- w |> filter(term %in% c("post_w", "tso_w", "eventual_w")) |> mutate(source = paste0("32_", model))
t <- rd("30_target_trial.csv")
if (!is.null(t)) rows$t <- t |> mutate(outcome = paste0(cfg$primary, "_change_per_yr")) |> mutate(source = paste0("30_", model))
cont <- bind_rows(rows)
out <- list()
if (nrow(cont)) {
  cont$sd <- ifelse(grepl("_change_per_yr$", cont$outcome), NA, sd_bl[sub("_change_per_yr$", "", cont$outcome)])
  if (!is.null(t)) {   # SD of the annualized change among trial participants (aggregate)
    tt <- tryCatch(load_derived("target_trial_data"), error = function(e) NULL)
    if (!is.null(tt)) cont$sd[grepl("_change_per_yr$", cont$outcome)] <- sd(tt$Y, na.rm = TRUE)
  }
  out$cont <- cont |> mutate(
    d = estimate / sd,
    near = ifelse(conf.low > 0, conf.low, ifelse(conf.high < 0, conf.high, 0)) / sd,
    RR = exp(0.91 * abs(d)), E_value = ev_fun(RR),
    E_value_CI = ifelse(near == 0, 1, ev_fun(exp(0.91 * abs(near))))) |>
    transmute(module = "37", source, outcome, term, estimate, conf.low, conf.high, d, RR, E_value, E_value_CI)
}
j <- rd("24_joint_models.csv"); dt <- rd("31_discrete_time_onset.csv")
hz <- bind_rows(
  if (!is.null(j)) j |> filter(grepl("^value\\(", term)) |> transmute(source = paste0("24_", model), outcome = "binge_onset", term, estimate, conf.low, conf.high),
  if (!is.null(dt)) dt |> filter(term == "mz") |> transmute(source = "31_discrete_time", outcome = "binge_onset", term, estimate, conf.low, conf.high))
if (nrow(hz)) out$hz <- hz |> mutate(RR = exp(estimate), E_value = ev_fun(RR),
                                     near = ifelse(conf.low > 0, conf.low, ifelse(conf.high < 0, conf.high, 0)),
                                     E_value_CI = ifelse(near == 0, 1, ev_fun(exp(near))), module = "37") |>
  select(module, source, outcome, term, estimate, conf.low, conf.high, RR, E_value, E_value_CI)
res <- bind_rows(out)
safe_write_csv(res, "37_evalues")
if (nrow(res)) key_write(res |> filter(outcome == cfg$primary | grepl("binge_onset|change_per_yr", outcome)) |>
                           transmute(model = source, outcome, term = paste(term, "E-value"), estimate = E_value, conf.low = E_value_CI), "37")
