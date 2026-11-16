# =============================================================================
# 23_mixture.R : latent class growth mixture models (lcmm::hlme)
# =============================================================================
# Fits 1..K classes for the primary marker (random intercept + slope within
# class, class-specific mean trajectories), selects by BIC, reports class sizes
# (suppressed < 10), class mean trajectories, and binge-onset hazard by class.
# Entropy is reported because poorly separated classes should not be treated
# as real subgroups. Starting values for k > 1 come from the 1-class model
# via gridsearch with random starts.
lp_log("23_mixture")
if (!need(c("lcmm", "survival"))) stop("lcmm and survival required")
hlme <- lcmm::hlme  # lcmm::gridsearch needs an unqualified hlme() call it can find from here
d <- img_data(); m <- cfg$primary
dd <- d |> filter(!is.na(.data[[m]])) |> mutate(y = .data[[m]], sid = as.integer(factor(id))) |>
  select(id, sid, fam, y, time, female, age_bl_c) |> as.data.frame()
fits <- list()
fits[["1"]] <- try_fit("hlme k=1", lcmm::hlme(y ~ time, random = ~ time, subject = "sid", ng = 1, data = dd, verbose = FALSE))
for (k in setdiff(cfg$mixture_k, 1)) {
  if (is.null(fits[["1"]])) break
  fits[[as.character(k)]] <- try_fit(paste0("hlme k=", k), {
    m1 <- fits[["1"]]
    if (isTRUE(cfg$fast)) lcmm::hlme(y ~ time, mixture = ~ time, random = ~ time, subject = "sid", ng = k, data = dd, B = m1, verbose = FALSE)
    else eval(bquote(lcmm::gridsearch(rep = 30, maxiter = 30, minit = m1,
           hlme(y ~ time, mixture = ~ time, random = ~ time, subject = "sid", ng = .(k), data = dd, verbose = FALSE))))  # unqualified: gridsearch cannot call lcmm::hlme
  })
}
fits <- fits[!vapply(fits, is.null, logical(1))]
entropy <- function(f) { p <- as.matrix(f$pprob[, grep("^prob", names(f$pprob)), drop = FALSE]); if (ncol(p) == 1) return(NA)
  1 - sum(-p * log(pmax(p, 1e-12))) / (nrow(p) * log(ncol(p))) }
sel <- bind_rows(lapply(names(fits), function(k) {
  f <- fits[[k]]; sz <- table(factor(f$pprob$class, levels = seq_len(as.integer(k))))
  tibble(module = "23", model = "hlme", outcome = m, k = as.integer(k), loglik = f$loglik, npm = length(f$best),
         BIC = f$BIC, entropy = entropy(f), conv = f$conv, smallest_class_n = min(sz),
         smallest_class_pct = round(100 * min(sz) / sum(sz), 1),
         degenerate = as.integer(k) > 1 && abs(f$loglik - fits[["1"]]$loglik) < 1e-3)
}))
# keep every k in the table; blank only the size cells that fall below min_cell
sel_out <- sel |> mutate(smallest_class_pct = ifelse(smallest_class_n < cfg$min_cell, NA, smallest_class_pct),
                         smallest_class_n = ifelse(smallest_class_n < cfg$min_cell, NA, smallest_class_n))
safe_write_csv(sel_out, "23_mixture_selection")
if (nrow(sel)) {
  # selectable: converged, not a copy of the 1-class fit, and every class holds
  # at least mixture_min_prop of the sample (and at least min_cell people)
  okk <- sel |> filter(conv == 1, !degenerate, smallest_class_n >= cfg$min_cell,
                       k == 1 | smallest_class_pct >= 100 * cfg$mixture_min_prop)
  best <- if (nrow(okk)) okk$k[which.min(okk$BIC)] else 1L
  f <- fits[[as.character(best)]]
  cls <- f$pprob |> select(sid, class) |> left_join(dd |> distinct(sid, id, fam), by = "sid")
  save_derived(cls, "mixture_classes")
  nd <- data.frame(time = seq(0, 6, by = 0.5))
  py <- lcmm::predictY(f, newdata = nd, var.time = "time")$pred
  traj <- tibble(module = "23", model = paste0("hlme_k", best), time = nd$time, as_tibble(py)) |>
    pivot_longer(-c(module, model, time), names_to = "class", values_to = "mean_trajectory")
  sizes <- cls |> count(class, name = "n")
  safe_write_csv(traj |> mutate(class = sub("Ypred_class", "", class)) |>
                   left_join(sizes |> mutate(class = as.character(class)), by = "class"), "23_mixture_trajectories", count_cols = "n")
  lp_log("  selected k = ", best, if (best == 1) " (no well-separated multi-class solution)" else "")
  ev_ok <- FALSE
  if (best > 1) {
    ev <- load_derived("events") |> inner_join(cls |> select(id, class), by = "id") |> mutate(class = factor(class))
    ev_ok <- all(table(ev$class[ev$event == 1]) >= cfg$min_cell) && length(unique(ev$class[ev$event == 1])) == nlevels(ev$class)
    if (!ev_ok) lp_log("  Cox by class skipped: a class has < ", cfg$min_cell, " onsets")
  }
  # drop a by-class Cox table left from an earlier run, so it is never read as current
  if (!(best > 1 && ev_ok)) unlink(file.path(cfg$dir_output, "23_onset_by_class.csv"))
  if (best > 1 && ev_ok) {
    cf <- try_fit("Cox by class", survival::coxph(survival::Surv(t_event, event) ~ class + female + age_bl, data = ev, cluster = fam))
    if (!is.null(cf)) {
      tb <- broom::tidy(cf, conf.int = TRUE, exponentiate = TRUE) |>
        transmute(module = "23", model = "cox_onset_by_class", term, HR = estimate, conf.low, conf.high, p.value,
                  n_id = nrow(ev), n_events = sum(ev$event))
      safe_write_csv(tb, "23_onset_by_class", count_cols = c("n_id", "n_events"))
    }
  }
}
