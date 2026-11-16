# =============================================================================
# 21_fpca.R : functional PCA of sparse brain trajectories (fdapace / PACE)
# and multivariate FPCA (Happ & Greven 2018), with FPC scores as predictors of
# binge onset.
# =============================================================================
# Each participant contributes 1-4 scans at irregular ages; PACE estimates the
# mean and covariance surfaces by pooling across participants and predicts
# individual FPC scores by conditional expectation. Leading FPCs summarise
# level (FPC1) and rate/shape of change (FPC2+).
#
# Leakage guard: scores used to predict onset are re-predicted from scans
# acquired at or before the last binge-negative assessment (t_left), so a
# post-onset change cannot "predict" the onset that caused it.
lp_log("21_fpca")
if (!need(c("fdapace", "survival"))) stop("fdapace and survival required")
d <- img_data(); ev <- load_derived("events")
fp <- list(); fp_ids <- list(); scores_pre <- list(); tabs <- list(); fve <- list()
to_lists <- function(dm, var) {
  sp <- split(dm, dm$id)
  sp <- sp[vapply(sp, nrow, 1L) >= 1]
  list(id = names(sp), Ly = lapply(sp, function(z) z[[var]]), Lt = lapply(sp, function(z) z$age))
}
for (m in names(cfg$brain)) {
  v <- paste0(m, "_z")
  # ages rounded to 0.1 y: fdapace builds its covariance on the grid of unique
  # observed times, and continuous visit ages make that grid ~N_obs wide.
  dm <- d |> filter(!is.na(.data[[v]])) |> mutate(age = round(age, 1)) |> arrange(id, age) |> distinct(id, age, .keep_all = TRUE)
  L <- to_lists(dm, v)
  op <- list(dataType = "Sparse", methodSelectK = "FVE", FVEthreshold = cfg$fpca_fve, nRegGrid = 41, verbose = FALSE)
  fit <- try_fit(paste("FPCA", m), fdapace::FPCA(L$Ly, L$Lt, op))
  if (is.null(fit)) next
  if (fit$selectK < 2) {   # keep at least level + change components
    op$methodSelectK <- 2L
    fit2 <- try_fit(paste("FPCA", m, "K=2"), fdapace::FPCA(L$Ly, L$Lt, op)); if (!is.null(fit2)) fit <- fit2
  }
  K <- min(3, fit$selectK)
  fp[[m]] <- fit; fp_ids[[m]] <- L$id
  tabs[[m]] <- tibble(module = "21", marker = m, age = fit$workGrid, mean = fit$mu,
                      !!!setNames(lapply(seq_len(K), function(k) fit$phi[, k]), paste0("phi", seq_len(K))))
  fve[[m]] <- tibble(module = "21", marker = m, component = seq_along(fit$lambda), eigenvalue = fit$lambda,
                     cum_fve = fit$cumFVE[seq_along(fit$lambda)], n_id = length(L$id))
  # pre-onset scores
  pre <- dm |> left_join(ev |> select(id, t_left), by = "id") |> filter(id %in% ev$id, is.na(t_left) | time <= t_left + 1e-8)
  Lp <- to_lists(pre, v)
  sc <- tryCatch(predict(fit, Lp$Ly, Lp$Lt, K = K)$scores, error = function(e) {
    lp_log("  predict.FPCA failed (", conditionMessage(e), "); using full-data scores"); fit$xiEst[match(Lp$id, L$id), 1:K, drop = FALSE] })
  scores_pre[[m]] <- tibble(id = Lp$id, !!!setNames(lapply(seq_len(K), function(k) as.numeric(scale(sc[, k]))), paste0(m, "_fpc", seq_len(K))))
}
save_derived(fp, "fpca_fits"); save_derived(scores_pre, "fpca_scores_preonset")
safe_write_csv(bind_rows(tabs), "21_fpca_mean_eigenfunctions")
safe_write_csv(bind_rows(fve), "21_fpca_variance_explained", count_cols = "n_id")

## ---- Multivariate FPCA (Happ & Greven): eigen-decomposition of the joint
## covariance of univariate scores across markers ---------------------------
mk <- intersect(names(cfg$brain), names(fp))
if (length(mk) >= 2) {
  S <- Reduce(function(a, b) inner_join(a, b, by = "id"), lapply(mk, function(m) {
    f <- fp[[m]]; K <- min(3, f$selectK)
    tibble(id = fp_ids[[m]], !!!setNames(lapply(seq_len(K), function(k) f$xiEst[, k]), paste0(m, "_", seq_len(K))))
  }))
  Z <- as.matrix(S[, -1]); e <- eigen(cov(Z), symmetric = TRUE)
  nm <- min(4, ncol(Z))
  mf_var <- tibble(module = "21", component = paste0("MFPC", seq_len(nm)), eigenvalue = e$values[seq_len(nm)],
                   prop_var = e$values[seq_len(nm)] / sum(e$values), n_id = nrow(Z))
  # multivariate eigenfunctions per marker: psi_m(t) = Phi_m(t) %*% c_m
  ef <- bind_rows(lapply(seq_len(nm), function(j) bind_rows(lapply(mk, function(m) {
    f <- fp[[m]]; K <- min(3, f$selectK); cols <- grep(paste0("^", m, "_"), colnames(Z))
    tibble(module = "21", component = paste0("MFPC", j), marker = m, age = f$workGrid,
           psi = as.numeric(f$phi[, 1:K, drop = FALSE] %*% e$vectors[cols, j]))
  }))))
  safe_write_csv(mf_var, "21_mfpca_variance_explained", count_cols = "n_id")
  safe_write_csv(ef, "21_mfpca_eigenfunctions")
}

## ---- FPC scores (pre-onset) as predictors of binge onset ------------------------
cx <- list()
for (m in names(scores_pre)) {
  dd <- ev |> inner_join(scores_pre[[m]], by = "id")
  f <- as.formula(paste0("survival::Surv(t_event, event) ~ ", paste(grep("_fpc", names(dd), value = TRUE), collapse = " + "),
                         " + female + age_bl"))
  cf <- try_fit(paste("Cox on FPC", m), survival::coxph(f, data = dd, cluster = fam))
  if (!is.null(cf)) cx[[m]] <- broom::tidy(cf, conf.int = TRUE, exponentiate = TRUE) |>
    transmute(module = "21", model = "cox_preonset_fpc", marker = m, term, HR = estimate, conf.low, conf.high,
              p.value, n_id = nrow(dd), n_events = sum(dd$event))
}
cx <- bind_rows(cx)
safe_write_csv(cx, "21_cox_on_fpc_scores", count_cols = c("n_id", "n_events"))
if (nrow(cx)) key_write(cx |> filter(grepl("fpc", term)) |> select(model, marker, term, estimate = HR, conf.low, conf.high, p.value, n_id, n_events), "21")
