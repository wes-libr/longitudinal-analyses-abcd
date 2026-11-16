# =============================================================================
# 13_sem_growth.R : latent growth, latent change score and RI-CLPM (lavaan)
# =============================================================================
# All models: FIML for missing waves, cluster-robust SEs by family
# (cluster = "fam", MLR). Imaging waves are two years apart; loadings are in
# years (0, 2, 4, 6) so slopes are per year, matching the LMMs.
#  LGCM        linear growth in the primary marker
#  LGCM_TIC    + sex, baseline age, caregiver education, eventual binge onset
#  LGCM_group  multigroup by sex (equal vs free slope means, LRT)
#  LCSM        dual-change latent change score model
#  BLCSM       bivariate LCSM: marker <-> CBCL externalizing (coupling)
#  RICLPM      random-intercept cross-lagged panel: marker <-> externalizing
lp_log("13_sem_growth")
if (!need("lavaan")) stop("lavaan required")
d <- img_data(); m <- cfg$primary
ev <- d |> distinct(id, eventual)
w <- img_wide(d |> mutate(y = .data[[m]], x = ext), c("y", "x")) |> left_join(ev, by = "id")
# rescale so variances are comparable (lavaan convergence)
sy <- sd(w$y_0, na.rm = TRUE); sx <- sd(w$x_0, na.rm = TRUE)
for (k in c(0, 2, 4, 6)) { w[[paste0("y_", k)]] <- w[[paste0("y_", k)]] / sy; w[[paste0("x_", k)]] <- w[[paste0("x_", k)]] / sx }
fitsem <- function(model, data = w, ...) lavaan::sem(model, data = data, missing = "fiml", estimator = "MLR", cluster = "fam", ...)
gr <- function(model, data = w, ...) lavaan::growth(model, data = data, missing = "fiml", estimator = "MLR", cluster = "fam", ...)

lgcm <- "i =~ 1*y_0 + 1*y_2 + 1*y_4 + 1*y_6
         s =~ 0*y_0 + 2*y_2 + 4*y_4 + 6*y_6"
lgcm_tic <- paste(lgcm, "\n i + s ~ female + age_bl + edu_bl + eventual")
lcsm <- "
  ly0 =~ 1*y_0; ly2 =~ 1*y_2; ly4 =~ 1*y_4; ly6 =~ 1*y_6
  ly2 ~ 1*ly0; ly4 ~ 1*ly2; ly6 ~ 1*ly4
  d2 =~ 1*ly2; d4 =~ 1*ly4; d6 =~ 1*ly6
  y_0 ~~ ve*y_0; y_2 ~~ ve*y_2; y_4 ~~ ve*y_4; y_6 ~~ ve*y_6
  y_0 + y_2 + y_4 + y_6 ~ 0*1
  ly0 ~ 1; ly0 ~~ ly0
  ys =~ 1*d2 + 1*d4 + 1*d6
  ys ~ 1; ys ~~ ys; ys ~~ ly0
  d2 ~ beta*ly0; d4 ~ beta*ly2; d6 ~ beta*ly4
  d2 ~ 0*1; d4 ~ 0*1; d6 ~ 0*1; d2 ~~ 0*d2; d4 ~~ 0*d4; d6 ~~ 0*d6
  ly2 ~~ 0*ly2; ly4 ~~ 0*ly4; ly6 ~~ 0*ly6"
blcsm <- paste(lcsm, "
  lx0 =~ 1*x_0; lx2 =~ 1*x_2; lx4 =~ 1*x_4; lx6 =~ 1*x_6
  lx2 ~ 1*lx0; lx4 ~ 1*lx2; lx6 ~ 1*lx4
  dx2 =~ 1*lx2; dx4 =~ 1*lx4; dx6 =~ 1*lx6
  x_0 ~~ vx*x_0; x_2 ~~ vx*x_2; x_4 ~~ vx*x_4; x_6 ~~ vx*x_6
  x_0 + x_2 + x_4 + x_6 ~ 0*1
  lx0 ~ 1; lx0 ~~ lx0
  xs =~ 1*dx2 + 1*dx4 + 1*dx6
  xs ~ 1; xs ~~ xs; xs ~~ lx0
  dx2 ~ bx*lx0; dx4 ~ bx*lx2; dx6 ~ bx*lx4
  dx2 ~ 0*1; dx4 ~ 0*1; dx6 ~ 0*1; dx2 ~~ 0*dx2; dx4 ~~ 0*dx4; dx6 ~~ 0*dx6
  lx2 ~~ 0*lx2; lx4 ~~ 0*lx4; lx6 ~~ 0*lx6
  d2 ~ g_xy*lx0; d4 ~ g_xy*lx2; d6 ~ g_xy*lx4
  dx2 ~ g_yx*ly0; dx4 ~ g_yx*ly2; dx6 ~ g_yx*ly4
  ly0 + ys ~~ lx0 + xs")
riclpm <- "
  RIy =~ 1*y_0 + 1*y_2 + 1*y_4 + 1*y_6
  RIx =~ 1*x_0 + 1*x_2 + 1*x_4 + 1*x_6
  wy0 =~ 1*y_0; wy2 =~ 1*y_2; wy4 =~ 1*y_4; wy6 =~ 1*y_6
  wx0 =~ 1*x_0; wx2 =~ 1*x_2; wx4 =~ 1*x_4; wx6 =~ 1*x_6
  y_0 ~~ 0*y_0; y_2 ~~ 0*y_2; y_4 ~~ 0*y_4; y_6 ~~ 0*y_6
  x_0 ~~ 0*x_0; x_2 ~~ 0*x_2; x_4 ~~ 0*x_4; x_6 ~~ 0*x_6
  wy2 ~ ay*wy0 + cxy*wx0; wy4 ~ ay*wy2 + cxy*wx2; wy6 ~ ay*wy4 + cxy*wx4
  wx2 ~ ax*wx0 + cyx*wy0; wx4 ~ ax*wx2 + cyx*wy2; wx6 ~ ax*wx4 + cyx*wy4
  wy0 ~~ wx0; wy2 ~~ wx2; wy4 ~~ wx4; wy6 ~~ wx6
  RIy ~~ RIx; RIy ~~ RIy; RIx ~~ RIx
  RIy + RIx ~~ 0*wy0 + 0*wx0"

fits <- list(
  LGCM     = try_fit("LGCM", gr(lgcm)),
  LGCM_TIC = try_fit("LGCM_TIC", gr(lgcm_tic)),
  LCSM     = try_fit("LCSM", fitsem(lcsm)),
  BLCSM    = try_fit("BLCSM", fitsem(blcsm)),
  RICLPM   = try_fit("RICLPM", fitsem(riclpm)))
grp_free <- try_fit("LGCM by sex (free)", gr(lgcm, group = "female"))
grp_eq   <- try_fit("LGCM by sex (equal slope mean)", gr(lgcm, group = "female", group.equal = "means"))

pe <- bind_rows(lapply(names(fits), function(k) {
  f <- fits[[k]]; if (is.null(f)) return(NULL)
  pe <- lavaan::parameterEstimates(f) |> as_tibble()
  if (!"label" %in% names(pe)) pe$label <- ""
  pe |> filter(op %in% c("~", "~1", "~~"), !(op == "~~" & lhs == rhs & se == 0), !(op == "~1" & se == 0)) |>
    transmute(module = "13", model = k, outcome = m, term = paste(lhs, op, rhs), label, estimate = est,
              std.error = se, conf.low = ci.lower, conf.high = ci.upper, p.value = pvalue)
}))
fitidx <- bind_rows(lapply(names(fits), function(k) {
  f <- fits[[k]]; if (is.null(f)) return(NULL)
  fm <- lavaan::fitMeasures(f, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "aic", "bic"))
  tibble(module = "13", model = k, outcome = m, !!!as.list(round(fm, 4)), n_id = lavaan::lavInspect(f, "nobs"))
}))
if (!is.null(grp_free) && !is.null(grp_eq)) {
  lrt <- lavaan::lavTestLRT(grp_eq, grp_free)
  pe <- bind_rows(pe, tibble(module = "13", model = "LGCM_group_sex", outcome = m,
                             term = "equal growth means across sex (LRT)", estimate = lrt$`Chisq diff`[2], p.value = lrt$`Pr(>Chisq)`[2]))
}
note <- sprintf("Marker and externalizing scaled by baseline SD (marker SD = %.4g, ext SD = %.4g).", sy, sx)
safe_write_csv(pe |> mutate(note = note), "13_sem_parameters")
safe_write_csv(fitidx, "13_sem_fit", count_cols = "n_id")
