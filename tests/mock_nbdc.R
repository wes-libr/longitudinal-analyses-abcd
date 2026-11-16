# =============================================================================
# tests/mock_nbdc.R : simulated ABCD 7.0 data in NBDCtools' create_dataset_abcd()
# output shape, plus stand-ins for the NBDCtools metadata functions.
# The real filter_events_abcd()/filter_empty_rows() are sourced from the
# NBDCtools package source when available (tests/NBDCtools_filter.R).
# =============================================================================
# Known truth built into the simulation (checked by tests/run_mock.R):
#   * thickness declines ~0.015 mm/yr; site, family, participant, scanner effects
#   * binge-onset hazard rises with lower current thickness: alpha = -0.4 per 0.1 mm
#   * after binge onset, thickness declines an extra 0.010 mm/yr
#   * imaging attrition depends on externalizing (MAR)
#   * regional: 6 Desikan ROIs x 2 hemispheres; lateral orbitofrontal thins an
#     extra 0.020 mm/yr after onset (the other ROIs follow global thickness)
#   * puberty: earlier maturers (timing +1 SD) have log HR +0.3 for onset
#   * period: thickness declines an extra 0.008 mm per year lived since
#     2020-03-15; binge-onset hazard is halved (log HR -0.7) in the restricted
#     pandemic era (2020-03-15 to 2021-06-30); remote visits in that era

simulate_abcd <- function(n_fam = 1300, seed = 1) {
  set.seed(seed)
  fs <- sample(1:3, n_fam, replace = TRUE, prob = c(.80, .17, .03))
  fam <- tibble(fam = sprintf("FAM%05d", seq_len(n_fam)), size = fs,
                site = sprintf("site%02d", sample(1:22, n_fam, TRUE)), u_fam = rnorm(n_fam, 0, .04))
  site_eff <- setNames(rnorm(22, 0, .02), sprintf("site%02d", 1:22))
  kids <- fam[rep(seq_len(n_fam), fam$size), ] |>
    mutate(participant_id = sprintf("sub-%08d", row_number()),
           sex = sample(1:2, n(), TRUE), age_bl = runif(n(), 9, 11), cal_bl = runif(n(), 2016.7, 2018.8),
           edu = sample(1:5, n(), TRUE, prob = c(.05, .10, .25, .30, .30)),   # 7.0 coding: 5 levels
           income = sample(1:3, n(), TRUE), wt = rlnorm(n(), 0, .4),
           b0 = rnorm(n(), 0, .10), b1 = rnorm(n(), 0, .006),
           e0 = rnorm(n(), 0, 1), ptime = rnorm(n(), 0, 1), z_area = .4 * b0 / .1 + rnorm(n(), 0, .9),
           z_scg = .4 * b0 / .1 + rnorm(n(), 0, .9), icv = rnorm(n(), 1.5e6, 1.2e5))
  n <- nrow(kids)
  c0 <- 2020 + 74 / 365.25; c1 <- 2021 + 180 / 365.25
  pand <- function(i, age) pmax(0, kids$cal_bl[i] + age - kids$age_bl[i] - c0)
  thk_true <- function(i, age, post_years = 0)
    -0.008 * pand(i, age) + 2.75 + site_eff[kids$site[i]] + kids$u_fam[i] + kids$b0[i] + 0.02 * (kids$sex[i] == 2) +
      (-0.015 + kids$b1[i]) * (age - 10) - 0.010 * post_years
  # binge onset: hazard on age scale with current true thickness
  onset <- rep(Inf, n)
  for (i in seq_len(n)) {
    for (a in seq(kids$age_bl[i] + 0.05, kids$age_bl[i] + 8, by = 0.05)) {
      calp <- kids$cal_bl[i] + a - kids$age_bl[i]
      h <- exp(-5.6 + 0.55 * (a - 12) - 4 * (thk_true(i, a) - 2.75) - 0.7 * (calp >= c0 & calp < c1) + 0.3 * kids$ptime[i])
      if (runif(1) < 1 - exp(-h * 0.05)) { onset[i] <- a; break }
    }
  }
  kids$onset <- onset
  ses <- c(sprintf("ses-%02dA", 0:7), sprintf("ses-%02dM", 0:6), "ses-00S")
  rows <- tidyr::expand_grid(i = seq_len(n), session_id = ses) |>
    mutate(k = as.numeric(substr(session_id, 5, 6)), typ = substr(session_id, 7, 7)) |>
    filter(!(typ == "M" & k < 1)) |>
    mutate(age = kids$age_bl[i] + k + ifelse(typ == "M", .5, 0) + rnorm(n(), 0, .12),
           age = ifelse(typ == "S", kids$age_bl[i] - .1, age))
  # attendance falls over time (lower for high externalizing)
  ext <- 50 + 10 * (kids$e0[rows$i] * .8 + rnorm(nrow(rows), 0, .6))
  att <- runif(nrow(rows)) < plogis(3 - 0.25 * rows$k - 0.3 * (ext - 50) / 10)
  rows <- rows[att | rows$session_id %in% c("ses-00A", "ses-00S"), ]
  ext <- ext[att | rows$session_id %in% c("ses-00A", "ses-00S")]
  is_img <- rows$typ == "A" & rows$k %% 2 == 0 & rows$k <= 6
  post <- pmax(0, rows$age - kids$onset[rows$i]); post[!is.finite(post)] <- 0
  thk <- thk_true(rows$i, rows$age, post)
  scanner <- paste0(kids$site[rows$i], ifelse(rows$k >= 4 & as.integer(substr(kids$site[rows$i], 5, 6)) %% 2 == 0, "_B", "_A"))
  scan_eff <- setNames(rnorm(44, 0, .02), sort(unique(c(paste0(sprintf("site%02d", 1:22), "_A"), paste0(sprintf("site%02d", 1:22), "_B")))))
  thk_obs <- thk + scan_eff[scanner] + rnorm(nrow(rows), 0, .03)
  qc <- ifelse(is_img, as.integer(runif(nrow(rows)) < .9), NA)
  # biennial imaging columns
  img_na <- function(x) ifelse(is_img, x, NA)
  agec <- rows$age - 10
  icv <- kids$icv[rows$i]
  pmean <- pmin(4, pmax(1, 1 + 0.45 * (rows$age - 9.5 + 0.8 * kids$ptime[rows$i]) + 0.25 * (kids$sex[rows$i] == 2) + rnorm(nrow(rows), 0, .2)))
  rois <- c("cmfrt", "lobfrt", "prcn", "sm", "ins", "ff")
  out <- tibble(
    participant_id = kids$participant_id[rows$i],
    session_id = rows$session_id,
    ab_g_dyn__design_site = ifelse(rows$typ == "S", NA, kids$site[rows$i]),
    ab_g_dyn__visit_age = ifelse(rows$typ == "S", NA, rows$age),
    ab_g_dyn__cohort_edu__cgs = factor(ifelse(rows$typ == "A", kids$edu[rows$i], NA)),
    ab_g_dyn__cohort_income__hhold__3lvl = factor(ifelse(rows$typ == "A", kids$income[rows$i], NA)),
    ab_g_dyn__visit_dtt = ifelse(rows$typ == "S", NA, format(as.Date("2016-01-01") + round((kids$cal_bl[rows$i] + rows$age - kids$age_bl[rows$i] - 2016) * 365.25), "%Y-%m-%d 10:00:00")),
    ab_g_dyn__visit_type = factor(ifelse(rows$typ == "S", NA, ifelse(!is_img & (kids$cal_bl[rows$i] + rows$age - kids$age_bl[rows$i]) >= c0 &
                                    (kids$cal_bl[rows$i] + rows$age - kids$age_bl[rows$i]) < c1 & runif(nrow(rows)) < .7, 2, 1))),
    mr_y_smri__thk__dsk_mean = img_na(thk_obs),
    mr_y_smri__area__dsk_sum = img_na(1.9e5 + 0.08 * (icv - 1.5e6) + 6000 * kids$z_area[rows$i] - 900 * agec + rnorm(nrow(rows), 0, 1500)),
    mr_y_smri__vol__dsk_sum = img_na(5.2e5 + 0.25 * (icv - 1.5e6) + 20000 * (thk - 2.75) / .1 - 4000 * agec + rnorm(nrow(rows), 0, 5000)),
    mr_y_smri__vol__aseg__scgv_sum = img_na(6.0e4 + 0.02 * (icv - 1.5e6) + 2500 * kids$z_scg[rows$i] - 300 * agec + rnorm(nrow(rows), 0, 600)),
    mr_y_smri__vol__aseg__hc__lh_sum = img_na(4200 + 300 * kids$z_scg[rows$i] + 20 * agec + rnorm(nrow(rows), 0, 80)),
    mr_y_smri__vol__aseg__hc__rh_sum = img_na(4300 + 300 * kids$z_scg[rows$i] + 20 * agec + rnorm(nrow(rows), 0, 80)),
    mr_y_smri__vol__aseg__ag__lh_sum = img_na(1500 + 120 * kids$z_scg[rows$i] + 15 * agec + rnorm(nrow(rows), 0, 50)),
    mr_y_smri__vol__aseg__ag__rh_sum = img_na(1550 + 120 * kids$z_scg[rows$i] + 15 * agec + rnorm(nrow(rows), 0, 50)),
    mr_y_smri__vol__aseg__icv_sum = img_na(icv + rnorm(nrow(rows), 0, 5000)),
    mr_y_qc__incl__smri__t1_indicator = factor(qc),
    mr_y_adm__info__dev_serial = ifelse(is_img, scanner, NA),
    su_y_tlfb__alc__binge_ud = ifelse(rows$typ %in% c("A", "M"),
                                      ifelse(rows$age > kids$onset[rows$i], rbinom(nrow(rows), 1, .75) * rpois(nrow(rows), 2) + (rows$age > kids$onset[rows$i]) * 0, 0), NA),
    su_y_dyn__alc__lt_indicator = factor(ifelse(rows$typ == "A", as.integer(rows$age > kids$onset[rows$i] - 1.5), NA)),
    ph_y_pds__f_categ = factor(ifelse(rows$typ == "A" & kids$sex[rows$i] == 2, pmin(5, pmax(1, round(1 + (rows$age - 9) * .6 + rnorm(nrow(rows), 0, .5)))), NA)),
    ph_y_pds__m_categ = factor(ifelse(rows$typ == "A" & kids$sex[rows$i] == 1, pmin(5, pmax(1, round((rows$age - 9) * .55 + rnorm(nrow(rows), 0, .5)))), NA)),
    mh_p_cbcl__synd__ext_tscore = ifelse(rows$typ == "A", round(ext), NA),
    ph_y_pds__f_mean = ifelse(rows$typ == "A" & kids$sex[rows$i] == 2, round(pmean * 4) / 4, NA),
    ph_y_pds__m_mean = ifelse(rows$typ == "A" & kids$sex[rows$i] == 1, round(pmean * 4) / 4, NA)
  )
  for (r in rois) for (h in c("lh", "rh"))
    out[[paste0("mr_y_smri__thk__dsk__", r, "__", h, "_mean")]] <-
      img_na(thk_obs + 0.1 * (match(r, rois) - 3) - ifelse(r == "lobfrt", 0.020 * post, 0) + rnorm(nrow(rows), 0, .05))
  # first positive TLFB must be >0 days: make post-onset positives at least 1 day
  pos <- !is.na(out$su_y_tlfb__alc__binge_ud) & rows$age > kids$onset[rows$i]
  out$su_y_tlfb__alc__binge_ud[pos] <- pmax(out$su_y_tlfb__alc__binge_ud[pos], rbinom(sum(pos), 1, .75))
  # Static session (session_id = NA) carries the ab_g_stc__ fields, as in 7.0
  stat <- tibble(participant_id = kids$participant_id, session_id = NA_character_,
                 ab_g_stc__design_id__fam = kids$fam, ab_g_stc__cohort_sex = factor(kids$sex),
                 ab_g_stc__cohort_acsps = kids$wt,
                 # twins: in 40% of two-child families the children share a birth id
                 ab_g_stc__design_id__birth = ifelse(kids$size == 2 & as.integer(substr(kids$fam, 4, 8)) %% 5 < 2,
                                                     paste0("B", kids$fam), paste0("B", kids$participant_id)))
  out <- bind_rows(out |> filter(substr(session_id, 7, 7) != "S"), stat) |>
    mutate(session_id = factor(session_id, levels = ses))
  attr(out, "truth") <- kids
  out
}

## ---- NBDCtools stand-ins -----------------------------------------------------
.mock <- new.env()
get_latest_release_abcd <- function() "7.0"
get_id_cols_abcd <- function(...) c("participant_id", "session_id")
get_id_cols_hbcd <- function(...) c("participant_id", "session_id", "run_id")
get_dd_abcd <- function(release = "latest", vars = NULL, tables = NULL) {
  nm <- setdiff(names(.mock$data), c("participant_id", "session_id"))
  dd <- tibble(name = nm, table_name = ifelse(grepl("^mr_y_smri__thk__dsk__", nm), "mr_y_smri__thk__dsk",
                                              sub("^(([a-z]+_){2}[a-z]+).*", "\\1", nm)), label = nm)
  if (!is.null(vars)) dd <- dd |> filter(name %in% vars)
  if (!is.null(tables)) dd <- dd |> filter(table_name %in% tables)
  dd
}
get_sessions_abcd <- function(release = "latest")
  tibble(session_id = levels(.mock$data$session_id), label = levels(.mock$data$session_id))
create_dataset_abcd <- function(dir_data, vars, ...) {
  d <- .mock$data[, c("participant_id", "session_id", intersect(vars, names(.mock$data)))]
  d[rowSums(!is.na(d[, -(1:2)])) > 0, ]
}
