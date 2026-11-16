# =============================================================================
# config.R : ABCD 7.0 longitudinal pipeline (brain trajectories x binge onset)
# =============================================================================
# Edit this file only. Every module reads `cfg`.
#
# Blinding model (same as the music-ed pipeline): participant-level data and
# model objects go to cfg$dir_derived, which must be OUTSIDE the shared repo.
# The repo receives aggregate tables only, written through safe_write_csv()
# with no cell below 10.

cfg <- list(

  ## ---- Paths ----------------------------------------------------------------
  # NBDC tabulated phenotype directory (one .parquet per table).
  dir_data = Sys.getenv(
    "ABCD_DATA_PATH",
    "/Users/wesleythompson/Dropbox/aaa_mac/research/brains/nbdc-datasets/abcd_7-0_phenotype/7_0/tabulated/release_abcd"),
  # Host-only cache for panels and fitted models (participant level).
  dir_derived = Sys.getenv(
    "LP_DERIVED_PATH",
    "/Users/wesleythompson/Dropbox/aaa_mac/research/brains/nbdc-datasets/abcd_7-0_phenotype/7_0/derived/longitudinal_pipeline"),
  # Aggregate outputs (inside the repo).
  dir_output = "output",

  ## ---- Release ------------------------------------------------------------
  # Used for get_dd_abcd()/get_sessions_abcd() validation. create_dataset_abcd()
  # is called WITHOUT `release=` because NBDCtools 1.1.0 aborts in the
  # factor-transform step when release is explicit; the default ("latest")
  # resolves to 7.0 with NBDCtoolsData >= 5.0.0 (see music-ed 01_data.R).
  release = "7.0",

  ## ---- Brain markers (biennial imaging: ses-00A/02A/04A/06A) ---------------
  # name = analysis label, value = 7.0 variable (or a vector summed L+R).
  brain = list(
    thk     = "mr_y_smri__thk__dsk_mean",          # mean cortical thickness (mm)
    area    = "mr_y_smri__area__dsk_sum",          # total cortical surface area
    cortvol = "mr_y_smri__vol__dsk_sum",           # total cortical volume
    scgv    = "mr_y_smri__vol__aseg__scgv_sum",    # subcortical gray volume
    hippo   = c("mr_y_smri__vol__aseg__hc__lh_sum", "mr_y_smri__vol__aseg__hc__rh_sum"),
    amyg    = c("mr_y_smri__vol__aseg__ag__lh_sum", "mr_y_smri__vol__aseg__ag__rh_sum")
  ),
  icv          = "mr_y_smri__vol__aseg__icv_sum",
  qc_t1        = "mr_y_qc__incl__smri__t1_indicator",  # keep == 1
  scanner      = "mr_y_adm__info__dev_serial",
  # Markers adjusted for ICV in the models (volumes/area; not thickness).
  icv_adjust   = c("area", "cortvol", "scgv", "hippo", "amyg"),
  # Primary marker(s) for the joint model and the single-outcome modules.
  primary      = "thk",
  joint_multi  = c("thk", "scgv"),

  ## ---- Binge drinking (TLFB) ---------------------------------------------
  # Onset = first assessment (annual or mid-year) with >= 1 binge day.
  binge_vars   = c("su_y_tlfb__alc__binge_ud",
                   "su_y_tlfb__alc__binge__1mo_ud",
                   "su_y_tlfb__alc__binge__3mo_ud"),
  alc_lifetime = "su_y_dyn__alc__lt_indicator",
  # Event-time convention for the interval (last negative, first positive]:
  #   "right"    = age at the first positive assessment (default)
  #   "midpoint" = midpoint of the interval (sensitivity)
  onset_time   = "right",

  ## ---- Covariates ----------------------------------------------------------
  site   = "ab_g_dyn__design_site",
  family = "ab_g_stc__design_id__fam",
  age    = "ab_g_dyn__visit_age",
  sex    = "ab_g_stc__cohort_sex",
  edu    = "ab_g_dyn__cohort_edu__cgs",
  income = "ab_g_dyn__cohort_income__hhold__3lvl",
  weight = "ab_g_stc__cohort_acsps",                 # raked propensity weight
  pds    = c(f = "ph_y_pds__f_categ", m = "ph_y_pds__m_categ"),
  ext    = "mh_p_cbcl__synd__ext_tscore",            # second process (RI-CLPM, LCSM)

  ## ---- Temporal trends (age-period-cohort) ---------------------------------
  visit_dt     = "ab_g_dyn__visit_dtt",               # visit date-time -> calendar period
  visit_type   = "ab_g_dyn__visit_type",              # 1 on-site, 2 remote, 3 hybrid
  covid_start  = "2020-03-15",                        # start of pandemic period
  covid_end    = "2021-06-30",                        # end of the restricted (remote-heavy) era

  ## ---- Regional, family and puberty fields (modules 32, 34, 36) ------------
  roi_table    = "mr_y_smri__thk__dsk",               # Desikan cortical thickness table
  roi_pattern  = "^mr_y_smri__thk__dsk__.+__(lh|rh)_mean$",
  birth_id     = "ab_g_stc__design_id__birth",        # shared by twins/triplets
  pds_mean     = c(f = "ph_y_pds__f_mean", m = "ph_y_pds__m_mean"),

  ## ---- Modules -------------------------------------------------------------
  modules = c(
    "02_describe_missing", "10_change", "11_lmm", "12_glmm_gee",
    "13_sem_growth", "20_gamm_age", "21_fpca", "22_multivariate",
    "23_mixture", "24_joint_model", "25_onset_anchored", "26_missing_data",
    "28_temporal_trends", "30_target_trial", "31_discrete_time_onset",
    "32_within_family", "33_harmonization", "34_regional", "35_dynamic_prediction",
    "36_puberty_time", "37_evalues", "27_synthesis"),

  ## ---- Module settings -----------------------------------------------------
  min_cell        = 10L,
  seed            = 20261114L,
  scanner_re      = TRUE,     # crossed random intercept for scanner serial
  mixture_k       = 1:4,      # latent classes tried in 23_mixture
  mixture_min_prop = 0.05,    # smallest class share for a mixture solution to be selectable (23)
  fpca_fve        = 0.95,     # FPCA fraction of variance explained
  jm_iter         = 6500L,    # JMbayes2 MCMC iterations per chain
  jm_burnin       = 1500L,
  jm24_iter       = 20000L,   # longer chains for the flagship joint models (24) so Rhat < 1.05
  jm24_burnin     = 5000L,
  jm24_thin       = 5L,       # keeps 3,000 draws per chain and the saved fits small
  jm_cores        = 1L,       # chains in sequence: parallel chains (3) ran 16x slower on the Mac (memory), 2026-11-15
  mi_m            = 20L,      # imputations in 26_missing_data
  delta_grid      = c(0, -0.1, -0.2, -0.3, -0.5),  # MNAR shifts, SD units
  tt_boot         = 200L,     # family-cluster bootstrap replicates (30_target_trial g-formula)
  cv_folds        = 5L,       # family-level folds for 35_dynamic_prediction
  dp_tstart       = 4,        # landmark (years since baseline) for dynamic prediction; onsets cluster at ages 14-16, so year 4 (age ~14) gives a window with enough events
  dp_dt           = 2,        # prediction horizon (years) after the landmark
  auto_install    = TRUE,     # install missing packages before the run (CRAN; longCombat from GitHub)
  stop_if_missing = TRUE,     # stop before any module if a needed package is still missing
  fast            = FALSE     # TRUE = fewer iterations (testing only)
)
