# ABCD 7.0 longitudinal pipeline: brain trajectories and binge-drinking onset

This pipeline builds on the longitudinal.dev tutorials, with the fixes from
`CHANGES.md` applied. It loads data with NBDCtools and runs a chain of
modules around one substantive model: how structural brain trajectories relate
to the onset of binge drinking. Each tutorial method appears as a module, plus
new ones (GAMM, FPCA/MFPCA, multivariate models, RI-CLPM, growth mixtures,
joint models and missing-data sensitivity).

## Run

```r
setwd("path/to/longitudinal-analyses-abcd")   # set dir_data and dir_derived in config.R first
source("run_all.R")          # real data: NBDCtools + NBDCtoolsData >= 5.0.0
source("tests/run_mock.R")   # simulated data, no ABCD files needed
```

The LaTeX manual is `abcd-analyses/manuscripts/longitudinal_manual/main.tex` (compiled `main.pdf` beside it).

Edit `config.R` only: paths, markers, event definition, modules and settings.
`cfg$modules` controls which modules run, and a module that fails is logged
without stopping the run (`output/99_run_status.csv`). Before any module runs,
`lp_check_packages()` installs missing packages (`cfg$auto_install`) and stops
if one is still missing (`cfg$stop_if_missing`), so no module is skipped silently.
To rerun selected modules: `.lp_args <- list(modules = c("11_lmm", "27_synthesis")); source("run_all.R")`.

## Blinding

Same rules as the music-ed pipeline:

- Panels and fitted models go to `cfg$dir_derived`, which must be outside
  `abcd-analyses/`. `lp_setup()` refuses to run otherwise.
- `output/` receives only aggregate tables, written by `safe_write_csv()`,
  which blanks any row with a count below 10.
- FPC scores, latent-class memberships and posterior draws stay in
  `dir_derived`.

## How NBDCtools is used (01_data.R)

| Step | NBDCtools |
|---|---|
| Release check | `get_latest_release_abcd()` must return `cfg$release` (7.0) |
| Variable validation | `get_dd_abcd(release = "7.0", vars = ...)`: stops if an essential variable is absent, otherwise drops and logs |
| Session table | `get_sessions_abcd(release = "7.0")` written to `output/00_sessions_table.csv` |
| Load | `create_dataset_abcd(vars, categ_to_factor, add_labels, value_to_na = TRUE)`, deliberately with no `release=` (the NBDCtools 1.1.0 bug documented in music-ed `01_data.R`) |
| Nonresponse | `value_to_na` keeps each code as a tagged NA. `haven::na_tag()` tabulates them (`00_nonresponse_codes.csv`) so refusals and don't-knows are kept apart from attrition |
| Panels | `filter_events_abcd(c("annual","even"))` for the imaging waves, `c("annual","mid_year"), connect = "or"` for binge ascertainment, `"annual"` for the annual panel |
| Cleaning | `filter_empty_rows()` |

The 7.0 static fields (family, sex, raked weight) come from the Static session
(`session_id = NA`) and are carried to every visit before any event filtering.
Site, caregiver education and income are fixed at baseline, so families stay
nested in one site.

## Variables (all verified against the 7.0 dictionary)

- **Markers:** mean cortical thickness, total cortical area and cortical volume,
  subcortical gray volume, and hippocampus and amygdala (left + right). These
  come from imaging waves 00A/02A/04A/06A with `mr_y_qc__incl__smri__t1_indicator == 1`.
  The models include ICV for volumes and area, and a crossed random effect for
  scanner serial.
- **Binge onset:** the first annual or mid-year assessment with
  `su_y_tlfb__alc__binge_ud > 0`. The onset interval runs from the last negative
  to the first positive assessment. Participants who binge at baseline are
  excluded. Time is years since baseline, from the actual visit age.
- **Covariates:** sex, baseline age, caregiver education, income, pubertal
  stage (youth PDS category and mean score), CBCL externalizing, and the raked
  propensity weight `ab_g_stc__cohort_acsps`.
- **Regional and family fields:** every Desikan thickness ROI listed by
  `get_dd_abcd(tables = "mr_y_smri__thk__dsk")`, and the birth ID
  `ab_g_stc__design_id__birth` that identifies twins.
- **Calendar time:** visit date `ab_g_dyn__visit_dtt` and visit type
  `ab_g_dyn__visit_type`.

## Modules

| Module | Methods |
|---|---|
| 02_describe_missing | Imaging patterns, monotone vs intermittent, family-clustered GEE for missing the next wave |
| 10_change | Paired difference scores; residualized change with binge onset before the last scan |
| 11_lmm | Random intercept and slope, time-invariant covariates, within/between puberty, cohort vs ageing (baseline age × time). All nested site/family/participant + scanner |
| 12_glmm_gee | Binge GLMM (glmmTMB), hurdle NB for binge days, family-clustered GEE, raked-weight GEE for population prevalence |
| 13_sem_growth | LGCM, LGCM with covariates, multigroup by sex, dual-change LCSM, bivariate LCSM, RI-CLPM (marker ↔ externalizing); FIML, cluster-robust by family |
| 20_gamm_age | gamm4 nonlinear age curves and growth rates by sex, cohort check s(age_bl), eventual bingers vs never-bingers |
| 21_fpca | PACE FPCA per marker, Happ–Greven multivariate FPCA, Cox on pre-onset FPC scores |
| 22_multivariate | Stacked multivariate LMM (correlated intercepts and slopes across markers); parallel-process LGCM |
| 23_mixture | Latent class growth mixtures 1..K (BIC, entropy, class sizes), onset hazard by class |
| 24_joint_model | **Flagship:** JMbayes2 joint model of marker trajectory and binge onset (value, value + slope, multivariate, midpoint-onset sensitivity), compared with a naive time-varying Cox |
| 25_onset_anchored | Piecewise LMM: pre-onset differences of eventual bingers, then level shift and slope change after onset; sensitivity analyses |
| 26_missing_data | MAR (likelihood, IPW, MI) vs MNAR (delta-adjusted MI tipping point, pattern-mixture CCMV/NCMV) |
| 28_temporal_trends | Age-period-cohort: calendar-year smooth, years of pandemic exposure at each scan, birth cohort, COVID era and remote visits in the binge GLMM, Cox model and joint model with a time-varying pandemic era |
| 30_target_trial | Target-trial emulation of binge initiation (nested trials at annual visits): IPTW/IPCW marginal structural model and parametric g-formula with family bootstrap |
| 31_discrete_time_onset | Discrete-time (cloglog, offset log interval length) hazard of onset, with each interval's calendar exposure to the pandemic and the last scan before the interval |
| 32_within_family | Sibling and co-twin designs (within-between/Mundlak) for the onset-anchored model; discordant-family counts |
| 33_harmonization | Longitudinal ComBat (longCombat) vs scanner random effect |
| 34_regional | All Desikan thickness ROIs: piecewise onset model with BH-FDR; PCA composite of ROIs as a joint-model marker |
| 35_dynamic_prediction | Joint-model individualized risk: time-dependent AUC and Brier score, family-level K-fold CV, vs covariates-only Cox |
| 36_puberty_time | Pubertal clock (gamm4 age vs PDS vs both, AIC), pubertal timing and brain slope, Cox with timing and time-varying PDS |
| 37_evalues | E-values for onset effects (25, 30, 32) and hazard ratios (24, 31) |
| 27_synthesis | `99_key_estimates.csv`: headline estimates in one table (runs last) |

## Tests

`tests/mock_nbdc.R` simulates data in the shape `create_dataset_abcd()` returns:

- factors for categorical variables
- a Static session
- annual, mid-year and screener sessions
- scanner upgrades, QC failures and MAR attrition

It also builds in known effects. The run sources the real `filter_events_abcd()`
and `filter_empty_rows()`, from the installed package or from an NBDCtools
source checkout (`NBDCTOOLS_SRC`).

## Packages

NBDCtools, dplyr, tidyr, purrr, lme4, lmerTest, glmmTMB, geepack, lavaan, mgcv, gamm4,
fdapace, lcmm, JMbayes2, nlme, survival, mice, longCombat (GitHub: jcbeer/longCombat), broom,
broom.mixed and haven.

### Package problems

`run_all.R` checks packages before any module runs. It installs missing packages, updates a dependency that is too old, and probes lme4/lmerTest by fitting and converting a tiny model. If something was installed or updated, or a package is loaded in the session at an older version than the one on disk, the run stops with a "restart R" message. Restart R (Session > Restart R in RStudio) and source `run_all.R` again, or run it in a fresh process:

    Rscript -e 'source("run_all.R")'

If lmerTest still cannot be used, mixed models fall back to the lme4 fit with normal-approximation p-values and the log says so.

## Data and disclosure control

The repository contains no ABCD data. ABCD release 7.0 is available to qualified researchers through the NIH Brain Development Cohorts Data Hub. Set `ABCD_DATA_PATH` (the tabulated release) and `LP_DERIVED_PATH` (a folder **outside** the repository for participant-level panels and fitted models) before running; the run stops if `dir_derived` is inside the repository. Every table is written through `safe_write_csv()`, which blanks any row with a count below 10.

## Related resources

- Tutorials: [longitudinal.dev](https://longitudinal.dev)
- Hawes et al. (2025). Longitudinal analysis of the ABCD study. *Developmental Cognitive Neuroscience* 72:101518. [doi:10.1016/j.dcn.2025.101518](https://doi.org/10.1016/j.dcn.2025.101518)
- Li et al. (2025). Missing data approaches for longitudinal neuroimaging research. *Developmental Cognitive Neuroscience* 74:101563. [doi:10.1016/j.dcn.2025.101563](https://doi.org/10.1016/j.dcn.2025.101563)
- Causal inference pipeline for ABCD: [github.com/wes-libr/abcd-causal-inference](https://github.com/wes-libr/abcd-causal-inference)

## License

MIT; see `LICENSE`.
