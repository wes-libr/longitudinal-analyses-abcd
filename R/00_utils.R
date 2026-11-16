# =============================================================================
# 00_utils.R : shared helpers (packages, blinding, disclosure, model helpers)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
})

## ---- Packages ---------------------------------------------------------------
# NBDCtools is required for real data. Tests set options(lp.mock = TRUE) and
# provide NBDCtools-compatible stand-ins (tests/mock_nbdc.R).
if (!isTRUE(getOption("lp.mock"))) {
  if (!requireNamespace("NBDCtools", quietly = TRUE))
    stop("Install NBDCtools (and NBDCtoolsData >= 5.0.0 for ABCD 7.0).", call. = FALSE)
  library(NBDCtools)
}

need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    lp_log("  missing package(s): ", paste(miss, collapse = ", "), " -> module skipped")
    return(FALSE)
  }
  TRUE
}

## ---- Package check / install (run before any module) -----------------------
# Which packages each module needs. lp_check_packages() installs missing CRAN
# packages (and longCombat from GitHub) when cfg$auto_install is TRUE, then
# reports any module that would still be skipped, so nothing drops silently.
lp_module_packages <- list(
  "02_describe_missing" = "geepack", "10_change" = "geepack", "11_lmm" = c("lme4", "lmerTest"),
  "12_glmm_gee" = c("glmmTMB", "geepack"), "13_sem_growth" = "lavaan", "20_gamm_age" = c("gamm4", "mgcv"),
  "21_fpca" = c("fdapace", "survival"), "22_multivariate" = c("glmmTMB", "lavaan"), "23_mixture" = c("lcmm", "survival"),
  "24_joint_model" = c("JMbayes2", "nlme", "survival"), "25_onset_anchored" = "lme4",
  "26_missing_data" = c("lme4", "geepack", "mice"), "28_temporal_trends" = c("lme4", "glmmTMB", "survival", "gamm4", "JMbayes2", "nlme"),
  "30_target_trial" = "geepack", "31_discrete_time_onset" = "geepack", "32_within_family" = "lme4",
  "33_harmonization" = c("lme4", "longCombat"), "34_regional" = c("lme4", "JMbayes2", "nlme", "survival"),
  "35_dynamic_prediction" = c("JMbayes2", "nlme", "survival"), "36_puberty_time" = c("lme4", "survival", "mgcv", "gamm4"),
  "37_evalues" = character(0), "27_synthesis" = character(0))
lp_base_packages <- c("dplyr", "tidyr", "purrr", "tibble", "broom", "broom.mixed", "haven")
lp_github_packages <- c(longCombat = "jcbeer/longCombat")

# Fit a tiny lme4 model and convert it with lmerTest: catches lme4/lmerTest
# version mismatches ("Unable to extract deviance function from model fit")
# before any module runs. Deterministic data; does not touch the RNG.
lp_probe_lmertest <- function() {
  if (nzchar(Sys.getenv("LP_TEST_FORCE_LMERTEST_FAIL"))) return("forced failure (test)")
  if (!requireNamespace("lme4", quietly = TRUE) || !requireNamespace("lmerTest", quietly = TRUE))
    return("lme4 or lmerTest not installed")
  d <- data.frame(g = factor(rep(1:20, each = 4)), x = rep(0:3, 20))
  d$y <- d$x + sin(as.integer(d$g)) + cos(seq_len(nrow(d)) * 1.7)
  tryCatch({
    f <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d))
    m <- lmerTest::as_lmerModLmerTest(f)
    if (!"Pr(>|t|)" %in% colnames(summary(m)$coefficients)) "no p-values from lmerTest" else NA_character_
  }, error = function(e) conditionMessage(e))
}

lp_check_packages <- function(modules = cfg$modules, install = isTRUE(cfg$auto_install)) {
  repos <- "https://cloud.r-project.org"
  mock <- isTRUE(getOption("lp.mock"))
  want <- unique(c(lp_base_packages, unlist(lp_module_packages[modules], use.names = FALSE)))
  have <- function(p) requireNamespace(p, quietly = TRUE)
  restart <- character(0)

  ## 1. missing packages
  miss <- want[!vapply(want, have, logical(1))]
  if (length(miss) && install && !mock) {
    cran <- setdiff(miss, names(lp_github_packages))
    if (length(cran)) {
      lp_log("installing from CRAN: ", paste(cran, collapse = ", "))
      utils::install.packages(cran, repos = repos)
    }
    gh <- intersect(miss, names(lp_github_packages))
    if (length(gh)) {
      if (!have("remotes")) utils::install.packages("remotes", repos = repos)
      for (g in gh) { lp_log("installing from GitHub: ", lp_github_packages[[g]]); remotes::install_github(lp_github_packages[[g]], upgrade = "never") }
    }
    miss <- want[!vapply(want, have, logical(1))]
  }

  ## 2. packages that still fail to load: log why; if a dependency is too old
  ##    ("namespace 'X' a is already loaded, but >= b is required" or
  ##    "namespace 'X' a is being loaded, but >= b is required"), install X
  for (p in miss) {
    why <- tryCatch({ loadNamespace(p); "loads now (restart R?)" }, error = function(e) conditionMessage(e))
    lp_log("  ", p, " load error: ", why)
    lks <- list.files(.libPaths()[1], pattern = paste0("^00LOCK"), full.names = TRUE)
    if (length(lks)) lp_log("  stale install lock(s) in ", .libPaths()[1], ": ", paste(basename(lks), collapse = ", "))
    dep <- if (grepl("namespace", why) && grepl("is required", why))
      sub("^.*namespace[^A-Za-z0-9.]+([A-Za-z][A-Za-z0-9.]*)[^A-Za-z0-9.].*$", "\\1", why) else character(0)
    if (length(dep) && dep != why) {
      if (install && !mock) { lp_log("  updating dependency ", dep, " required by ", p); utils::install.packages(dep, repos = repos) }
      restart <- c(restart, sprintf("%s needs a newer %s", p, dep))
    }
  }

  ## 3. loaded-but-outdated namespaces: an update was installed while an older
  ##    copy is loaded in this R session; only a restart picks up the new one
  for (p in intersect(c(want, "reformulas", "Matrix", "TMB", "Rcpp", "nlme"), loadedNamespaces())) {
    inst <- tryCatch(utils::packageVersion(p), error = function(e) NULL)
    if (!is.null(inst) && inst != package_version(getNamespaceVersion(p)))
      restart <- c(restart, sprintf("%s loaded %s but %s installed", p, getNamespaceVersion(p), as.character(inst)))
  }

  ## 4. lme4/lmerTest compatibility probe (modules 11, 25, 28, 32, 34, 36)
  lmt_needed <- any(vapply(lp_module_packages[modules], function(x) "lme4" %in% x, logical(1)))
  lmt <- if (lmt_needed && !any(c("lme4", "lmerTest") %in% miss)) lp_probe_lmertest() else NA_character_
  if (!is.na(lmt)) {
    lp_log("  lmerTest probe failed: ", lmt)
    ## reinstall at most once per lme4/lmerTest version pair, so a pair that is
    ## genuinely incompatible does not cause a reinstall-and-restart loop
    vers <- paste0("lme4 ", utils::packageVersion("lme4"), ", lmerTest ", utils::packageVersion("lmerTest"))
    mark <- file.path(cfg$dir_output, ".lmertest_reinstalled")
    tried <- file.exists(mark) && identical(readLines(mark, warn = FALSE)[1], vers)
    if (install && !mock && !tried) {
      lp_log("  reinstalling lmerTest (and lme4) to match")
      utils::install.packages(c("lme4", "lmerTest"), repos = repos)
      vers_new <- paste0("lme4 ", utils::packageDescription("lme4", fields = "Version"),
                         ", lmerTest ", utils::packageDescription("lmerTest", fields = "Version"))
      writeLines(vers_new, mark)
      restart <- c(restart, "lme4/lmerTest reinstalled")
    } else if (tried) lp_log("  already reinstalled for ", vers, "; not retrying")
  }
  options(lp.lmertest_ok = is.na(lmt))

  blocked <- names(Filter(function(p) any(p %in% miss), lp_module_packages[modules]))
  if (length(miss)) lp_log("WARNING: still missing ", paste(miss, collapse = ", "),
                           " -> affected modules: ", paste(blocked, collapse = ", "))
  if (length(restart)) lp_log("RESTART R NEEDED: ", paste(unique(restart), collapse = "; "))
  if (!length(miss) && !length(restart) && is.na(lmt)) lp_log("all packages for the selected modules are installed and load")
  invisible(list(missing = miss, affected_modules = blocked, restart = unique(restart), lmertest_ok = is.na(lmt)))
}

## ---- Logging ----------------------------------------------------------------
lp_log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  message(msg)
  if (!is.null(.lp_env$logfile)) cat(msg, "\n", file = .lp_env$logfile, append = TRUE)
}
.lp_env <- new.env()

## ---- Paths and blinding guard -----------------------------------------------
lp_setup <- function(cfg, project_dir = getwd()) {
  dir.create(cfg$dir_output, showWarnings = FALSE, recursive = TRUE)
  dir.create(cfg$dir_derived, showWarnings = FALSE, recursive = TRUE)
  der <- normalizePath(cfg$dir_derived, winslash = "/", mustWork = FALSE)
  # Walk up to the shared repo root (the folder containing analyses/ or
  # longitudinal-pipeline/) and refuse to cache inside it.
  root <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:5) {
    if (dir.exists(file.path(root, "analyses")) ||
        basename(root) %in% c("abcd-analyses")) break
    root <- dirname(root)
  }
  if (startsWith(paste0(der, "/"), paste0(normalizePath(project_dir, "/", FALSE), "/")) ||
      startsWith(paste0(der, "/"), paste0(root, "/")))
    stop("BLINDING VIOLATION: cfg$dir_derived (", der, ") is inside the shared repo (",
         root, "). Point it outside abcd-analyses/.", call. = FALSE)
  .lp_env$logfile <- file.path(cfg$dir_output, "run_log.txt")
  set.seed(cfg$seed)
  invisible(TRUE)
}

save_derived <- function(obj, name) saveRDS(obj, file.path(cfg$dir_derived, paste0(name, ".rds")))
load_derived <- function(name) {
  f <- file.path(cfg$dir_derived, paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing derived object ", name, ": run 01_data first.", call. = FALSE)
  readRDS(f)
}

## ---- Disclosure control (same rule as music-ed: no cell < 10) ---------------
safe_write_csv <- function(df, name, count_cols = character(0),
                           key_cols = c("module", "model", "outcome", "marker",
                                        "term", "wave", "level", "variable"),
                           min_cell = cfg$min_cell) {
  df <- as.data.frame(df)
  path <- file.path(cfg$dir_output, paste0(name, ".csv"))
  count_cols <- intersect(count_cols, names(df))
  value_cols <- setdiff(names(df), c(count_cols, intersect(key_cols, names(df))))
  small <- rep(FALSE, nrow(df))
  for (cc in count_cols) {
    x <- suppressWarnings(as.numeric(df[[cc]]))
    small <- small | (!is.na(x) & x < min_cell)
  }
  if (any(small)) for (cc in c(count_cols, value_cols)) {
    if (is.factor(df[[cc]])) df[[cc]] <- as.character(df[[cc]])
    df[[cc]][small] <- if (is.numeric(df[[cc]])) NA else "<10"
  }
  utils::write.csv(df, path, row.names = FALSE)
  lp_log("  wrote ", basename(path), if (any(small)) sprintf(" [%d rows suppressed]", sum(small)) else "")
  invisible(df)
}

## ---- Small helpers ----------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Years since baseline, parsed from an NBDC session label (ses-04A -> 4).
years_from_session <- function(s) num(sub("^ses-(\\d+).*$", "\\1", as.character(s)))

# Nested random-effects terms used throughout (site / family / participant),
# optionally crossed with scanner. `slope` adds a participant random slope.
re_terms <- function(slope = NULL, scanner = cfg$scanner_re) {
  p <- if (is.null(slope)) "(1 | site:fam:id)" else paste0("(1 + ", slope, " | site:fam:id)")
  paste(c("(1 | site)", "(1 | site:fam)", p, if (scanner) "(1 | scanner)"), collapse = " + ")
}

# Standardize a marker using baseline-wave mean/SD (keeps trajectories on one
# scale across markers; FIXED constants, not recomputed per subset).
std_by_baseline <- function(x, is_bl) (x - mean(x[is_bl], na.rm = TRUE)) / sd(x[is_bl], na.rm = TRUE)

# Tidy a fitted model into a standard aggregate table.
tidy_fixed <- function(fit, module, model, outcome) {
  out <- tryCatch({
    if (inherits(fit, c("merMod", "glmmTMB", "lmerModLmerTest")) && requireNamespace("broom.mixed", quietly = TRUE))
      broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
    else if (inherits(fit, "geeglm")) {
      s <- summary(fit)$coefficients
      tibble(term = rownames(s), estimate = s$Estimate, std.error = s$Std.err,
             conf.low = s$Estimate - 1.96 * s$Std.err, conf.high = s$Estimate + 1.96 * s$Std.err,
             p.value = s[["Pr(>|W|)"]])
    } else broom::tidy(fit, conf.int = TRUE)
  }, error = function(e) tibble(term = NA_character_, estimate = NA_real_))
  if (!"p.value" %in% names(out) && all(c("estimate", "std.error") %in% names(out)))
    out$p.value <- 2 * pnorm(-abs(out$estimate / out$std.error))
  if (all(c("estimate", "std.error") %in% names(out))) {
    if (!"conf.low" %in% names(out)) out$conf.low <- out$estimate - 1.96 * out$std.error
    if (!"conf.high" %in% names(out)) out$conf.high <- out$estimate + 1.96 * out$std.error
  }
  for (v in c("std.error", "conf.low", "conf.high", "p.value")) if (!v %in% names(out)) out[[v]] <- NA_real_
  out <- out[, intersect(c("component", "term", "estimate", "std.error", "conf.low",
                           "conf.high", "p.value"), names(out))]
  tibble(module = module, model = model, outcome = outcome, out)
}

# Variance components from any lme4/glmmTMB fit.
tidy_vc <- function(fit, module, model, outcome) {
  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  if (is.null(vc) && inherits(fit, "glmmTMB")) {
    v <- glmmTMB::VarCorr(fit)$cond
    vc <- bind_rows(lapply(names(v), function(g) tibble(grp = g, var1 = rownames(v[[g]]), var2 = NA, vcov = diag(v[[g]]), sdcor = sqrt(diag(v[[g]])))))
  }
  if (is.null(vc)) return(NULL)
  tibble(module = module, model = model, outcome = outcome, as_tibble(vc))
}

# Run one model, log and swallow failures so run_all continues.
try_fit <- function(label, expr) {
  t0 <- Sys.time()
  res <- tryCatch(expr, error = function(e) { lp_log("  FAILED ", label, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) lp_log("  fit ", label, sprintf(" (%.1fs)", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res
}

# Imaging panel with binge-onset information attached.
#   eventual = participant has a binge onset during follow-up
#   post     = scan acquired after the onset assessment
#   tso      = years since onset (0 before onset / never-bingers)
img_data <- function() {
  ev <- load_derived("events")
  load_derived("panel_img") |>
    left_join(ev |> select(id, event, t_event), by = "id") |>
    mutate(eventual = as.integer(!is.na(event) & event == 1),
           post     = as.integer(eventual == 1 & time > t_event),
           tso      = ifelse(post == 1, time - t_event, 0),
           age_bl_c = age_bl - 10,
           female_f = factor(female),
           fam      = paste(site, fam, sep = "_"))   # family nested in baseline site
}
icv_term <- function(marker) if (sub("_z$", "", marker) %in% cfg$icv_adjust) " + icv_z" else ""

# Long panel -> wide by imaging wave for SEM modules (columns <var>_<wave_yr>).
img_wide <- function(d, vars, waves = c(0, 2, 4, 6)) {
  d |> filter(wave_yr %in% waves) |>
    select(id, fam, female, age_bl, edu_bl, wave_yr, all_of(vars)) |>
    pivot_wider(names_from = wave_yr, values_from = all_of(vars), names_sep = "_")
}

# Headline estimates collected by 27_synthesis (one small file per module).
key_write <- function(df, module) safe_write_csv(df |> mutate(module = module), paste0("zz_key_", module))

# Count participants/observations by group for descriptive tables.
n_table <- function(d, ...) d |> group_by(...) |> summarise(n_obs = n(), n_id = n_distinct(id), .groups = "drop")
