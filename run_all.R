# =============================================================================
# run_all.R : ABCD 7.0 longitudinal pipeline
#   setwd("<abcd-analyses>/longitudinal-pipeline/pipeline"); source("run_all.R")
# =============================================================================
rm(list = setdiff(ls(), c(".mock", ".lp_args", grep("^(get_|create_dataset|filter_|simulate_)", ls(), value = TRUE))))
source("config.R")
if (exists(".lp_args")) cfg <- modifyList(cfg, .lp_args)
source("R/00_utils.R")
lp_setup(cfg)
lp_log("=== run start: ", format(Sys.time()), "; modules: ", paste(cfg$modules, collapse = ", "))
pkg <- lp_check_packages()
if (length(pkg$restart) && isTRUE(cfg$stop_if_missing))
  stop("Packages were updated or an outdated copy is loaded (", paste(pkg$restart, collapse = "; "),
       "). Restart R (RStudio: Session > Restart R; R.app: quit and reopen), then source(\"run_all.R\") again; or run from a terminal: ",
       "Rscript -e 'source(\"run_all.R\")'", call. = FALSE)
if (length(pkg$missing) && isTRUE(cfg$stop_if_missing))
  stop("Missing packages: ", paste(pkg$missing, collapse = ", "), " -- see the load error(s) logged above; to override settings use .lp_args <- list(...) before source()", call. = FALSE)
if (!isTRUE(pkg$lmertest_ok))
  lp_log("WARNING: lmerTest unusable; mixed models use lme4 fits with normal-approximation p-values")

source("R/01_data.R")
status <- tibble(module = character(), ok = logical(), seconds = numeric(), error = character())
for (m in cfg$modules) {
  t0 <- Sys.time()
  err <- tryCatch({ source(file.path("R", paste0(m, ".R")), local = new.env()); NA_character_ },
                  error = function(e) conditionMessage(e))
  status <- add_row(status, module = m, ok = is.na(err),
                    seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), error = err)
  lp_log(m, if (is.na(err)) " OK" else paste(" FAILED:", err))
}
safe_write_csv(status, "99_run_status")
lp_log("=== run end; ", sum(status$ok), "/", nrow(status), " modules OK")
