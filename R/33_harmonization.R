# =============================================================================
# 33_harmonization.R : longitudinal ComBat vs scanner random effect
# =============================================================================
# Sites upgraded scanners during follow-up, so scanner is partly confounded
# with time. Two ways to handle it are compared for every marker:
#   (a) raw data, scanner as a crossed random intercept (the pipeline default)
#   (b) longitudinal ComBat (Beer et al. 2020; longCombat): empirical-Bayes
#       removal of scanner-specific additive and multiplicative effects while
#       preserving time, sex and baseline age effects, with a participant
#       random intercept; then the same LMM without a scanner term.
# Agreement of the slopes supports either choice; disagreement flags scanner-
# time confounding to investigate before the substantive models are read.
# Harmonized values are saved to dir_derived for use by other modules.
lp_log("33_harmonization")
if (!need(c("lme4", "longCombat"))) stop("lme4 and longCombat required (remotes::install_github('jcbeer/longCombat'))")
d <- img_data() |> mutate(id_n = as.integer(factor(id)))
mk <- names(cfg$brain)
dc <- d |> select(id, id_n, fam, site, time, female, age_bl_c, icv_z, scanner, all_of(mk)) |> filter(if_all(all_of(mk), ~ !is.na(.x)))
small <- dc |> count(scanner) |> filter(n < 2) |> pull(scanner)
dc <- dc |> filter(!scanner %in% small) |> as.data.frame()
lc <- try_fit("longCombat", longCombat::longCombat(idvar = "id_n", timevar = "time", batchvar = "scanner", features = mk,
                                                  formula = "time + female + age_bl_c", ranef = "(1|id_n)", data = dc, verbose = FALSE))
if (is.null(lc)) stop("longCombat failed")
hz <- bind_cols(dc, lc$data_combat[, paste0(mk, ".combat")])
save_derived(hz, "harmonized_longcombat")
res <- list()
for (m in mk) {
  cov <- paste0(" + age_bl_c + female", icv_term(m))
  fa <- try_fit(paste(m, "raw + scanner RE"), lme4::lmer(as.formula(paste0(m, " ~ time", cov, " + (1 | site) + (1 | site:fam) + (1 + time | site:fam:id) + (1 | scanner)")),
                                                         data = hz, control = lme4::lmerControl(calc.derivs = FALSE)))
  fb <- try_fit(paste(m, "longCombat"), lme4::lmer(as.formula(paste0("`", m, ".combat` ~ time", cov, " + (1 | site) + (1 | site:fam) + (1 + time | site:fam:id)")),
                                                   data = hz, control = lme4::lmerControl(calc.derivs = FALSE)))
  for (k in list(list("scanner_random_effect", fa), list("longCombat", fb))) if (!is.null(k[[2]])) {
    b <- lme4::fixef(k[[2]]); se <- sqrt(diag(as.matrix(vcov(k[[2]]))))
    res[[paste(m, k[[1]])]] <- tibble(module = "33", model = k[[1]], outcome = m, term = "time (per year)",
                                      estimate = b[["time"]], std.error = se[["time"]],
                                      conf.low = b[["time"]] - 1.96 * se[["time"]], conf.high = b[["time"]] + 1.96 * se[["time"]],
                                      residual_sd = sigma(k[[2]]), n_id = n_distinct(hz$id), n_scanners = n_distinct(hz$scanner))
  }
}
out <- bind_rows(res)
safe_write_csv(out, "33_harmonization", count_cols = c("n_id", "n_scanners"))
if (nrow(out)) key_write(out |> filter(outcome == cfg$primary) |> select(model, outcome, term, estimate, conf.low, conf.high, n_id), "33")
